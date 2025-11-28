nextflow.enable.dsl = 2

params.assemblies = "$projectDir/data/tutorial_assemblies/*_trinity.Trinity.fasta"
params.vsg_db = "$projectDir/data/blastdb/concatAnTattb427.fa"
params.notvsg_db = "$projectDir/data/blastdb/vsgseq2NOTvsgs.fa"
params.vsgome = "$projectDir/data/blastdb/concatAnTattb427.fa"
params.genome = "$projectDir/data/genome/TriTrypDB-67_TbruceiEATRO1125_Genome.fasta"
params.decoys = "$projectDir/data/genome/TriTrypDB-67_TbruceiEATRO1125_Genome_decoys.txt"
params.full_vsg_db = ""
params.cores = "4"
params.trinitymem = "20"
params.cdslength = "300"
params.partial = false
params.cdhit_id = "0.94"
params.cdhit_as = "0.94"
params.threshold = "100000"
params.today = new Date().format('ddMMMYY')
params.outdir = "results/${params.today}"
params.samplesheet = "$projectDir/data/reads/samples.csv"
params.help = ""
params.mode = "full"
params.conda_yml = "$projectDir/vsgseq2.yml"
params.scripts = "$projectDir/bin/"

def validateParams() {
    def errors = []

    try {
        if (file(params.assemblies).isEmpty()) {
            errors << "The assemblies path '${params.assemblies}' does not exist or is empty."
        }
        
        def vsg_db = file(params.vsg_db)
        if (!vsg_db.exists()) {
            errors << "The VSG database path '${params.vsg_db}' does not exist."
        }
        
        def notvsg_db = file(params.notvsg_db)
        if (!notvsg_db.exists()) {
            errors << "The NOTVSG database path '${params.notvsg_db}' does not exist."
        }
        
        def vsgome = file(params.vsgome)
        if (!vsgome.exists()) {
            errors << "The VSGome path '${params.vsgome}' does not exist."
        }
        
        if (params.full_vsg_db) {
            def full_vsg_db = file(params.full_vsg_db)
            if (!full_vsg_db.exists()) {
                errors << "The full VSG database path '${params.full_vsg_db}' does not exist."
            }
        }
        
        def samplesheet = file(params.samplesheet)
        if (!samplesheet.exists()) {
            errors << "The samplesheet path '${params.samplesheet}' does not exist."
        }
        
        if (!["full", "analyse"].contains(params.mode)) {
            errors << "Invalid mode '${params.mode}'. Allowed values are: full and analyse."
        }
        if (params.cores <= 0) {
            errors << "The cores '${params.cores}' must be greater than 0."
        }
        if (params.trinitymem.toInteger() <= 0) {
            errors << "The Trinity memory '${params.trinitymem}' must be greater than 0."
        }
        if (params.cdslength.toInteger() <= 0) {
        errors << "The CDS length '${params.cdslength}' must be greater than 0."
        }
        if (params.cdhit_id.toFloat() < 0.0 || params.cdhit_id.toFloat() > 1.0) {
            errors << "The CD-HIT identity threshold '${params.cdhit_id}' must be between 0.0 and 1.0."
        }
        if (params.cdhit_as.toFloat() < 0.0 || params.cdhit_as.toFloat() > 1.0) {
            errors << "The CD-HIT alignment coverage '${params.cdhit_as}' must be between 0.0 and 1.0."
        }

    } catch (Exception e) {
        errors << "Error validating parameters: ${e.message}"
    }

    if (errors) {
        log.info("Parameter validation failed with the following errors:")
        errors.each { error -> 
            System.err.println("ERROR: ${error}")
        }
        System.exit(1)
    }
}

include { TRIM } from './modules/trim'
include { ASSEMBLE } from './modules/assemble'
include { ORF } from './modules/orf'
include { BLAST } from './modules/blast'
include { CONCATENATE_VSGS } from './modules/concatenate'
include { CONCATENATED_CDHIT } from './modules/concatenated_cdhit'
include { INDEX } from './modules/index'
include { QUANTIFY } from './modules/quantify'
include { MULTIQC } from './modules/multiqc'
include { SUMMARISE } from './modules/summarise'

workflow help {
if (params.help) {
    help = """VSGSEQ2.nf: A pipeline for analysing VSGSeq data
             |
             |Required arguments:
             |
             |  --assemblies Location of assemblies
             |                [default: ${params.assemblies}]
             |  --vsg_db    Location of VSGdb
             |                [default: ${params.vsg_db}]
             |  --notvsg_db Location of NOTVSGdb
             |                [default: ${params.notvsg_db}]
             |  --vsgome    Location of VSGome
             |                [default: ${params.vsgome}]
             |  --full_vsg_db    Location of a database to add into the VSGome (such as data/blastdb/concatAnTattb427_full.fa). 
             |                    Default will only include the assembled VSGome.
             |  --mode    The mode to run the pipeline in. Options are full, analyse.
             |                [default: ${params.mode}]
             |  --outdir        VSGSeq output directory. 
             |                [default: ${params.outdir}]
             |  --samplesheet  Define the path to the samplesheet.
             |                [default: ${params.samplesheet}]
             |                
             |
             |Optional arguments:
             |
             |  --requestedcpus  Define number of cores VSGSeq2 will use.
             |                [default: ${params.requestedcpus}]
             |  --cores  Define number of cores Trinity and other tools will use.
             |                [default: ${params.cores}]
             |  --trinitymem  Define memory amount for Trinity in Gb.
             |                [default: ${params.trinitymem} Gb]
             |  --cdslength    Define minimum CDS length (amino acids).
             |                [default: ${params.cdslength}]
             |  --partial       Only complete ORFs are produced as standard. It is possible to allow the production of partial ORFs by adding the --partial flag.
             |                [default: ${params.partial}]
             |  --cdhit_id       Define sequence identity threshold - how much the alignment has to match (0.0 - 1.0).
             |                [default: ${params.cdhit_id}]
             |  --cdhit_as       Define alignment coverage for the shorter sequence (0.0 - 1.0).
             |                [default: ${params.cdhit_as}]
             |  --help         Print this message.""".stripMargin()
    println(help)
    exit(0)
}
}

workflow vsgseq2 {
    main:

    validateParams()

    ch_samplesheet = Channel.fromPath(params.samplesheet)

    ch_reads = ch_samplesheet
        .splitCsv(header:true)
        .map { row -> 
        def meta = [:]
        meta.id = row.sample
        meta.prefix = row.sample
        meta.single_end = row.r2.toString().isEmpty()
        
        def reads = []
        def r1 = file(row.r1, checkIfExists: true)
        def r2 = meta.single_end ? null : file(row.r2, checkIfExists: true)

        reads.add(r1)
        if (r2) reads.add(r2)

        [ meta, reads ]
    }

    ch_blast_dbs = Channel.fromPath("${params.vsg_db}*")
        .mix(Channel.fromPath("${params.notvsg_db}*"))
        .collect()

    if (params.mode == "full") {
        TRIM(
            ch_reads,
            params.cores
            )
        ASSEMBLE(
            TRIM.out,
            params.cores,
            params.trinitymem
            )
        ORF(
            ASSEMBLE.out,
            params.cdslength,
            params.partial
            )
        BLAST(
            ORF.out,
            ch_blast_dbs
            )
        CONCATENATE_VSGS(
            (BLAST.out.vsgs).collect(),
            params.full_vsg_db
            )
        CONCATENATED_CDHIT(
            CONCATENATE_VSGS.out,
            params.cdhit_id,
            params.cdhit_as
            )
        INDEX(
            CONCATENATE_VSGS.out,
            params.cores,
            params.genome,
            params.decoys
            )
        QUANTIFY(
            INDEX.out,
            params.cores,
            TRIM.out
            )
        MULTIQC(
            (QUANTIFY.out.quants).collect())
        SUMMARISE(
            (QUANTIFY.out.quants).collect(),
            params.threshold,
            (BLAST.out.vsgs).collect(),
            CONCATENATED_CDHIT.out.clstr,
            CONCATENATE_VSGS.out)
        
    }
    else if (params.mode == "analyse"){
        TRIM(
            ch_reads,
            params.cores
            )
        assemblies_ch = 
        Channel.fromPath(
            params.assemblies,
            checkIfExists: true
            )
        ORF(
            assemblies_ch,
            params.cdslength,
            params.partial
            )
        BLAST(
            ORF.out,
            ch_blast_dbs
            )
        CONCATENATE_VSGS(
            (BLAST.out.vsgs).collect(),
            params.full_vsg_db
            )
        CONCATENATED_CDHIT(
            CONCATENATE_VSGS.out,
            params.cdhit_id,
            params.cdhit_as
            )
        INDEX(
            CONCATENATE_VSGS.out,
            params.cores,
            params.genome,
            params.decoys
            )
        QUANTIFY(
            INDEX.out,
            params.cores,
            TRIM.out)
        MULTIQC(
            (QUANTIFY.out.quants).collect()
            )
        SUMMARISE(
            (QUANTIFY.out.quants).collect(),
            params.threshold,
            (BLAST.out.vsgs).collect(),
            CONCATENATED_CDHIT.out.clstr,
            CONCATENATE_VSGS.out
            )
    }  
    else {
        log.error("Invalid mode selected. Please select one of the following: full, assemble, predictvsgs, quantify, analyse.")
        System.exit(1)
    }

    emit:
        summary_champions = SUMMARISE.out.champ_vsgs
}

workflow {
    help ()
    vsgseq2()

    workflow.onComplete {
        log.info("\nWorkflow completed!")
    }

    workflow.onError {
        log.error("Workflow execution failed: ${workflow.errorMessage}")
    }
}