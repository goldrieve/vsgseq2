nextflow.enable.dsl = 2

params.assemblies = "$projectDir/data/tutorial_assemblies/*_trinity.Trinity.fasta"
params.vsg_db = "$projectDir/data/blastdb/concatAnTattb427.fa"
params.notvsg_db = "$projectDir/data/blastdb/vsgseq2NOTvsgs.fa"
params.vsgome = "$projectDir/data/blastdb/concatAnTattb427.fa"
params.full_vsg_db = ""
params.requestedcpus = 4
params.cores = "4"
params.trinitymem = "20"
params.cdslength = "300"
params.cdhit_id = "0.94"
params.cdhit_as = "0.94"
params.threshold = "100000"
def timestamp = new java.text.SimpleDateFormat("yyyyMMdd_HHmmss").format(new Date())
params.outdir = "results/${timestamp}"
params.samplesheet = "$projectDir/data/reads/samples.csv"
params.help = ""
params.mode = "full"

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
             |  --cdhit_id       Define sequence identity threshold - how much the alignment has to match (0.0 - 1.0).
             |                [default: ${params.cdhit_id}]
             |  --cdhit_as       Define alignment coverage for the shorter sequence (0.0 - 1.0).
             |                [default: ${params.cdhit_as}]
             |  --threshold       Define the lowest number of reads a sample must have mapped to the VSGome to include in the filtered tpm.csv.
             |                [default: ${params.threshold}]
             |  --help         Print this message.""".stripMargin()
    // Print the help with the stripped margin and exit
    println(help)
    exit(0)
}

// Add a step to validate parameters
def validateParams() {
    def errors = []

    // Validate file existence using file() method
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
        
        // Validate mode
        if (!["full", "analyse"].contains(params.mode)) {
            errors << "Invalid mode '${params.mode}'. Allowed values are: full and analyse."
        }
        if (params.requestedcpus <= 0) {
            errors << "The requested CPUs '${params.requestedcpus}' must be greater than 0."
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
            // Print error message without the log file reference
            System.err.println("ERROR: ${error}")
        }
        System.exit(1)
    }
}

// Call the validation function before the workflow starts
validateParams()

include { TRIM } from './modules/trim'
include { ASSEMBLE } from './modules/assemble'
include { ORF } from './modules/orf'
include { INDIVIDUAL_CDHIT } from './modules/cdhit'
include { BLAST } from './modules/blast'
include { CONCATENATE_VSGS } from './modules/concatenate'
include { CONCATENATED_CDHIT } from './modules/concatenated_cdhit'
include { INDEX } from './modules/index'
include { QUANTIFY } from './modules/quantify'
include { MULTIQC } from './modules/multiqc'
include { SUMMARISE } from './modules/summarise'
include { SUMMARISEVSGOME } from './modules/summarise_vsgome'

ch_samplesheet = Channel.fromPath(params.samplesheet)

ch_reads = ch_samplesheet
    .splitCsv(header:true)
    .map { row -> 
        def meta = [:]
        meta.id = row.sample    // first column = prefix
        meta.prefix = row.sample
        meta.single_end = row.r2.toString().isEmpty()
        
        def reads = []
        def r1 = file(row.r1, checkIfExists: true)
        def r2 = meta.single_end ? null : file(row.r2, checkIfExists: true)

        reads.add(r1)
        if (r2) reads.add(r2)

        [ meta, reads ]
    }

workflow {
    if (params.mode == "old_full") {
        trimmed_reads_ch = TRIM(ch_reads, params.cores)
        assemble_ch = ASSEMBLE(trimmed_reads_ch, params.cores, params.trinitymem)
        orf_ch = ORF(assemble_ch, params.cdslength)
        indivcdhit_ch = INDIVIDUAL_CDHIT(orf_ch.fasta, params.cdhit_id, params.cdhit_as)
        blast_ch = BLAST(indivcdhit_ch, params.vsg_db, params.notvsg_db)
        population_ch = CONCATENATE_VSGS((blast_ch.vsgs).collect(), params.full_vsg_db)
        catcdhit_ch = CONCATENATED_CDHIT(population_ch, params.cdhit_id, params.cdhit_as)
        index_ch = INDEX(catcdhit_ch, params.cores)
        quant_ch = QUANTIFY(index_ch, params.cores, trimmed_reads_ch)
        multiqc_ch = MULTIQC((quant_ch.quants).collect())
        summarise_ch = SUMMARISE((quant_ch.quants).collect(), params.threshold, (blast_ch.vsgs).collect(), catcdhit_ch.clstr)
    }  
    else if (params.mode == "assemble"){
        trimmed_reads_ch = TRIM(ch_reads, params.cores)
        assemble_ch = ASSEMBLE(trimmed_reads_ch, params.cores, params.trinitymem)
    }
    else if (params.mode == "predictvsgs"){
        assemblies_ch = Channel.fromPath(params.assemblies, checkIfExists: true)
        orf_ch = ORF(assemblies_ch, params.cdslength)
        indivcdhit_ch = INDIVIDUAL_CDHIT(orf_ch.fasta, params.cdhit_id, params.cdhit_as)
        blast_ch = BLAST(indivcdhit_ch, params.vsg_db, params.notvsg_db)
        population_ch = CONCATENATE_VSGS((blast_ch.vsgs).collect(), params.full_vsg_db)
        catcdhit_ch = CONCATENATED_CDHIT(population_ch, params.cdhit_id, params.cdhit_as)
    }
    else if (params.mode == "quantify"){
        index_ch = INDEX(params.vsgome, params.cores)
        quant_ch = QUANTIFY(index_ch, params.cores, ch_reads)
        multiqc_ch = MULTIQC((quant_ch.quants).collect())
        summarise_ch = SUMMARISEVSGOME((quant_ch.quants).collect())
    }
    else if (params.mode == "old_analyse"){
        assemblies_ch = Channel.fromPath(params.assemblies, checkIfExists: true)
        orf_ch = ORF(assemblies_ch, params.cdslength)
        indivcdhit_ch = INDIVIDUAL_CDHIT(orf_ch.fasta, params.cdhit_id, params.cdhit_as)
        blast_ch = BLAST(indivcdhit_ch, params.vsg_db, params.notvsg_db)
        population_ch = CONCATENATE_VSGS((blast_ch.vsgs).collect(), params.full_vsg_db)
        catcdhit_ch = CONCATENATED_CDHIT(population_ch, params.cdhit_id, params.cdhit_as)
        index_ch = INDEX(catcdhit_ch, params.cores)
        quant_ch = QUANTIFY(index_ch, params.cores, ch_reads)
        multiqc_ch = MULTIQC((quant_ch.quants).collect())
        summarise_ch = SUMMARISE((quant_ch.quants).collect(), params.threshold, (blast_ch.vsgs).collect(), catcdhit_ch.clstr)
    }
    else if (params.mode == "full") {
        trimmed_reads_ch = TRIM(ch_reads, params.cores)
        assemble_ch = ASSEMBLE(trimmed_reads_ch, params.cores, params.trinitymem)
        orf_ch = ORF(assemble_ch, params.cdslength)
        blast_ch = BLAST(orf_ch, params.vsg_db, params.notvsg_db)
        population_ch = CONCATENATE_VSGS((blast_ch.vsgs).collect(), params.full_vsg_db)
        catcdhit_ch = CONCATENATED_CDHIT(population_ch, params.cdhit_id, params.cdhit_as)
        index_ch = INDEX(population_ch, params.cores)
        quant_ch = QUANTIFY(index_ch, params.cores, trimmed_reads_ch)
        multiqc_ch = MULTIQC((quant_ch.quants).collect())
        summarise_ch = SUMMARISE((quant_ch.quants).collect(), params.threshold, (blast_ch.vsgs).collect(), catcdhit_ch.clstr, population_ch)
    }
    else if (params.mode == "analyse"){
        assemblies_ch = Channel.fromPath(params.assemblies, checkIfExists: true)
        orf_ch = ORF(assemblies_ch, params.cdslength)
        blast_ch = BLAST(orf_ch, params.vsg_db, params.notvsg_db)
        population_ch = CONCATENATE_VSGS((blast_ch.vsgs).collect(), params.full_vsg_db)
        catcdhit_ch = CONCATENATED_CDHIT(population_ch, params.cdhit_id, params.cdhit_as)
        index_ch = INDEX(population_ch, params.cores)
        quant_ch = QUANTIFY(index_ch, params.cores, ch_reads)
        multiqc_ch = MULTIQC((quant_ch.quants).collect())
        summarise_ch = SUMMARISE((quant_ch.quants).collect(), params.threshold, (blast_ch.vsgs).collect(), catcdhit_ch.clstr, population_ch)
    }  
    else {
        log.error("Invalid mode selected. Please select one of the following: full, assemble, predictvsgs, quantify, analyse.")
        System.exit(1)
}
}

workflow.onComplete {
    log.info(workflow.success ? "\nDone!" : "Oops .. something went wrong")
}