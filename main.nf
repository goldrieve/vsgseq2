nextflow.enable.dsl = 2

params.assemblies = "$projectDir/results/assemblies/*_trinity.Trinity.fasta"
params.reads = "$projectDir/data/reads/*{1,2}.fq.gz"
params.vsg_db = "$projectDir/data/blastdb/concatAnTattb427.fa"
params.notvsg_db = "$projectDir/data/blastdb/NOTvsgs.fa"
params.vsgome = "$projectDir/data/blastdb/concatAnTattb427.fa"
params.full_vsg_db = ""
params.requestedcpus = 4
params.cores = "4"
params.trinitymem = "20"
params.cdslength = "300"
params.cdhitid = "0.98"
params.outdir = "results"
params.samplesheet = "samplesheet.csv"
params.help = ""
params.mode = "full"

if (params.help) {
    help = """VSGSEQ2.nf: A pipeline for analysing VSGSeq data
             |
             |Required arguments:
             |
             |  --assemblies Location of assemblies
             |                [default: ${params.assemblies}]
             |  --reads Location of reads, if not in reads dir
             |                [default: ${params.reads}]
             |  --vsg_db    Location of VSGdb
             |                [default: ${params.vsg_db}]
             |  --notvsg_db Location of NOTVSGdb
             |                [default: ${params.notvsg_db}]
             |  --vsgome    Location of VSGome
             |                [default: ${params.vsgome}]
             |  --full_vsg_db    Location of a database to add into the VSGome (such as data/blastdb/concatAnTattb427_full.fa). 
             |                    Default will only include the assembled VSGome.
             |  --mode    The mode to run the pipeline in. Options are full, assemble, predictvsgs, quantify, analyse.
             |                [default: ${params.mode}]
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
             |  --cdhitid       Define sequence identity threshold - how much the alignment has to match (0.0 - 1.0).
             |                [default: ${params.cdhitid}]
             |  --outdir        VSGSeq output directory. 
             |                [default: ${params.outdir}]
             |  --samplesheet  Define the path to the samplesheet.
             |                [default: ${params.samplesheet}]
             |  --help         Print this message.""".stripMargin()
    // Print the help with the stripped margin and exit
    println(help)
    exit(0)
}

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

ch_reads = ch_samplesheet.splitCsv(header:true).map {

    r1 = it['r1']
    r2 = it['r2']

    is_singleEnd = r2.toString()=='' ? true : false

    meta = [id: it['sample'], single_end: is_singleEnd]

    r2.toString()=='' ? [meta, [r1]] : [meta, [r1, r2]]

}

workflow {
    if (params.mode == "full") {
        trimmed_reads_ch = TRIM(ch_reads, params.cores)
        assemble_ch = ASSEMBLE(trimmed_reads_ch, params.cores, params.trinitymem)
        orf_ch = ORF(assemble_ch, params.cdslength)
        indivcdhit_ch = INDIVIDUAL_CDHIT(orf_ch.fasta, params.cdhitid)
        blast_ch = BLAST(indivcdhit_ch, params.vsg_db, params.notvsg_db)
        population_ch = CONCATENATE_VSGS((blast_ch.vsgs).collect(), params.full_vsg_db)
        catcdhit_ch = CONCATENATED_CDHIT(population_ch, params.cdhitid)
        index_ch = INDEX(catcdhit_ch, params.cores)
        quant_ch = QUANTIFY(index_ch, params.cores, trimmed_reads_ch)
        multiqc_ch = MULTIQC((quant_ch.quants).collect())
        summarise_ch = SUMMARISE((quant_ch.quants).collect(), (blast_ch.vsgs).collect())
    }  
    else if (params.mode == "assemble"){
        trimmed_reads_ch = TRIM(ch_reads, params.cores)
        assemble_ch = ASSEMBLE(trimmed_reads_ch, params.cores, params.trinitymem)
    }
    else if (params.mode == "predictvsgs"){
        assemblies_ch = Channel.fromPath(params.assemblies, checkIfExists: true)
        orf_ch = ORF(assemblies_ch, params.cdslength)
        indivcdhit_ch = INDIVIDUAL_CDHIT(orf_ch.fasta, params.cdhitid)
        blast_ch = BLAST(indivcdhit_ch, params.vsg_db, params.notvsg_db)
        population_ch = CONCATENATE_VSGS((blast_ch.vsgs).collect(), params.full_vsg_db)
        catcdhit_ch = CONCATENATED_CDHIT(population_ch, params.cdhitid)
    }
    else if (params.mode == "quantify"){
        index_ch = INDEX(params.vsgome, params.cores)
        quant_ch = QUANTIFY(index_ch, params.cores, ch_reads)
        multiqc_ch = MULTIQC((quant_ch.quants).collect())
        summarise_ch = SUMMARISEVSGOME((quant_ch.quants).collect())
    }
    else if (params.mode == "analyse"){
        assemblies_ch = Channel.fromPath(params.assemblies, checkIfExists: true)
        orf_ch = ORF(assemblies_ch, params.cdslength)
        indivcdhit_ch = INDIVIDUAL_CDHIT(orf_ch.fasta, params.cdhitid)
        blast_ch = BLAST(indivcdhit_ch, params.vsg_db, params.notvsg_db)
        population_ch = CONCATENATE_VSGS((blast_ch.vsgs).collect(), params.full_vsg_db)
        catcdhit_ch = CONCATENATED_CDHIT(population_ch, params.cdhitid)
        index_ch = INDEX(catcdhit_ch, params.cores)
        quant_ch = QUANTIFY(index_ch, params.cores, ch_reads)
        multiqc_ch = MULTIQC((quant_ch.quants).collect())
        summarise_ch = SUMMARISE((quant_ch.quants).collect(), (blast_ch.vsgs).collect())
    }
    else {
        log.error("Invalid mode selected. Please select one of the following: full, assemble, predictvsgs, quantify, analyse.")
        System.exit(1)
}
}

workflow.onComplete {
    log.info(workflow.success ? "\nDone!" : "Oops .. something went wrong")
}