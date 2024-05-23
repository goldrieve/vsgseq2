#!/usr/bin/env nextflow

/*
 * pipeline input parameters
 */
params.reads = "$projectDir/data/reads/*{1,2}.fq.gz"
params.cores = "4"
params.trinitymem = "20"
params.outdir = "results"
params.help = ""
params.requestedcpus = 4

log.info """\

   |======================================|
   | V S G S E Q 2 - N F - A S S E M B L E|
   |======================================|
    """
    .stripIndent()

/*
 * define the `INDEX` process that creates a binary index
 * given the transcriptome file
 */

if ( params.help ) {
    help = """VSGSEQ2.nf: A pipeline for analysing VSGSeq data
             |   
             |Required arguments:
             |
             |  --reads Location of reads, if not in reads dir
             |                [default: ${params.reads}]
             |
             |Optional arguments:
             |
             |  --requestedcpus  Define number of cores VSGSeq2 will use.
             |                [default: ${params.requestedcpus}]
             |  --cores  Define number of cores Trinity will use.
             |                [default: ${params.cores}]
             |  --trinitymem    Define mem Trinity will use.
             |                [default: ${params.trinitymem}]
             |  --outdir        VSGSeq outdir. 
             |                [default: ${params.outdir}]
	         |  --help         Print this message.""".stripMargin()
    // Print the help with the stripped margin and exit
    println(help)
    exit(0)
}

process TRIM {
    publishDir params.outdir, mode:'copy'
    input:
    tuple val(sample_id), path(reads)

    output:
    path "${reads[0].baseName.replace("_1.fq","")}_trimmed_{1,2}.fq.gz", emit: trimmed

    script:
    """
    trimmomatic PE ${reads[0]} ${reads[1]} ${reads[0].baseName.replace("_1.fq","")}_trimmed_1.fq.gz ${reads[0].baseName}_unpaired.fq.gz ${reads[0].baseName.replace("_1.fq","")}_trimmed_2.fq.gz ${reads[0].baseName}_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:5 LEADING:5 TRAILING:5 MINLEN:25
    """
}

process ASSEMBLE {
    tag "ASSEMBLE on $sample_id"
    publishDir params.outdir, pattern: '*trinity.Trinity.fasta', mode:'copy'
    cpus = params.requestedcpus
    
    input:
    path trimmed
    val cores 
    val trinitymem

    output:
    path "${sample_id}_trinity.Trinity.fasta", emit: trinity_fasta

    script:
    """
    Trinity --seqType fq --left ${reads[0]}  --right ${reads[1]} --CPU ${cores} --max_memory ${trinitymem}G --no_path_merging --min_kmer_cov 2 --no_parallel_norm_stats --output ${sample_id}_trinity
    rm -r ${sample_id}_trinity
    """
}

workflow {
    Channel
        .fromFilePairs(params.reads, checkIfExists: true)
        .set { read_pairs_ch }

    trim_ch = TRIM(read_pairs_ch, params.cores, params.trinitymem)
    assemble_ch = ASSEMBLE(trim_ch, params.cores, params.trinitymem)
}

workflow.onComplete {
    log.info ( workflow.success ? "\nDone!" : "Oops .. something went wrong" )
}
