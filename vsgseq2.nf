#!/usr/bin/env nextflow

/*
 * pipeline input parameters
 */
params.reads = "$projectDir/data/reads/*{1,2}.fq.gz"
params.vsg_db = "concatAnTattb427.fa"
params.notvsg_db = "NOTvsgs.fa"
params.cores = "4"
params.trinitymem = "20"
params.cdslength = "100"
params.cdhitid = "0.98"
params.outdir = "results"
params.help = ""
params.requestedcpus = 4

log.info """\

   |====================|
   | V S G S E Q 2 - N F|
   |====================|
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
             |  --vsg_db    Location of VSGdb
             |                [default: ${params.vsg_db}]
             |  --NOTVSG_db Location of NOTVSGdb
             |                [default: ${params.notvsg_db}]
             |
             |Optional arguments:
             |
             |  --requestedcpus  Define number of cores VSGSeq2 will use.
             |                [default: ${params.requestedcpus}]
             |  --cores  Define number of cores Trinity will use.
             |                [default: ${params.cores}]
             |  --trinitymem    Define mem Trinity will use.
             |                [default: ${params.trinitymem}]
             |  --cdslength    Define minimium CDS length (amino acids).
             |                [default: ${params.cdslength}]
             |  --cdhitid       Define sequence identiy threshold - how much the alignment has to match (0.0 - 1.0).
             |                [default: ${params.cdhitid}]
             |  --outdir        VSGSeq outdir. 
             |                [default: ${params.outdir}]
	         |  --help         Print this message.""".stripMargin()
    // Print the help with the stripped margin and exit
    println(help)
    exit(0)
}

process ASSEMBLE {
    tag "ASSEMBLE on $sample_id"
    publishDir params.outdir, pattern: '*trinity.Trinity.fasta', mode:'copy'
    cpus = params.requestedcpus
    
    input:
    tuple val(sample_id), path(reads)
    val cores 
    val trinitymem

    output:
    path "${sample_id}*{1,2}.fq.gz.P.qtrim.gz", emit: trimmed
    path "${sample_id}", emit: dir
    path "${sample_id}_trinity.Trinity.fasta", emit: trinity_fasta

    script:
    """
    Trinity --seqType fq --left ${reads[0]}  --right ${reads[1]} --CPU ${cores} --max_memory ${trinitymem}G --no_path_merging --no_normalize_reads --trimmomatic --output ${sample_id}_trinity
    mv ${sample_id}_trinity/${sample_id}*.fq.gz.P.qtrim.gz ./
    rm -r ${sample_id}_trinity
    mkdir ${sample_id}
    """
}

process ORF {
    tag "ORF on $sample_id"
    publishDir params.outdir, mode:'copy'

    input:
    path (sample_id)
    val cdslength

    output:
    path "${sample_id}_renamed.fasta"

    script:
    """
    TransDecoder.LongOrfs -m ${cdslength} -t ${sample_id}
    TransDecoder.Predict --single_best_only -t ${sample_id}
    seqkit replace ${sample_id}.transdecoder.cds -p "(.+)" -r "${sample_id}{nr}" > ${sample_id}_sample.fasta
    seqkit replace ${sample_id}_sample.fasta -p "Trinity.fasta" -r '' > ${sample_id}_renamed.fasta
    """
}

process BLAST_VSG {
    tag "BLAST_VSG on $sample_id"

    input:
    path (sample_id)
    val vsg_db

    output:
    path "${sample_id}_VSGBLAST.fasta"

    script:
    """
    blastn -db ../../../data/blastdb/$vsg_db -query ${sample_id} -outfmt "6 qseqid" -qcov_hsp_perc 30 -perc_identity 90 -evalue 1.0e-10 -max_target_seqs 1 -out ${sample_id}_BLAST_VSG.txt
    seqkit grep -f ${sample_id}_BLAST_VSG.txt ${sample_id} -o ${sample_id}_VSGBLAST.fasta
    """
}

process BLAST_NOT {
    tag "BLAST_NOT on $sample_id"

    input:
    path (sample_id)
    val notvsg_db

    output:
    path "${sample_id}_NOTBLAST.fasta"

    script:
    """
    blastn -db ../../../data/blastdb/$notvsg_db -query ${sample_id} -outfmt "6 qseqid" -qcov_hsp_perc 30 -perc_identity 90 -evalue 1.0e-10 -max_target_seqs 1 -out ${sample_id}_BLAST_NOT.txt
    seqkit grep -v -f ${sample_id}_BLAST_NOT.txt ${sample_id} -o ${sample_id}_NOTBLAST.fasta
    """
}

process INDIVIDUAL_CDHIT {
    tag "INDIVIDUAL_CDHIT on $sample_id"

    input:
    path (sample_id)
    val cdhitid

    output:
    path "${sample_id}_sample.fasta"

    script:
    """
    cd-hit-est -i ${sample_id} -o ${sample_id}_sample.fasta -d 0 -c ${cdhitid} -n 8 -G 1 -g 1 -s 0.0 -aL 0.0 -M 50
    """
}

process CONCATENATE_VSGS {
    tag "CONCATENATE_VSGS on $sample_id"

    input:
    val vsg_db
    path (sample_id)

    output:
    path "concatenated_vsgs.fasta"

    script:
    """
    cat ../../../data/blastdb/$vsg_db ${sample_id} > concatenated_vsgs.fasta
    """
}

process CONCATENATED_CDHIT {
    tag "CONCATENATED_CDHIT is running"
    publishDir params.outdir, mode:'copy'

    input:
    path (vsgs)
    val cdhitid

    output:
    path "VSGome.fasta"

    script:
    """
    cd-hit-est -i ${vsgs} -o VSGome.fasta -d 0 -c ${cdhitid} -n 8 -G 1 -g 1 -s 0.0 -aL 0.0 -M 50
    """
}

process INDEX {
    tag "INDEX is running"

    input:
    path (vsgs)
    val cores

    output:
    path "salmon_index"


    script:
    """
    cat $vsgs ../../../data/genome/TriTrypDB-67_TbruceiEATRO1125_Genome.fasta > gentrome.fasta
    salmon index --threads $cores -t gentrome.fasta -d ../../../data/genome//TriTrypDB-67_TbruceiEATRO1125_Genome_decoys.txt -i salmon_index --gencode
    """
}

process QUANTIFY {
    tag "QUANTIFY on $trimmed"
    publishDir params.outdir, mode:'copy'

    input:
    path (index)
    val cores
    path (trimmed)
    path (trinity)

    output:
    path "${trinity}_quant"


    script:
    """
    salmon quant --libType A --threads $cores -i $index -1 ${trimmed[0]} -2 ${trimmed[1]}  --validateMappings -o ${trinity}_quant
    """
}


workflow {
    Channel
        .fromFilePairs(params.reads, checkIfExists: true)
        .set { read_pairs_ch }

    assemble_ch = ASSEMBLE(read_pairs_ch, params.cores, params.trinitymem)
    longorf_ch = ORF(assemble_ch.trinity_fasta, params.cdslength)
    blastvsg_ch = BLAST_VSG(longorf_ch, params.vsg_db)
    blastnot_ch = BLAST_NOT(blastvsg_ch, params.notvsg_db)
    indivcdhit_ch = INDIVIDUAL_CDHIT(blastnot_ch, params.cdhitid)
    population_ch = CONCATENATE_VSGS(params.vsg_db, (indivcdhit_ch).collect())
    catcdhit_ch = CONCATENATED_CDHIT(population_ch, params.cdhitid)
    index_ch = INDEX(catcdhit_ch, params.cores)
    quant_ch = QUANTIFY(index_ch, params.cores, assemble_ch.trimmed, assemble_ch.dir)
}

workflow.onComplete {
    log.info ( workflow.success ? "\nDone!" : "Oops .. something went wrong" )
}
