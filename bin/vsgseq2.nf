#!/usr/bin/env nextflow

/*
 * pipeline input parameters
 */
params.reads = "$projectDir/data/reads/*{1,2}.fq.gz"
params.vsg_db = "concatAnTattb427.fa"
params.notvsg_db = "NOTvsgs.fa"
params.cores = "4"
params.trinitymem = "20"
params.cdslength = "300"
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

process TRIM {
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
    tuple val(sample_id), path(reads)
    val cores 
    val trinitymem

    output:
    path "${sample_id}*{1,2}.fq.gz.P.qtrim.gz", emit: trimmed
    path "${sample_id}", emit: dir
    path "${sample_id}_trinity.Trinity.fasta", emit: trinity_fasta

    script:
    """
    Trinity --seqType fq --left ${reads[0]}  --right ${reads[1]} --CPU ${cores} --max_memory ${trinitymem}G --no_path_merging --min_kmer_cov 2 --no_parallel_norm_stats --output ${sample_id}_trinity
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
    path "${sample_id}_ORF.fasta"

    script:
    """
    TransDecoder.LongOrfs -m ${cdslength} -t ${sample_id}
    TransDecoder.Predict --no_refine_starts --single_best_only -t ${sample_id}
    seqkit replace ${sample_id}.transdecoder.cds -p "(.+)" -r "${sample_id}{nr}" > ${sample_id}_sample.fasta
    seqkit replace ${sample_id}_sample.fasta -p "Trinity.fasta" -r '' > ${sample_id}_ORF.fasta
    """
}

process INDIVIDUAL_CDHIT {
    tag "INDIVIDUAL_CDHIT on $sample_id"
    publishDir params.outdir, mode:'copy'

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

process BLAST {
    tag "BLAST on $sample_id"
    publishDir params.outdir, mode:'copy'

    input:
    path (sample_id)
    val vsg_db
    val notvsg_db

    output:
    path "${sample_id}.xml", emit: vsgblast
    path "${sample_id}.nonVSG.xml", emit: notvsgblast


    script:
    """
    blastn -db ../../../data/blastdb/$vsg_db -query ${sample_id} -outfmt 5 -out ${sample_id}.xml
    blastn -db ../../../data/blastdb/$notvsg_db -query ${sample_id} -outfmt 5 -out ${sample_id}.nonVSG.xml
    """
}

process PROCESS_BLAST {
    tag "PROCESS_BLAST on $sample_id"
    publishDir params.outdir, mode:'copy'

    input:
    path (dir)
    path (vsgblast)
    path (notvsgblast)

    output:
    path "${sample_id}_VSGs.fasta"

    script:
    """
    python ../../../bin/process_vsgs.py ${dir}
    """
}

process CONCATENATE_VSGS {
    tag "CONCATENATE_VSGS on $sample_id"
    publishDir params.outdir, mode:'copy'

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
    cd-hit-est -i ${vsgs} -o VSGome.fasta -d 0 -c ${cdhitid} -n 8 -G 1 -g 1 -s 0.0 -aL 0.0 -M 5000
    """
}

process INDEX {
    tag "INDEX is running"
    publishDir params.outdir, mode:'copy'

    input:
    path (vsgs)
    val cores

    output:
    path "salmon_index"


    script:
    """
    cat $vsgs ../../../data/genome/TriTrypDB-67_TbruceiEATRO1125_Genome.fasta > gentrome.fasta
    salmon index --threads $cores -t gentrome.fasta -d ../../../data/genome/TriTrypDB-67_TbruceiEATRO1125_Genome_decoys.txt -i salmon_index --gencode
    """
}

process QUANTIFY {
    tag "QUANTIFY on $trimmed"
    publishDir params.outdir, mode:'copy'
    cpus = params.requestedcpus

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
    indivcdhit_ch = INDIVIDUAL_CDHIT(longorf_ch, params.cdhitid)
    blast_ch = BLAST(indivcdhit_ch, params.vsg_db, params.notvsg_db)
    process_blast_ch = PROCESS_BLAST (assemble_ch.dir, blast_ch.vsgblast, blast_ch.notvsgblast)
    population_ch = CONCATENATE_VSGS(params.vsg_db, (process_blast_ch).collect())
    catcdhit_ch = CONCATENATED_CDHIT(population_ch, params.cdhitid)
    index_ch = INDEX(catcdhit_ch, params.cores)
    quant_ch = QUANTIFY(index_ch, params.cores, assemble_ch.trimmed, assemble_ch.dir)
}

workflow.onComplete {
    log.info ( workflow.success ? "\nDone!" : "Oops .. something went wrong" )
}