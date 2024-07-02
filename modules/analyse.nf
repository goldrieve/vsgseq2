#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


process INDIVIDUAL_CDHIT {
    input:
    path (assemblies)
    val cdhitid

    output:
    path "${assemblies.baseName.replace("_ORF","")}_cdhit.fasta"

    script:
    """
    cd-hit-est -i ${assemblies} -o ${assemblies.baseName.replace("_ORF","")}_cdhit.fasta -d 0 -c ${cdhitid} -n 8 -G 1 -g 1 -s 0.0 -aL 0.0 -M 50
    """
}

process BLAST {
    publishDir "${params.outdir}/VSGs", mode:'copy', pattern: '*_VSGs.fasta'
    input:
    path (assemblies)
    val (vsg_db)
    val (notvsg_db)

    output:
    path "${assemblies.baseName.replace("_cdhit","")}.xml", emit: vsgblast
    path "${assemblies.baseName.replace("_cdhit","")}_nonVSG.xml", emit: notvsgblast
    path "${assemblies.baseName.replace("_cdhit","")}_VSGs.fasta", emit: vsgs

    script:
    """
    blastn -db $vsg_db -query ${assemblies} -outfmt 5 -out ${assemblies.baseName.replace("_cdhit","")}.xml
    blastn -db $notvsg_db -query ${assemblies} -outfmt 5 -out ${assemblies.baseName.replace("_cdhit","")}_nonVSG.xml
    python $projectDir/bin/process_vsgs.py ${assemblies.baseName.replace("_cdhit","")}
    """
}

process CONCATENATE_VSGS {
    publishDir "${params.outdir}/VSGs", mode:'copy'
    input:
    path (vsgome)

    output:
    path "concatenated_vsgs.fasta"

    script:
    """
    cat ${vsgome} > concatenated_vsgs.fasta
    """
}

process CONCATENATED_CDHIT {
    publishDir "${params.outdir}/VSGs", mode:'copy'

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
    input:
    path (vsgs)
    val cores

    output:
    path "salmon_index"


    script:
    """
    cat $vsgs $projectDir/data/genome/TriTrypDB-67_TbruceiEATRO1125_Genome.fasta > gentrome.fasta
    salmon index --threads $cores -t gentrome.fasta -d $projectDir/data/genome/TriTrypDB-67_TbruceiEATRO1125_Genome_decoys.txt -i salmon_index --gencode
    """
}

process QUANTIFY {
    cpus = params.requestedcpus
    input:
    path (index)
    val cores
    tuple val(sample_id), path(reads)

    output:
    path "${reads[0].baseName.replace("_trimmed_1.fq","")}_quant", emit: quants


    script:
    """
    salmon quant --libType A --threads $cores -i $index -1 ${reads[0]} -2 ${reads[1]}  --validateMappings -o ${reads[0].baseName.replace("_trimmed_1.fq","")}_quant
    """
}

process MULTIQC {
    publishDir "${params.outdir}/summary", mode:'copy'
    cpus = params.requestedcpus

    input:
    path (quants)

    output:
    path "multiqc_report.html"


    script:
    """
    multiqc .
    """
}

process SUMMARISE {
    publishDir "${params.outdir}/summary", mode:'copy'
    input:
    val quants
    val vsgs

    output:
    path "tpm.csv"
    path "vsg_count.csv"


    script:
    """
    Rscript $projectDir/bin/summarise_quant.R "${quants}"
    Rscript $projectDir/bin/summarise_vsgs.R "${vsgs}"
    """
}


workflow {
    Channel
        .fromFilePairs(params.reads, checkIfExists: true)
        .set { files_ch }

    assemblies_ch = Channel.fromPath(params.assemblies, checkIfExists: true)
    longorf_ch = ORF(assemblies_ch, params.cdslength)
    indivcdhit_ch = INDIVIDUAL_CDHIT(longorf_ch.fasta, params.cdhitid)
    blast_ch = BLAST(indivcdhit_ch, params.vsg_db, params.notvsg_db)
    population_ch = CONCATENATE_VSGS((blast_ch.vsgs).collect())
    catcdhit_ch = CONCATENATED_CDHIT(population_ch, params.cdhitid)
    index_ch = INDEX(catcdhit_ch, params.cores)
    quant_ch = QUANTIFY(index_ch, params.cores, files_ch)
    multiqc_ch = MULTIQC((quant_ch.quants).collect())
    summarise_ch = SUMMARISE((quant_ch.quants).collect(), (blast_ch.vsgs).collect())
}

workflow.onComplete {
    log.info ( workflow.success ? "\nDone!" : "Oops .. something went wrong" )
}
