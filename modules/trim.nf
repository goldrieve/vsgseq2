process TRIM {
    publishDir "${params.outdir}/trimmed", mode:'copy'
    input:
    tuple val(meta), path(reads)
    val cores

    output:
        path "*trimmed*.fq.gz", emit: trimmed_reads


    script:
        if(meta.single_end){
        """
        trimmomatic SE -threads ${cores} ${reads[0]} ${reads[0].baseName.replace(".fq","")}_trimmed.fq.gz SLIDINGWINDOW:4:5 LEADING:5 TRAILING:5 MINLEN:25
        """
    } else {
        """
        trimmomatic PE -threads ${cores} ${reads[0]} ${reads[1]} ${reads[0].baseName.replace("_1.fq","")}_trimmed_1.fq.gz ${reads[0].baseName}_unpaired_1.fq.gz ${reads[0].baseName.replace("_1.fq","")}_trimmed_2.fq.gz ${reads[0].baseName}_unpaired_2.fq.gz SLIDINGWINDOW:4:5 LEADING:5 TRAILING:5 MINLEN:25
        """
    }
}