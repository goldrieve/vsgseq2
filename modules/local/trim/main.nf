process TRIM {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/trimmomatic:0.39--f4a7af7a8bfa5641"

    input:
        tuple val(meta), path(reads)
        val (cores)

    output:
        tuple val(meta), path("*_trimmed*.fq.gz"), emit: trimmed_reads
        tuple val("${task.process}"), val('trimmomatic'), eval("trimmomatic -version"), topic: versions, emit: versions_trimmomatic

    script:
        def args   = task.ext.args ?: ''
        def prefix = task.ext.prefix ?: "${meta.id}"
        if (meta.single_end) {
            """
            trimmomatic SE \\
                -threads ${cores} \\
                ${reads[0]} \\
                ${prefix}_trimmed.fq.gz \\
                SLIDINGWINDOW:4:5 LEADING:5 TRAILING:5 MINLEN:25 ${args}
            """
        } else {
            """
            trimmomatic PE \\
                -threads ${cores} \\
                ${reads[0]} ${reads[1]} \\
                ${prefix}_trimmed_1.fq.gz ${prefix}_unpaired_1.fq.gz \\
                ${prefix}_trimmed_2.fq.gz ${prefix}_unpaired_2.fq.gz \\
                SLIDINGWINDOW:4:5 LEADING:5 TRAILING:5 MINLEN:25 ${args}
            """
        }
}
