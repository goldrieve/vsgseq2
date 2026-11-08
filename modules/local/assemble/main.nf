process ASSEMBLE {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/salmon_trinity:f0bfbacfba988caf"

    input:
        tuple val(meta), path(trimmed)
        val (cores)
        val (trinitymem)

    output:
        tuple val(meta), path("${meta.id}_trinity.Trinity.fasta"), emit: trinity_fasta
        tuple val("${task.process}"), val('trinity'), eval("Trinity --version | sed -n 's/^Trinity version: Trinity-v//p'"), topic: versions, emit: versions_trinity

    script:
        def args   = task.ext.args ?: ''
        def prefix = "${meta.id}_trinity"
        def reads  = meta.single_end ? "--single ${trimmed}" : "--left ${trimmed[0]} --right ${trimmed[1]}"
        """
        Trinity \\
            --seqType fq \\
            ${reads} \\
            --CPU ${cores} \\
            --max_memory ${trinitymem}G \\
            --no_path_merging \\
            --min_kmer_cov 2 \\
            --no_parallel_norm_stats \\
            --output ${prefix} \\
            ${args}

        rm -r ${prefix}
        """
}
