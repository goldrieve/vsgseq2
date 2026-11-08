process QUANTIFY {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/salmon:1.10.3--6ea6d56a4b0c26b8"

    input:
        path (index)
        val (cores)
        tuple val(meta), path(trimmed)

    output:
        tuple val(meta), path("${meta.id}_quant"), emit: quants
        tuple val("${task.process}"), val('salmon'), eval("salmon --version | sed 's/salmon //'"), topic: versions, emit: versions_salmon

    script:
        def reads = meta.single_end ? "-r ${trimmed}" : "-1 ${trimmed[0]} -2 ${trimmed[1]}"
        """
        salmon quant \\
            --libType A \\
            --threads ${cores} \\
            -i ${index} \\
            ${reads} \\
            --validateMappings \\
            -o ${meta.id}_quant
        """
}
