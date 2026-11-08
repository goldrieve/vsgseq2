process MULTIQC {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/multiqc:1.21--5801f97837e41c3b"

    input:
        path (quants)

    output:
        path ("multiqc_report.html"), emit: report
        path ("multiqc_data"), emit: data
        tuple val("${task.process}"), val('multiqc'), eval("multiqc --version | sed 's/.* //'"), topic: versions, emit: versions_multiqc

    script:
        def args = task.ext.args ?: ''
        """
        multiqc . ${args}
        """
}
