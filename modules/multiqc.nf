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
