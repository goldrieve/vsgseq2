process MULTIQC {
    publishDir "${params.outdir}/summary", mode:'copy'

    input:
        path quants

    output:
        path "multiqc_report.html"
        path "multiqc_data"


    script:
        """
        multiqc .
        """
}
