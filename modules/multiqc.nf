process MULTIQC {
    conda "${params.conda_yml}"
    container 'goldrieve/vsgseq2:latest' 
    publishDir "${params.outdir}/summary", mode:'copy'

    input:
        path (quants)

    output:
        path "multiqc_report.html"


    script:
        """
        multiqc .
        """
}
