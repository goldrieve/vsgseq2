process MULTIQC {
    conda "$projectDir/vsgseq2.yml"
    container 'goldrieve/vsgseq2:latest' 
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
