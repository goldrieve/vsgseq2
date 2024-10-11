process SUMMARISE {
    conda "$projectDir/vsgseq2.yml"
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
