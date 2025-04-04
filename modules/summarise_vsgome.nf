process SUMMARISEVSGOME {
    conda "$projectDir/vsgseq2.yml"
    container 'goldrieve/vsgseq2:latest'
    publishDir "${params.outdir}/summary", mode:'copy'
    input:
    val quants

    output:
    path "tpm.csv"


    script:
    """
    Rscript $projectDir/bin/summarise_quant.R "${quants}"
    """
}