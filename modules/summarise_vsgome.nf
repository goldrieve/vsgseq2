process SUMMARISEVSGOME {
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
