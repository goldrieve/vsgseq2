process SUMMARISE {
    conda "$projectDir/vsgseq2.yml"
    container 'goldrieve/vsgseq2:latest'
    publishDir "${params.outdir}/summary", mode:'copy'
    input:
    val quants
    val threshold
    val vsgs
    val clstr

    output:
    path "tpm.csv"
    path "num_reads.csv"
    path "total_read_counts.csv"
    path "vsg_count.csv"
    path "filtered_tpm.csv"
    path "filtered_tpm_clusters.csv"


    script:
    """
    Rscript $projectDir/bin/summarise_quant.R "${quants}" ${threshold}
    Rscript $projectDir/bin/summarise_vsgs.R "${vsgs}"
    python $projectDir/bin/add_cluster.py filtered_tpm.csv ${clstr} filtered_tpm_clusters.csv
    """
}
