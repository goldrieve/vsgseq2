process SUMMARISE {
    conda "${params.conda_yml}"
    container 'goldrieve/vsgseq2:latest'
    
    publishDir "${params.outdir}/summary/tpm", mode:'copy', pattern: "tpm.csv"
    publishDir "${params.outdir}/summary/read_counts", mode:'copy', pattern: "num_reads.csv"
    publishDir "${params.outdir}/summary/read_counts", mode:'copy', pattern: "total_read_counts.csv"
    publishDir "${params.outdir}/summary/vsgs", mode:'copy', pattern: "vsg_count.csv"
    publishDir "${params.outdir}/summary/tpm", mode:'copy', pattern: "filtered_tpm.csv"
    publishDir "${params.outdir}/summary/cluster", mode:'copy', pattern: "filtered_tpm_clusters.csv"
    publishDir "${params.outdir}/summary/length", mode:'copy', pattern: "length.csv"
    publishDir "${params.outdir}/summary/cluster", mode:'copy', pattern: "filtered_tpm_clusters_length.csv"
    publishDir "${params.outdir}/summary/tpm", mode:'copy', pattern: "cluster_tpm.csv"
    publishDir "${params.outdir}/summary/cluster", mode:'copy', pattern: "cluster_champion.csv"
    publishDir "${params.outdir}/summary/cluster", mode:'copy', pattern: "champion_vsgs.fasta"
    
    input:
        val quants
        val threshold
        val vsgs
        val clstr
        val fasta

    output:
        path "tpm.csv"
        path "num_reads.csv"
        path "total_read_counts.csv"
        path "vsg_count.csv"
        path "filtered_tpm.csv"
        path "filtered_tpm_clusters.csv"
        path "length.csv"
        path "filtered_tpm_clusters_length.csv"
        path "cluster_tpm.csv"
        path "cluster_champion.csv"
        path "champion_vsgs.fasta", emit: champion_vsgs


    script:
        """
        Rscript $projectDir/bin/summarise_quant.R "${quants}" ${threshold}
        Rscript $projectDir/bin/summarise_vsgs.R "${vsgs}"
        python $projectDir/bin/add_cluster.py filtered_tpm.csv ${clstr} filtered_tpm_clusters.csv
        python $projectDir/bin/length.py "${fasta}" length.csv
        python $projectDir/bin/merge_length_tpm.py filtered_tpm.csv length.csv filtered_tpm_clusters_length.csv
        python $projectDir/bin/sum_cluster.py filtered_tpm_clusters.csv cluster_tpm.csv cluster_champion.csv "${fasta}" champion_vsgs.fasta
        """
}
