process SUMMARISE {
    
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
        path quants
        val threshold
        path vsgs
        path clstr
        path fasta

    output:
        path "tpm.csv", emit: tpm
        path "num_reads.csv", emit: num_reads
        path "total_read_counts.csv", emit: total_reads
        path "vsg_count.csv", emit: vsg_count
        path "filtered_tpm.csv", emit: filtered_tpm
        path "filtered_tpm_clusters.csv", emit: tpm_clusters
        path "length.csv", emit: length
        path "filtered_tpm_clusters_length.csv", emit: cluster_tpm_length
        path "cluster_tpm.csv", emit: cluster_tpm
        path "cluster_champion.csv", emit: cluster_champ
        path "champion_vsgs.fasta", emit: champ_vsgs


    script:
        """
        Rscript ${params.scripts}summarise_quant.R "${quants}" ${threshold}
        Rscript ${params.scripts}summarise_vsgs.R "${vsgs}"
        python ${params.scripts}add_cluster.py filtered_tpm.csv ${clstr} filtered_tpm_clusters.csv
        python ${params.scripts}length.py "${fasta}" length.csv
        python ${params.scripts}merge_length_tpm.py filtered_tpm.csv length.csv filtered_tpm_clusters_length.csv
        python ${params.scripts}sum_cluster.py filtered_tpm_clusters.csv cluster_tpm.csv cluster_champion.csv "${fasta}" champion_vsgs.fasta
        """
}
