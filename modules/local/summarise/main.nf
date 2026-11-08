process SUMMARISE {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/biopython_pandas_python_r-base_pruned:9905093b1f3b0d9c"

    input:
        path (quants)
        val (threshold)
        path (vsgs)
        path (clstr)
        path (fasta)

    output:
        path ("tpm.csv"), emit: tpm
        path ("num_reads.csv"), emit: num_reads
        path ("total_read_counts.csv"), emit: total_reads
        path ("vsg_count.csv"), emit: vsg_count
        path ("filtered_tpm.csv"), emit: filtered_tpm
        path ("filtered_tpm_clusters.csv"), emit: tpm_clusters
        path ("length.csv"), emit: length
        path ("filtered_tpm_clusters_length.csv"), emit: cluster_tpm_length
        path ("cluster_tpm.csv"), emit: cluster_tpm
        path ("cluster_champion.csv"), emit: cluster_champ
        path ("champion_vsgs.fasta"), emit: champ_vsgs
        tuple val("${task.process}"), val('r-base'), eval("Rscript -e 'cat(as.character(getRversion()))'"), topic: versions, emit: versions_r
        tuple val("${task.process}"), val('pandas'), eval("python -c 'import pandas; print(pandas.__version__)'"), topic: versions, emit: versions_pandas

    script:
        def quants_array = quants instanceof List ? "[${quants.join(', ')}]" : "[${quants}]"
        def vsgs_array   = vsgs instanceof List ? "[${vsgs.join(', ')}]" : "[${vsgs}]"
        """
        summarise_quant.R "${quants_array}" ${threshold}
        summarise_vsgs.R "${vsgs_array}"
        add_cluster.py filtered_tpm.csv ${clstr} filtered_tpm_clusters.csv
        length.py "${fasta}" length.csv
        merge_length_tpm.py filtered_tpm.csv length.csv filtered_tpm_clusters_length.csv
        sum_cluster.py filtered_tpm_clusters.csv cluster_tpm.csv cluster_champion.csv "${fasta}" champion_vsgs.fasta
        """
}
