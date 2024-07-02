process ASSEMBLE {
    publishDir "${params.outdir}/assemblies", pattern: '*trinity.Trinity.fasta', mode:'copy'
    cpus = params.requestedcpus
    
    input:
    tuple val(meta), path(trimmed)
    val cores 
    val trinitymem

    output:
    path "*_trinity.Trinity.fasta", emit: trinity_fasta

    script:
    if(meta.single_end){
        """
        Trinity --seqType fq --single ${trimmed} --CPU ${cores} --max_memory ${trinitymem}G --no_path_merging --min_kmer_cov 2 --no_parallel_norm_stats --output ${trimmed.baseName.replace("_single.fq","")}_trinity
        rm -r ${trimmed.baseName.replace("_single.fq","")}_trinity
        """
    } else {
        """
        Trinity --seqType fq --left ${trimmed[0]}  --right ${trimmed[1]} --CPU ${cores} --max_memory ${trinitymem}G --no_path_merging --min_kmer_cov 2 --no_parallel_norm_stats --output ${trimmed[0].baseName.replace("_1.fq","")}_trinity
        rm -r ${trimmed[0].baseName.replace("_1.fq","")}_trinity
        """
    }
}