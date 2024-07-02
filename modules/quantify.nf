process QUANTIFY {
    cpus = params.requestedcpus
    input:
    path (index)
    val cores
    tuple val(meta), path(trimmed)

    output:
    path "*_quant", emit: quants

    script:
    if(meta.single_end){
        """
        salmon quant --libType A --threads $cores -i $index -r ${trimmed} --validateMappings -o ${trimmed.baseName.replace(".fq","")}_quant
        """
    } else {
        """
        salmon quant --libType A --threads $cores -i $index -1 ${trimmed[0]} -2 ${trimmed[1]}  --validateMappings -o ${trimmed[0].baseName.replace("_trimmed_1.fq","")}_quant
        """
    }
}