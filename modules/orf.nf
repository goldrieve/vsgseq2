process ORF {
    conda "$projectDir/vsgseq2.yml"
    container 'goldrieve/vsgseq2:latest'
    input:
    path (assemblies)
    val cdslength
    val partial

    output:
    path "${assemblies.baseName.replace("_trinity.Trinity","")}_ORF.fasta", emit: fasta

    script:
    """
    TransDecoder.LongOrfs ${partial} -m ${cdslength} -t ${assemblies}
    TransDecoder.Predict --no_refine_starts --single_best_only -t ${assemblies}
    seqkit replace ${assemblies}.transdecoder.cds -p "(.+)" -r "${assemblies}{nr}" > ${assemblies}_transdecoder.fasta
    seqkit replace ${assemblies}_transdecoder.fasta -p "Trinity.fasta" -r '' > ${assemblies.baseName.replace("_trinity.Trinity","")}_ORF.fasta
    """
}
