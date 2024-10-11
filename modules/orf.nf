process ORF {
    conda "$projectDir/vsgseq2.yml" 
    input:
    path (assemblies)
    val cdslength

    output:
    path "${assemblies.baseName.replace("_trinity.Trinity","")}_ORF.fasta", emit: fasta

    script:
    """
    TransDecoder.LongOrfs --complete_orfs_only -m ${cdslength} -t ${assemblies}
    TransDecoder.Predict --no_refine_starts --single_best_only -t ${assemblies}
    seqkit replace ${assemblies}.transdecoder.cds -p "(.+)" -r "${assemblies}{nr}" > ${assemblies}_transdecoder.fasta
    seqkit replace ${assemblies}_transdecoder.fasta -p "Trinity.fasta" -r '' > ${assemblies.baseName.replace("_trinity.Trinity","")}_ORF.fasta
    """
}
