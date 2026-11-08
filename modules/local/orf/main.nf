process ORF {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/seqkit_transdecoder:cd76ae97f7c6e5fd"

    input:
        tuple val(meta), path(assembly)
        val (cdslength)
        val (partial)

    output:
        tuple val(meta), path("${meta.id}_ORF.fasta"), emit: fasta
        tuple val("${task.process}"), val('transdecoder'), eval("TransDecoder.LongOrfs --version | sed 's/^TransDecoder.LongOrfs //'"), topic: versions, emit: versions_transdecoder
        tuple val("${task.process}"), val('seqkit'), eval("seqkit version | sed 's/^seqkit v//'"), topic: versions, emit: versions_seqkit

    script:
        """
        TransDecoder.LongOrfs ${partial ? '' : '--complete_orfs_only'} -m ${cdslength} -t ${assembly}
        TransDecoder.Predict --no_refine_starts --single_best_only -t ${assembly}
        seqkit replace ${assembly}.transdecoder.cds -p "(.+)" -r "${assembly}{nr}" > ${assembly}_transdecoder.fasta
        seqkit replace ${assembly}_transdecoder.fasta -p "Trinity.fasta" -r '' > ${meta.id}_ORF.fasta
        """
}
