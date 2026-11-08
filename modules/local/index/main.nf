process INDEX {
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/salmon:1.10.3--6ea6d56a4b0c26b8"

    input:
        path (vsgs)
        val (cores)
        path (genome)
        path (decoys)

    output:
        path ("salmon_index"), emit: index
        tuple val("${task.process}"), val('salmon'), eval("salmon --version | sed 's/salmon //'"), topic: versions, emit: versions_salmon

    script:
        """
        cat ${vsgs} ${genome} > gentrome.fasta
        salmon index --threads ${cores} -t gentrome.fasta -d ${decoys} -i salmon_index --gencode
        """
}
