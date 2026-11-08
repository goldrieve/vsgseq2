process CONCATENATED_CDHIT {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/cd-hit:4.8.1--847ceda1bc4ad3ba"

    input:
        path (vsgs)
        val (cdhit_id)
        val (cdhit_as)

    output:
        path ("VSGome.fasta"), emit: vsgome
        path ("VSGome.fasta.clstr"), emit: clstr
        tuple val("${task.process}"), val('cd-hit'), eval("cd-hit-est -h | sed -n '1s/.*CD-HIT version \\([0-9.]*\\).*/\\1/p'"), topic: versions, emit: versions_cdhit

    script:
        def args = task.ext.args ?: ''
        """
        cd-hit-est \\
            -i ${vsgs} \\
            -o VSGome.fasta \\
            -d 0 \\
            -c ${cdhit_id} \\
            -n 8 \\
            -G 0 \\
            -g 1 \\
            -s 0.0 \\
            -aS ${cdhit_as} \\
            -t 4 \\
            -M 2000 \\
            ${args}
        """
}
