process CONCATENATE_VSGS {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/coreutils:9.5--af6fea6f985629a9"

    input:
        path (vsgs, stageAs: 'vsgs/*')
        path (full_vsg_db, stageAs: 'full_vsg_db/*')

    output:
        path ("concatenated_vsgs.fasta"), emit: fasta
        tuple val("${task.process}"), val('coreutils'), eval("cat --version | sed -n '1s/.* //p'"), topic: versions, emit: versions_coreutils

    script:
        """
        cat ${vsgs} ${full_vsg_db} > concatenated_vsgs.fasta
        """
}
