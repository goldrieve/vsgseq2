process BLAST {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/blast_biopython_python:118ae71af8a11bc3"

    input:
        tuple val(meta), path(orfs)
        tuple val(vsg_db), path(vsg_db_files, stageAs: 'vsg_db/*')
        tuple val(notvsg_db), path(notvsg_db_files, stageAs: 'notvsg_db/*')

    output:
        tuple val(meta), path("${meta.id}_ORF.xml"), emit: vsgblast
        tuple val(meta), path("${meta.id}_ORF_nonVSG.xml"), emit: notvsgblast
        tuple val(meta), path("${meta.id}_ORF_VSGs.fasta"), emit: vsgs
        tuple val("${task.process}"), val('blast'), eval("blastn -version | sed -n 's/^blastn: //p'"), topic: versions, emit: versions_blast
        tuple val("${task.process}"), val('biopython'), eval("python -c 'import Bio; print(Bio.__version__)'"), topic: versions, emit: versions_biopython

    script:
        def prefix = "${meta.id}_ORF"
        """
        if [ "${orfs}" != "${prefix}.fasta" ]; then
            ln -s ${orfs} ${prefix}.fasta
        fi

        blastn -db vsg_db/${vsg_db} -query ${prefix}.fasta -outfmt 5 -out ${prefix}.xml
        blastn -db notvsg_db/${notvsg_db} -query ${prefix}.fasta -outfmt 5 -out ${prefix}_nonVSG.xml
        process_vsgs.py ${prefix}
        """
}
