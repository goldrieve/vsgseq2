process CONCATENATE_VSGS {
    publishDir "${params.outdir}/VSGs", mode:'copy'
    input:
    path (vsgome)

    output:
    path "concatenated_vsgs.fasta"

    script:
    """
    cat ${vsgome} > concatenated_vsgs.fasta
    """
}
