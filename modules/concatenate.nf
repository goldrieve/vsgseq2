process CONCATENATE_VSGS {
    conda "$projectDir/vsgseq2.yml" 
    publishDir "${params.outdir}/VSGs", mode:'copy'
    input:
    path (vsgome)
    val (full_vsg_db)

    output:
    path "concatenated_vsgs.fasta"

    script:
    """
    cat ${vsgome} ${full_vsg_db} > concatenated_vsgs.fasta
    """
}
