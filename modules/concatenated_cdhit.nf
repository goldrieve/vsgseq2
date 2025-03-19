process CONCATENATED_CDHIT {
    conda "$projectDir/vsgseq2.yml"
    container 'goldrieve/vsgseq2:latest'
    publishDir "${params.outdir}/VSGs", mode:'copy'

    input:
    path (vsgs)
    val cdhit_id
    val cdhit_as

    output:
    path "VSGome.fasta"

    script:
    """
    cd-hit-est -i ${vsgs} -o VSGome.fasta -d 0 -c ${cdhit_id} -n 8 -G 0 -g 1 -s 0.0 -aS ${cdhit_as}  -M 50
    """
}
