process CONCATENATED_CDHIT {
    conda "$projectDir/vsgseq2.yml"
    container 'goldrieve/vsgseq2:latest'
    publishDir "${params.outdir}/VSGs/VSGome", mode:'copy'

    input:
    path (vsgs)
    val cdhit_id
    val cdhit_as

    output:
    path "VSGome.fasta"
    path "VSGome.fasta.clstr"

    script:
    """
    cd-hit-est -i ${vsgs} -o VSGome.fasta -d 0 -c ${cdhit_id} -n 8 -G 0 -g 1 -s 0.0 -aS ${cdhit_as} -t 4 -M 2000
    """
}
