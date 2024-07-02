process CONCATENATED_CDHIT {
    publishDir "${params.outdir}/VSGs", mode:'copy'

    input:
    path (vsgs)
    val cdhitid

    output:
    path "VSGome.fasta"

    script:
    """
    cd-hit-est -i ${vsgs} -o VSGome.fasta -d 0 -c ${cdhitid} -n 8 -G 1 -g 1 -s 0.0 -aL 0.0 -M 5000
    """
}
