process INDIVIDUAL_CDHIT {
    input:
    path (assemblies)
    val cdhitid

    output:
    path "${assemblies.baseName.replace("_ORF","")}_cdhit.fasta"

    script:
    """
    cd-hit-est -i ${assemblies} -o ${assemblies.baseName.replace("_ORF","")}_cdhit.fasta -d 0 -c ${cdhitid} -n 8 -G 1 -g 1 -s 0.0 -aL 0.0 -M 50
    """
}