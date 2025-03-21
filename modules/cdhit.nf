process INDIVIDUAL_CDHIT {
    conda "$projectDir/vsgseq2.yml"
    container 'goldrieve/vsgseq2:latest'
    input:
    path (assemblies)
    val cdhit_id
    val cdhit_as

    output:
    path "${assemblies.baseName.replace("_ORF","")}_cdhit.fasta"

    script:
    """
    cd-hit-est -i ${assemblies} -o ${assemblies.baseName.replace("_ORF","")}_cdhit.fasta -d 0 -c ${cdhit_id} -n 8 -G 0 -g 1 -s 0.0 -aS ${cdhit_as} -M 50
    """
}