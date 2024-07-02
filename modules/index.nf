process INDEX {
    input:
    path (vsgs)
    val cores

    output:
    path "salmon_index"


    script:
    """
    cat $vsgs $projectDir/data/genome/TriTrypDB-67_TbruceiEATRO1125_Genome.fasta > gentrome.fasta
    salmon index --threads $cores -t gentrome.fasta -d $projectDir/data/genome/TriTrypDB-67_TbruceiEATRO1125_Genome_decoys.txt -i salmon_index --gencode
    """
}
