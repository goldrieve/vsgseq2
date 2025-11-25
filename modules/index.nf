process INDEX {
    
    input:
        path vsgs
        val cores
        path genome
        path decoys

    output:
        path "salmon_index"

    script:
        """
        cat $vsgs $genome > gentrome.fasta
        salmon index --threads $cores -t gentrome.fasta -d $decoys -i salmon_index --gencode
        """
}
