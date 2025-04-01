process BLAST {
    conda "$projectDir/vsgseq2.yml"
    container 'goldrieve/vsgseq2:latest'
    publishDir "${params.outdir}/VSGs", mode:'copy', pattern: '*_VSGs.fasta'
    
    input:
    path assemblies
    val vsg_db
    val notvsg_db

    output:
    path "*.xml", emit: vsgblast
    path "*_nonVSG.xml", emit: notvsgblast
    path "*_VSGs.fasta", emit: vsgs

    script:
    // Handle filename based on mode
    def basename = assemblies.simpleName.replace('_cdhit', '')
    basename = (params.mode == 'new_full' || params.mode == 'new_analyse') ? basename : basename + '_cdhit'
    """
    blastn -db ${vsg_db} -query ${assemblies} -outfmt 5 -out ${basename}.xml
    blastn -db ${notvsg_db} -query ${assemblies} -outfmt 5 -out ${basename}_nonVSG.xml
    python $projectDir/bin/process_vsgs.py ${basename}
    """
}