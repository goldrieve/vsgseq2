process BLAST {
    publishDir "${params.outdir}/VSGs", mode:'copy', pattern: '*_VSGs.fasta'
    
    input:
        path assemblies
        path "db_files/*"

    output:
        path "*.xml", emit: vsgblast
        path "*_nonVSG.xml", emit: notvsgblast
        path "*_VSGs.fasta", emit: vsgs

    script:
        def basename = assemblies.simpleName.replace('_cdhit', '')
        basename = (params.mode == 'full' || params.mode == 'analyse') ? basename : basename + '_cdhit'
        
        // Extract the database basename from the params
        def vsg_db_name = file(params.vsg_db).name.replaceAll(/\.(fa|fasta)$/, '')
        def notvsg_db_name = file(params.notvsg_db).name.replaceAll(/\.(fa|fasta)$/, '')
        """
        blastn -db db_files/${vsg_db_name}.fa -query ${assemblies} -outfmt 5 -out ${basename}.xml
        blastn -db db_files/${notvsg_db_name}.fa -query ${assemblies} -outfmt 5 -out ${basename}_nonVSG.xml
        python ${params.scripts}process_vsgs.py ${basename}
        """
}