process BLAST {
    publishDir "${params.outdir}/VSGs", mode:'copy', pattern: '*_VSGs.fasta'
    input:
    path (assemblies)
    val (vsg_db)
    val (notvsg_db)

    output:
    path "${assemblies.baseName.replace("_cdhit","")}.xml", emit: vsgblast
    path "${assemblies.baseName.replace("_cdhit","")}_nonVSG.xml", emit: notvsgblast
    path "${assemblies.baseName.replace("_cdhit","")}_VSGs.fasta", emit: vsgs

    script:
    """
    blastn -db $vsg_db -query ${assemblies} -outfmt 5 -out ${assemblies.baseName.replace("_cdhit","")}.xml
    blastn -db $notvsg_db -query ${assemblies} -outfmt 5 -out ${assemblies.baseName.replace("_cdhit","")}_nonVSG.xml
    python $projectDir/bin/process_vsgs.py ${assemblies.baseName.replace("_cdhit","")}
    """
}
