/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'

include { TRIM               } from '../modules/local/trim/main'
include { ASSEMBLE           } from '../modules/local/assemble/main'
include { ORF                } from '../modules/local/orf/main'
include { BLAST              } from '../modules/local/blast/main'
include { CONCATENATE_VSGS   } from '../modules/local/concatenate_vsgs/main'
include { CONCATENATED_CDHIT } from '../modules/local/concatenated_cdhit/main'
include { INDEX              } from '../modules/local/index/main'
include { QUANTIFY           } from '../modules/local/quantify/main'
include { MULTIQC            } from '../modules/local/multiqc/main'
include { SUMMARISE          } from '../modules/local/summarise/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow VSGSEQ2 {

    take:
    ch_samplesheet // channel: [ val(meta), [ reads ] ] one per samplesheet row
    outdir

    main:

    //
    // MODULE: Trim reads with Trimmomatic
    //
    TRIM (
        ch_samplesheet,
        params.cores
    )

    if (params.mode == 'full') {
        //
        // MODULE: Trinity assembly
        //
        ASSEMBLE (
            TRIM.out.trimmed_reads,
            params.cores,
            params.trinitymem
        )
        ch_assemblies = ASSEMBLE.out.trinity_fasta
    } else {
        ch_assemblies = channel
            .fromPath(params.assemblies, checkIfExists: true)
            .map { assembly ->
                def id = assembly.baseName.replace('_trinity.Trinity', '')
                if (!(id ==~ /^[A-Za-z0-9][A-Za-z0-9_-]*$/) || !(assembly.name ==~ /^[A-Za-z0-9._-]+$/)) {
                    error("Invalid assembly file name '${assembly.name}': sample names must start with a letter or digit and contain only letters, digits, underscores or hyphens, e.g. '<sample>_trinity.Trinity.fasta'")
                }
                [ [ id: id ], assembly ]
            }
    }

    //
    // MODULE: Predict ORFs with TransDecoder
    //
    ORF (
        ch_assemblies,
        params.cdslength,
        params.partial
    )

    //
    // MODULE: BLAST ORFs against the VSG and non-VSG databases
    //
    ch_vsg_db = channel.fromPath("${params.vsg_db}*", checkIfExists: true)
        .collect()
        .map { db_files -> [ file(params.vsg_db).name, db_files ] }
    ch_notvsg_db = channel.fromPath("${params.notvsg_db}*", checkIfExists: true)
        .collect()
        .map { db_files -> [ file(params.notvsg_db).name, db_files ] }

    BLAST (
        ORF.out.fasta,
        ch_vsg_db,
        ch_notvsg_db
    )

    ch_vsgs = BLAST.out.vsgs
        .map { _meta, vsgs -> vsgs }
        .collect(sort: { a, b -> a.name <=> b.name })

    //
    // MODULE: Concatenate VSGs from all samples, optionally with an existing VSG database
    //
    CONCATENATE_VSGS (
        ch_vsgs,
        params.full_vsg_db ? file(params.full_vsg_db, checkIfExists: true) : []
    )

    //
    // MODULE: Cluster VSGs with cd-hit-est
    //
    CONCATENATED_CDHIT (
        CONCATENATE_VSGS.out.fasta,
        params.cdhit_id,
        params.cdhit_as
    )

    //
    // MODULE: Salmon index of VSGs with the genome as decoys
    //
    INDEX (
        CONCATENATE_VSGS.out.fasta,
        params.cores,
        file(params.genome, checkIfExists: true),
        file(params.decoys, checkIfExists: true)
    )

    //
    // MODULE: Salmon quantification
    //
    QUANTIFY (
        INDEX.out.index.collect(),
        params.cores,
        TRIM.out.trimmed_reads
    )

    ch_quants = QUANTIFY.out.quants
        .map { _meta, quant -> quant }
        .collect(sort: { a, b -> a.name <=> b.name })

    //
    // MODULE: MultiQC report of Salmon quantification
    //
    MULTIQC (
        ch_quants
    )

    //
    // MODULE: Summarise expression, VSG counts and clusters
    //
    SUMMARISE (
        ch_quants,
        params.threshold,
        ch_vsgs,
        CONCATENATED_CDHIT.out.clstr,
        CONCATENATE_VSGS.out.fasta
    )

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(topic_versions.versions_file)
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name:  'vsgseq2_software_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    emit:
    champ_vsgs = SUMMARISE.out.champ_vsgs // channel: [ path(champion_vsgs.fasta) ]
    versions   = ch_collated_versions     // channel: [ path(vsgseq2_software_versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
