/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { EXTRACT                } from '../modules/local/extract'
include { MAKE_FRAGMENT          } from '../modules/local/snapatac2/make_fragment'
include { SNAPATAC2              } from '../modules/local/snapatac2/snapatac2'
include { MULTIQC                } from '../modules/local/multiqc_sgr'

include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { BWA_INDEX              } from '../modules/nf-core/bwa/index/main'
include { BWA_MEM                } from '../modules/nf-core/bwa/mem/main'
include { QUALIMAP_BAMQC         } from '../modules/nf-core/qualimap/bamqc/main'

include { paramsSummaryMap       } from 'plugin/nf-validation'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_scatac_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow scatac {

    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
    ch_fasta = [ [], params.fasta ]
    ch_gtf = [ [], params.gtf ]


    // fastqc
    if (params.run_fastqc) {
        FASTQC (
            ch_samplesheet
        )
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
        ch_versions = ch_versions.mix(FASTQC.out.versions.first())
    }

    EXTRACT (
        ch_samplesheet,
        "${projectDir}/assets/",
        params.protocol,
    )
    ch_versions = ch_versions.mix(EXTRACT.out.versions.first())
    ch_multiqc_files = ch_multiqc_files.mix(EXTRACT.out.json.collect{it[1]})

    if (params.bwa_index) {
        bwa_index = [ [], params.bwa_index ]
    } else {
        BWA_INDEX (
            ch_fasta,
        )
        bwa_index = BWA_INDEX.out.index
    }

    BWA_MEM (
        EXTRACT.out.out_reads,
        bwa_index,
        ch_fasta,
        true,
    )

    QUALIMAP_BAMQC(
        BWA_MEM.out.bam,
        [],
    )
    ch_multiqc_files = ch_multiqc_files.mix(QUALIMAP_BAMQC.out.results.collect{it[1]})
    ch_versions = ch_versions.mix(QUALIMAP_BAMQC.out.versions.first())

    MAKE_FRAGMENT(
        BWA_MEM.out.bam
    )
    ch_versions = ch_versions.mix(MAKE_FRAGMENT.out.versions.first())

    // snapatac2 
    SNAPATAC2 (
        MAKE_FRAGMENT.out.fragment,
        ch_fasta,
        ch_gtf,
    )
    ch_multiqc_files = ch_multiqc_files.mix(SNAPATAC2.out.json.collect{it[1]})
    ch_versions = ch_versions.mix(SNAPATAC2.out.versions.first())

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'nf_core_pipeline_software_mqc_versions.yml', sort: true, newLine: true)
        .set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config                     = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config              = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
    ch_multiqc_logo                       = params.multiqc_logo ? Channel.fromPath(params.multiqc_logo, checkIfExists: true) : Channel.empty()
    summary_params                        = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary                   = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: false))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        "${projectDir}/multiqc_sgr/singleron_logo.png",
        "${projectDir}/multiqc_sgr/",
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
