/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { BAM_TO_SAF             } from '../subworkflows/local/bam_to_saf'
include { BAM_TO_STRUCTURE       } from '../subworkflows/local/bam_to_structure'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_glpopgen_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow GLPOPGEN {

    take:
    ch_samplesheet // channel: [ val(meta), path(bam), path(bai) ] read in from --input; meta = [ id, population ]
    multiqc_config
    multiqc_logo
    multiqc_methods_description
    outdir

    main:

    def ch_versions = channel.empty()
    def ch_multiqc_files = channel.empty()

    //
    // Reference and (optional) ancestral genomes: reused by every population
    //
    def ch_reference = channel.value([
        file(params.fasta, checkIfExists: true),
        file(params.fasta_fai ?: "${params.fasta}.fai", checkIfExists: true)
    ])
    def ch_ancestral = params.ancestral_fasta
        ? channel.value([
            file(params.ancestral_fasta, checkIfExists: true),
            file(params.ancestral_fai ?: "${params.ancestral_fasta}.fai", checkIfExists: true)
        ])
        : channel.value([ [], [] ])

    //
    // All samples together, sorted by sample ID so the order (and -resume cache) is stable: [ meta, [ bams ], [ bais ] ]
    //
    def ch_all_bams = ch_samplesheet
        .map { meta, bam, bai -> [ meta.id, bam, bai ] }
        .toSortedList { a, b -> a[0] <=> b[0] }
        .map { samples ->
            [ [ id: 'all_samples', samples: samples.collect { s -> s[0] } ], samples.collect { s -> s[1] }, samples.collect { s -> s[2] } ]
        }

    //
    // SUBWORKFLOW: Population structure across all samples (genotype likelihoods -> PCA, admixture)
    //
    BAM_TO_STRUCTURE (
        ch_all_bams,
        ch_reference,
        params.admix_k_min,
        params.admix_k_max,
        params.admix_reps
    )

    //
    // Sample IDs of the all-samples Beagle columns, one per line (first line = Ind0)
    //
    BAM_TO_STRUCTURE.out.beagle
        .collectFile(storeDir: "${outdir}/angsd/gl_all") { meta, _beagle ->
            [ "${meta.id}.samples.txt", meta.samples.join('\n') + '\n' ]
        }

    //
    // Group samples by population: [ meta, [ bams ], [ bais ] ]
    //
    def ch_population_bams = ch_samplesheet
        .map { meta, bam, bai -> [ meta.population, meta.id, bam, bai ] }
        .groupTuple()
        .map { population, ids, bams, bais ->
            [ [ id: population, population: population, samples: ids ], bams, bais ]
        }

    //
    // SUBWORKFLOW: Site allele frequency likelihoods per population with ANGSD
    //
    BAM_TO_SAF (
        ch_population_bams,
        BAM_TO_STRUCTURE.out.beagle,
        ch_reference,
        ch_ancestral,
        params.gl_model,
        params.saf_method
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

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name: 'nf_core_'  +  'glpopgen_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    //
    // MODULE: MultiQC
    //
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    def ch_summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def ch_workflow_summary = channel.value(paramsSummaryMultiqc(ch_summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    def ch_multiqc_custom_methods_description = multiqc_methods_description
        ? file(multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    def ch_methods_description = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))
    MULTIQC(
        ch_multiqc_files.flatten().collect().map { files ->
            [
                [id: 'glpopgen'],
                files,
                multiqc_config
                    ? file(multiqc_config, checkIfExists: true)
                    : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
                multiqc_logo ? file(multiqc_logo, checkIfExists: true) : [],
                [],
                [],
            ]
        }
    )
    emit:multiqc_report = MULTIQC.out.report.map { _meta, report -> [report] }.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
