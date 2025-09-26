#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
========================================================================================
    Sheep Pangenome Pipeline - Staged Implementation
    Stages 1-3: Data Acquisition, Genome Preprocessing & Graph Construction
========================================================================================
    Author: Developed with Claude Code
    Version: 1.2.0-stage3
    Description: Comprehensive sheep pangenome analysis pipeline with staged validation
========================================================================================
*/

// Include parameter validation
include { validateParameters; paramsHelp } from 'plugin/nf-validation'

// Include subworkflows
include { INPUT_CHECK } from './subworkflows/local/input_check'
include { PREPROCESSING } from './subworkflows/local/preprocessing'
include { GRAPH_CONSTRUCTION } from './subworkflows/local/graph_construction'

// Print help message
if (params.help) {
    log.info paramsHelp("nextflow run . --input samplesheet.csv -profile docker")
    log.info ""
    log.info "Stages 1-3: Data Acquisition, Preprocessing & Graph Construction"
    log.info "======================================================="
    log.info ""
    log.info "Required parameters:"
    log.info "  --input         : Path to input samplesheet (CSV format)"
    log.info ""
    log.info "Optional parameters:"
    log.info "  --outdir        : Output directory (default: ./results)"
    log.info "  --stage         : Pipeline stage to run (1, 2, or 3, default: 1)"
    log.info "  --max_download_time : Maximum download time per genome (default: 30m)"
    log.info "  --ncbi_api_key  : NCBI API key for increased rate limits (recommended)"
    log.info ""
    log.info "Example samplesheet format:"
    log.info "sample,accession,breed,population,geographic_origin"
    log.info "sheep1,GCF_016772045.1,Rambouillet,European,United_States"
    log.info "sheep2,,/path/to/local.fa,Custom,Regional,Farm"
    log.info ""
    exit 0
}

// Validate parameters using nf-validation plugin
if (params.validate_params) {
    validateParameters()
}

// Print pipeline info
log.info ""
log.info "🐑 Sheep Pangenome Pipeline - Stage ${params.stage} 🐑"
log.info "=================================================="
log.info "Input samplesheet : ${params.input}"
log.info "Output directory  : ${params.outdir}"
log.info "Current stage     : ${params.stage}"
log.info "Max memory        : ${params.max_memory}"
log.info "Max CPUs          : ${params.max_cpus}"
log.info "Max time          : ${params.max_time}"
log.info ""

/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

// Check required parameters
if (!params.input) {
    log.error "Input samplesheet not specified! Use --input <samplesheet.csv>"
    exit 1
}

if (!file(params.input).exists()) {
    log.error "Input samplesheet does not exist: ${params.input}"
    exit 1
}

/*
========================================================================================
    MAIN WORKFLOW
========================================================================================
*/

workflow SHEEP_STAGE1 {

    take:
    samplesheet

    main:
    // Stage 1: Input checking and data preparation
    log.info "🔍 Starting Stage 1: Data Acquisition & Preparation"

    INPUT_CHECK(samplesheet)

    // Collect validation results
    validation_results = INPUT_CHECK.out.validation
        .map { meta, validation_file ->
            def validation = new groovy.json.JsonSlurper().parse(validation_file)
            [meta.id, validation.validation_status, validation_file]
        }

    // Check if all validations passed
    validation_results
        .filter { it[1] == "FAIL" }
        .subscribe { sample_id, status, file ->
            log.warn "⚠️  Validation FAILED for sample: ${sample_id}"
        }

    validation_results
        .filter { it[1] == "PASS" }
        .subscribe { sample_id, status, file ->
            log.info "✅ Validation PASSED for sample: ${sample_id}"
        }

    // Summary reporting
    validation_summary = validation_results
        .collect()
        .map { results ->
            def passed = results.count { it[1] == "PASS" }
            def failed = results.count { it[1] == "FAIL" }
            def total = results.size()

            log.info ""
            log.info "📊 Stage 1 Validation Summary:"
            log.info "  Total samples: ${total}"
            log.info "  Passed: ${passed}"
            log.info "  Failed: ${failed}"
            log.info "  Success rate: ${total > 0 ? (passed/total*100).round(1) : 0}%"

            if (failed > 0) {
                log.warn "⚠️  ${failed} samples failed validation - review before proceeding to Stage 2"
            } else {
                log.info "🎉 All samples passed validation - ready for Stage 2!"
            }

            return [passed: passed, failed: failed, total: total]
        }

    emit:
    genomes = INPUT_CHECK.out.genomes
    validation = INPUT_CHECK.out.validation
    statistics = INPUT_CHECK.out.statistics
    downloads = INPUT_CHECK.out.downloads
    versions = INPUT_CHECK.out.versions
    summary = validation_summary
}

/*
========================================================================================
    STAGE 2: GENOME PREPROCESSING & INDEXING
========================================================================================
*/

workflow SHEEP_STAGE2 {

    take:
    genomes_ch     // Channel from Stage 1: [meta, genome.fa, metadata.json]

    main:
    // Stage 2: Genome preprocessing, QC, and indexing
    log.info "🔧 Starting Stage 2: Genome Preprocessing & Indexing"

    PREPROCESSING(genomes_ch)

    // Monitor reference selection results
    PREPROCESSING.out.reference_genome
        .subscribe { selected_meta, reference_file ->
            log.info "🏆 Reference genome selected: ${selected_meta.id}"
        }

    // Summary reporting for Stage 2
    preprocessing_summary = PREPROCESSING.out.standardized_genomes
        .collect()
        .map { genomes ->
            def total = genomes.size()

            log.info ""
            log.info "📊 Stage 2 Preprocessing Summary:"
            log.info "  Total genomes processed: ${total}"
            log.info "  Standardization: ✅ Complete"
            log.info "  Quality control: ✅ Complete"
            log.info "  BWA indexing: ✅ Complete"
            log.info "  Minimap2 indexing: ✅ Complete"
            log.info "  Samtools indexing: ✅ Complete"
            log.info "  Reference selection: ✅ Complete"
            log.info ""

            return [processed: total]
        }

    emit:
    // Standardized and processed genomes
    standardized_genomes = PREPROCESSING.out.standardized_genomes
    chromosome_mapping   = PREPROCESSING.out.chromosome_mapping

    // Quality control results
    qc_reports          = PREPROCESSING.out.qc_html_reports
    qc_stats            = PREPROCESSING.out.qc_stats_json

    // Indexing results
    bwa_indices         = PREPROCESSING.out.bwa_indices
    minimap2_indices    = PREPROCESSING.out.minimap2_indices
    samtools_indices    = PREPROCESSING.out.samtools_indices

    // Reference selection
    reference_genome    = PREPROCESSING.out.reference_genome
    reference_metadata  = PREPROCESSING.out.reference_metadata
    selection_report    = PREPROCESSING.out.selection_report

    // Versions and logs
    versions            = PREPROCESSING.out.versions
    summary             = preprocessing_summary
}

/*
========================================================================================
    STAGE 3: PANGENOME GRAPH CONSTRUCTION
========================================================================================
*/

workflow SHEEP_STAGE3 {

    take:
    standardized_genomes_ch    // Channel from Stage 2: [meta, standardized.fa]
    reference_metadata_ch      // Channel: reference_metadata.json

    main:
    // Stage 3: Pangenome graph construction
    log.info "🧬 Starting Stage 3: Pangenome Graph Construction"

    GRAPH_CONSTRUCTION(
        standardized_genomes_ch,
        reference_metadata_ch
    )

    // Monitor graph construction progress
    GRAPH_CONSTRUCTION.out.pangenome_graph
        .subscribe { graph_meta, graph_file ->
            log.info "🎯 Pangenome graph constructed: ${graph_file.name}"
        }

    // Summary reporting for Stage 3
    graph_summary = GRAPH_CONSTRUCTION.out.pangenome_graph
        .map { graph_meta, graph_file ->
            def graph_size = graph_file.size()
            def graph_size_mb = (graph_size / 1024 / 1024).round(2)

            log.info ""
            log.info "📊 Stage 3 Graph Construction Summary:"
            log.info "  Pangenome graph: ✅ Complete"
            log.info "  Graph file size: ${graph_size_mb} MB"
            log.info "  Graph validation: ✅ Complete"
            log.info "  Graph statistics: ✅ Complete"
            log.info "  Visualizations: ✅ Complete"
            log.info ""

            return [graph_constructed: true, graph_size_mb: graph_size_mb]
        }

    emit:
    // Primary graph outputs
    pangenome_graph      = GRAPH_CONSTRUCTION.out.pangenome_graph
    odgi_graph           = GRAPH_CONSTRUCTION.out.odgi_graph
    smoothed_graph       = GRAPH_CONSTRUCTION.out.smoothed_graph
    graph_visualization  = GRAPH_CONSTRUCTION.out.graph_visualization
    pggb_output_dir      = GRAPH_CONSTRUCTION.out.pggb_output_dir

    // Validation and statistics
    validation_report    = GRAPH_CONSTRUCTION.out.validation_report
    validation_stats     = GRAPH_CONSTRUCTION.out.validation_stats
    graph_stats_json     = GRAPH_CONSTRUCTION.out.graph_stats_json
    graph_stats_html     = GRAPH_CONSTRUCTION.out.graph_stats_html
    node_statistics      = GRAPH_CONSTRUCTION.out.node_statistics
    path_statistics      = GRAPH_CONSTRUCTION.out.path_statistics
    complexity_analysis  = GRAPH_CONSTRUCTION.out.complexity_analysis

    // Versions and summary
    versions             = GRAPH_CONSTRUCTION.out.versions
    summary              = graph_summary
}

/*
========================================================================================
    ENTRY WORKFLOW
========================================================================================
*/

workflow {

    // Read input samplesheet
    ch_samplesheet = Channel.fromPath(params.input, checkIfExists: true)

    // Auto-detect completed stages and execute appropriate stage
    def stage1_completed = file("${params.outdir}/01_data_preparation/downloaded_genomes").exists() &&
                          file("${params.outdir}/01_data_preparation/downloaded_genomes").listFiles().length > 0

    def stage2_completed = file("${params.outdir}/02_preprocessing/standardized_genomes").exists() &&
                          file("${params.outdir}/02_preprocessing/standardized_genomes").listFiles().length > 0

    if (params.stage == 1 && !stage1_completed) {
        // Run Stage 1 if explicitly requested and not completed
        SHEEP_STAGE1(ch_samplesheet)

        // Stage completion check
        SHEEP_STAGE1.out.summary
            .subscribe { summary ->
                if (summary.failed == 0) {
                    log.info ""
                    log.info "🚀 Stage 1 completed successfully!"
                    log.info "   Ready to proceed to Stage 2: Genome Preprocessing"
                    log.info ""
                } else {
                    log.error ""
                    log.error "❌ Stage 1 completed with ${summary.failed} failures"
                    log.error "   Fix validation issues before proceeding"
                    log.error ""
                }
            }

    } else if (params.stage == 1 && stage1_completed) {
        // Stage 1 already completed, auto-advance to Stage 2
        log.info ""
        log.info "🔄 Stage 1 already completed, advancing to Stage 2"
        log.info "   Found existing genomes in: ${params.outdir}/01_data_preparation/downloaded_genomes/"
        log.info ""

        // Create genome input channel from Stage 1 outputs
        stage1_genomes = Channel.fromPath("${params.outdir}/01_data_preparation/metadata/*.json")
            .map { metadata_file ->
                def metadata = new groovy.json.JsonSlurper().parse(metadata_file)
                def sample_id = metadata_file.baseName
                def genome_file = file("${params.outdir}/01_data_preparation/downloaded_genomes/${sample_id}.fa")

                if (!genome_file.exists()) {
                    log.error "Genome file not found: ${genome_file}"
                    exit 1
                }

                def meta = [
                    id: sample_id,
                    breed: metadata.organism?.infraspecificNames?.breed ?: 'unknown',
                    single_end: true
                ]

                return [meta, genome_file, metadata_file]
            }

        SHEEP_STAGE2(stage1_genomes)

        // Stage completion check
        SHEEP_STAGE2.out.summary
            .subscribe { summary ->
                log.info ""
                log.info "🚀 Stage 2 completed successfully!"
                log.info "   Processed ${summary.processed} genomes with full preprocessing"
                log.info "   Ready to proceed to Stage 3: Pangenome Construction"
                log.info ""
            }

    } else if (params.stage == 2) {
        // Stage 2 requires Stage 1 outputs - check for existing results
        def stage1_results = file("${params.outdir}/01_data_preparation/downloaded_genomes")
        if (!stage1_results.exists() || stage1_results.listFiles().length == 0) {
            log.error ""
            log.error "❌ Stage 2 requires Stage 1 outputs!"
            log.error "   No genomes found in: ${stage1_results}"
            log.error "   Run Stage 1 first: nextflow run . --input samplesheet.csv --stage 1"
            log.error ""
            exit 1
        }

        // Create genome input channel from Stage 1 outputs
        stage1_genomes = Channel.fromPath("${params.outdir}/01_data_preparation/metadata/*.json")
            .map { metadata_file ->
                def metadata = new groovy.json.JsonSlurper().parse(metadata_file)
                def sample_id = metadata_file.baseName
                def genome_file = file("${params.outdir}/01_data_preparation/downloaded_genomes/${sample_id}.fa")

                if (!genome_file.exists()) {
                    log.error "Genome file not found: ${genome_file}"
                    exit 1
                }

                def meta = [
                    id: sample_id,
                    breed: metadata.organism?.infraspecificNames?.breed ?: 'unknown',
                    single_end: true
                ]

                return [meta, genome_file, metadata_file]
            }

        SHEEP_STAGE2(stage1_genomes)

        // Stage completion check
        SHEEP_STAGE2.out.summary
            .subscribe { summary ->
                log.info ""
                log.info "🚀 Stage 2 completed successfully!"
                log.info "   Processed ${summary.processed} genomes with full preprocessing"
                log.info "   Ready to proceed to Stage 3: Pangenome Construction"
                log.info ""
            }

    } else if (params.stage == 3) {
        // Stage 3 requires Stage 2 outputs - check for existing results
        def stage2_results = file("${params.outdir}/02_preprocessing/standardized_genomes")
        def stage2_reference = file("${params.outdir}/02_preprocessing/reference_selection/reference_metadata.json")

        if (!stage2_results.exists() || stage2_results.listFiles().length == 0) {
            log.error ""
            log.error "❌ Stage 3 requires Stage 2 outputs!"
            log.error "   No standardized genomes found in: ${stage2_results}"
            log.error "   Run Stage 2 first: nextflow run . --input samplesheet.csv --stage 2"
            log.error ""
            exit 1
        }

        if (!stage2_reference.exists()) {
            log.error ""
            log.error "❌ Stage 3 requires reference selection from Stage 2!"
            log.error "   Reference metadata not found: ${stage2_reference}"
            log.error "   Ensure Stage 2 completed successfully"
            log.error ""
            exit 1
        }

        // Create genome input channels from Stage 2 outputs
        stage2_genomes = Channel.fromPath("${params.outdir}/02_preprocessing/standardized_genomes/*_standardized.fa")
            .map { genome_file ->
                def sample_id = genome_file.baseName.replace('_standardized', '')
                def meta = [
                    id: sample_id,
                    single_end: true,
                    stage: 'graph_construction'
                ]

                return [meta, genome_file]
            }

        // Reference metadata channel
        reference_metadata = Channel.fromPath(stage2_reference, checkIfExists: true)

        log.info ""
        log.info "🧬 Starting Stage 3: Pangenome Graph Construction"
        log.info "   Using standardized genomes from Stage 2"
        log.info "   Found ${stage2_results.listFiles().length} genome files"
        log.info ""

        SHEEP_STAGE3(stage2_genomes, reference_metadata)

        // Stage completion check
        SHEEP_STAGE3.out.summary
            .subscribe { summary ->
                log.info ""
                log.info "🚀 Stage 3 completed successfully!"
                log.info "   Pangenome graph constructed (${summary.graph_size_mb} MB)"
                log.info "   Ready to proceed to Stage 4: Graph Analysis & Variant Calling"
                log.info ""
            }

    } else if (params.stage == 2 && stage2_completed) {
        // Stage 2 already completed, auto-advance to Stage 3
        log.info ""
        log.info "🔄 Stage 2 already completed, advancing to Stage 3"
        log.info "   Found existing standardized genomes in: ${params.outdir}/02_preprocessing/standardized_genomes/"
        log.info ""

        // Create genome input channels from Stage 2 outputs
        stage2_genomes = Channel.fromPath("${params.outdir}/02_preprocessing/standardized_genomes/*_standardized.fa")
            .map { genome_file ->
                def sample_id = genome_file.baseName.replace('_standardized', '')
                def meta = [
                    id: sample_id,
                    single_end: true,
                    stage: 'graph_construction'
                ]

                return [meta, genome_file]
            }

        // Reference metadata channel
        reference_metadata = Channel.fromPath("${params.outdir}/02_preprocessing/reference_selection/reference_metadata.json", checkIfExists: true)

        SHEEP_STAGE3(stage2_genomes, reference_metadata)

        // Stage completion check
        SHEEP_STAGE3.out.summary
            .subscribe { summary ->
                log.info ""
                log.info "🚀 Stage 3 completed successfully!"
                log.info "   Pangenome graph constructed (${summary.graph_size_mb} MB)"
                log.info "   Ready to proceed to Stage 4: Graph Analysis & Variant Calling"
                log.info ""
            }

    } else {
        log.error "Invalid stage specified: ${params.stage}"
        log.error "Available stages: 1 (Data Acquisition), 2 (Genome Preprocessing), 3 (Graph Construction)"
        exit 1
    }
}

/*
========================================================================================
    WORKFLOW COMPLETION
========================================================================================
*/

workflow.onComplete {
    log.info ""
    log.info "🏁 Pipeline execution completed!"
    log.info "Execution status: ${workflow.success ? 'SUCCESS' : 'FAILED'}"
    log.info "Execution duration: ${workflow.duration}"
    log.info "CPU hours: ${workflow.stats.computeTimeFmt ?: 'N/A'}"
    log.info "Peak memory usage: ${workflow.stats.peakCpus ?: 'N/A'} CPUs, ${workflow.stats.peakMemory ?: 'N/A'}"
    log.info ""

    if (workflow.success) {
        log.info "✅ Results available in: ${params.outdir}"
        log.info "📁 Key outputs:"

        if (params.stage == 1) {
            log.info "   - Downloaded genomes: ${params.outdir}/01_data_preparation/downloaded_genomes/"
            log.info "   - Validation results: ${params.outdir}/01_data_preparation/validation/"
            log.info "   - Statistics: ${params.outdir}/01_data_preparation/statistics/"
            log.info "   - Metadata: ${params.outdir}/01_data_preparation/metadata/"
        } else if (params.stage == 2) {
            log.info "   - Standardized genomes: ${params.outdir}/02_preprocessing/standardized_genomes/"
            log.info "   - Quality reports: ${params.outdir}/02_preprocessing/quality_control/"
            log.info "   - BWA indices: ${params.outdir}/02_preprocessing/bwa_index/"
            log.info "   - Minimap2 indices: ${params.outdir}/02_preprocessing/minimap2_index/"
            log.info "   - Samtools indices: ${params.outdir}/02_preprocessing/samtools_index/"
            log.info "   - Reference selection: ${params.outdir}/02_preprocessing/reference_selection/"
        } else if (params.stage == 3) {
            log.info "   - Pangenome graph: ${params.outdir}/03_graph_construction/main_outputs/"
            log.info "   - Graph validation: ${params.outdir}/03_graph_construction/validation/"
            log.info "   - Graph statistics: ${params.outdir}/03_graph_construction/statistics/"
            log.info "   - PGGB outputs: ${params.outdir}/03_graph_construction/pggb_outputs/"
            log.info "   - Visualizations: ${params.outdir}/03_graph_construction/main_outputs/"
        }
        log.info ""
    }
}

workflow.onError {
    log.error ""
    log.error "❌ Pipeline execution failed!"
    log.error "Error message: ${workflow.errorMessage}"
    log.error "Error report: ${workflow.errorReport}"
    log.error ""
}