#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
========================================================================================
    Sheep Pangenome Pipeline - Staged Implementation
    Stage 1: Data Acquisition & Preparation
========================================================================================
    Author: Developed with Claude Code
    Version: 1.0.0-stage1
    Description: Comprehensive sheep pangenome analysis pipeline with staged validation
========================================================================================
*/

// Include parameter validation
include { validateParameters; paramsHelp } from 'plugin/nf-validation'

// Include subworkflows
include { INPUT_CHECK } from './subworkflows/local/input_check'

// Print help message
if (params.help) {
    log.info paramsHelp("nextflow run . --input samplesheet.csv -profile docker")
    log.info ""
    log.info "Stage 1: Data Acquisition & Preparation"
    log.info "======================================="
    log.info ""
    log.info "Required parameters:"
    log.info "  --input         : Path to input samplesheet (CSV format)"
    log.info ""
    log.info "Optional parameters:"
    log.info "  --outdir        : Output directory (default: ./results)"
    log.info "  --stage         : Pipeline stage to run (default: 1)"
    log.info "  --max_download_time : Maximum download time per genome (default: 30m)"
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
    ENTRY WORKFLOW
========================================================================================
*/

workflow {

    // Read input samplesheet
    ch_samplesheet = Channel.fromPath(params.input, checkIfExists: true)

    // Execute current stage
    if (params.stage == 1) {
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
    } else {
        log.error "Only Stage 1 is currently implemented. Use --stage 1"
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
        log.info "   - Downloaded genomes: ${params.outdir}/01_data_preparation/downloaded_genomes/"
        log.info "   - Validation results: ${params.outdir}/01_data_preparation/validation/"
        log.info "   - Statistics: ${params.outdir}/01_data_preparation/statistics/"
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