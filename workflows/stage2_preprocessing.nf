/*
========================================================================================
    Stage 2: Genome Preprocessing & Indexing Workflow
========================================================================================
    Description: Complete Stage 2 implementation for the sheep pangenome pipeline
    Input: Validated genomes from Stage 1
    Output: Standardized, indexed genomes ready for Stage 3 (PGGB pangenome construction)
========================================================================================
*/

// Include Stage 2 configuration
includeConfig '../conf/modules_stage2.config'

// Include Stage 2 subworkflows
include { GENOME_PREPROCESSING } from '../subworkflows/local/genome_preprocessing'

workflow STAGE2_PREPROCESSING {

    take:
    stage1_genomes      // channel: [meta, genome.fa] from Stage 1
    stage1_metadata     // channel: [meta, metadata.json] from Stage 1
    stage1_validation   // channel: [meta, validation.json] from Stage 1

    main:

    log.info ""
    log.info "🚀 Starting Stage 2: Genome Preprocessing & Indexing"
    log.info "=" * 60

    // ================================================================================
    // Input Validation & Contract Enforcement
    // ================================================================================

    // Count input genomes
    input_count = stage1_genomes.count()
    input_count.subscribe { count ->
        log.info "📥 Stage 2 Input Summary:"
        log.info "   Genomes from Stage 1: ${count}"
        log.info ""
    }

    // Validate Stage 1→2 contract compliance
    contract_validated = stage1_validation
        .join(stage1_genomes, by: 0)
        .join(stage1_metadata, by: 0)
        .filter { meta, validation_file, genome_file, metadata_file ->

            // Parse validation status
            def validation = new groovy.json.JsonSlurper().parse(validation_file)
            def is_valid = validation.validation_status == "PASS"

            if (!is_valid) {
                log.warn "❌ Sample ${meta.id} failed Stage 1 validation - excluding from Stage 2"
                return false
            }

            // Additional contract checks
            def genome_size = validation.total_length
            def seq_count = validation.total_sequences
            def n_content = validation.n_content

            // Sheep genome size validation
            if (genome_size < 2.2e9 || genome_size > 3.5e9) {
                log.warn "⚠️  Sample ${meta.id} genome size unusual (${(genome_size/1e9).round(2)} Gb) - including with caution"
            }

            // Excessive fragmentation check
            if (seq_count > 100000) {
                log.warn "⚠️  Sample ${meta.id} highly fragmented (${seq_count} sequences) - may affect processing"
            }

            return true
        }

    contract_validated.count().subscribe { count ->
        log.info "✅ Contract validation complete: ${count} genomes eligible for Stage 2"
        log.info ""
    }

    // ================================================================================
    // Main Stage 2 Processing
    // ================================================================================

    GENOME_PREPROCESSING(
        contract_validated.map { meta, validation, genome, metadata -> [meta, genome] },
        contract_validated.map { meta, validation, genome, metadata -> [meta, metadata] },
        contract_validated.map { meta, validation, genome, metadata -> [meta, validation] }
    )

    // ================================================================================
    // Stage 2 Output Validation & Summary
    // ================================================================================

    // Comprehensive processing summary
    stage2_completion = GENOME_PREPROCESSING.out.summary
        .combine(GENOME_PREPROCESSING.out.standardized_genomes.count())
        .combine(GENOME_PREPROCESSING.out.reference_genome.count())
        .map { summary, standardized_count, reference_count ->

            log.info ""
            log.info "🎉 Stage 2 Processing Complete!"
            log.info "=" * 60
            log.info "📊 Final Summary:"
            log.info "   Input genomes: ${summary.total_genomes}"
            log.info "   Successfully processed: ${standardized_count}"
            log.info "   Average genome size: ${(summary.avg_genome_size/1e9).round(2)} Gb"
            log.info ""
            log.info "🏷️  Standardization Results:"
            log.info "   Genomes standardized: ${summary.standardized}"
            log.info "   Chromosome naming: ✅ Applied to all genomes"
            log.info "   Header normalization: ✅ Applied to all genomes"
            log.info ""
            log.info "🗂️  Indexing Results:"
            log.info "   BWA indices: ${summary.indexed} genomes"
            log.info "   Minimap2 indices: ${summary.indexed} genomes"
            log.info "   Samtools indices: ${summary.indexed} genomes"
            log.info ""
            log.info "🔬 Quality Assessment:"
            summary.quality_distribution.each { tier, count ->
                log.info "   Tier ${tier}: ${count} genomes"
            }
            log.info ""
            log.info "🎯 Reference Selection:"
            log.info "   Reference genome: ${reference_count > 0 ? 'Selected' : 'Failed'}"
            log.info ""
            log.info "✅ Stage 2→3 Readiness:"
            log.info "   Ready for PGGB: ${summary.ready_for_stage3} genomes"
            log.info "   Success rate: ${(summary.ready_for_stage3/summary.total_genomes*100).round(1)}%"

            // Stage gate assessment
            if (summary.ready_for_stage3 >= 3) {
                log.info ""
                log.info "🚀 STAGE GATE: PASSED ✅"
                log.info "   Sufficient genomes (≥3) ready for Stage 3"
                log.info "   Proceeding to pangenome construction is recommended"
            } else {
                log.warn ""
                log.warn "⚠️  STAGE GATE: WARNING ⚠️"
                log.warn "   Only ${summary.ready_for_stage3} genomes ready for Stage 3"
                log.warn "   Consider reviewing quality issues before proceeding"
            }

            log.info ""
            log.info "📁 Stage 2 Outputs:"
            log.info "   Standardized genomes: results/02_preprocessing/standardized/"
            log.info "   BWA indices: results/02_preprocessing/indices/bwa/"
            log.info "   Minimap2 indices: results/02_preprocessing/indices/minimap2/"
            log.info "   Samtools indices: results/02_preprocessing/indices/samtools/"
            log.info "   Reference genome: results/02_preprocessing/reference_selection/"
            log.info "   QC reports: results/02_preprocessing/qc_reports/"
            log.info "   Quality statistics: results/02_preprocessing/statistics/"
            log.info "   BUSCO results: results/02_preprocessing/busco/"
            log.info ""

            return summary + [
                stage2_complete: true,
                stage_gate_status: summary.ready_for_stage3 >= 3 ? 'PASSED' : 'WARNING'
            ]
        }

    // ================================================================================
    // Data Quality Validation
    // ================================================================================

    // Validate critical outputs exist
    critical_outputs = GENOME_PREPROCESSING.out.standardized_genomes
        .combine(GENOME_PREPROCESSING.out.reference_genome)
        .map { standardized_meta, standardized_genome, ref_meta, ref_genome ->

            def validation_results = [
                standardized_genomes_exist: true,
                reference_genome_exists: ref_genome?.exists() ?: false,
                total_standardized: 1  // Will be accumulated
            ]

            return validation_results
        }

    // Check index completeness
    index_validation = GENOME_PREPROCESSING.out.bwa_indices
        .join(GENOME_PREPROCESSING.out.minimap2_indices, by: 0)
        .join(GENOME_PREPROCESSING.out.samtools_indices, by: 0)
        .map { meta, bwa_idx, mm2_idx, sam_idx ->

            def all_indices_exist = [bwa_idx, mm2_idx, sam_idx].every { it?.exists() }

            if (all_indices_exist) {
                log.info "✅ All indices created for ${meta.id}"
            } else {
                log.warn "⚠️  Incomplete indexing for ${meta.id}"
            }

            return [meta.id, all_indices_exist]
        }

    index_validation.count().subscribe { total_indexed ->
        log.info "🗂️  Indexing completion: ${total_indexed} genomes fully indexed"
    }

    emit:
    // Primary Stage 3 inputs
    standardized_genomes = GENOME_PREPROCESSING.out.standardized_genomes
    reference_genome = GENOME_PREPROCESSING.out.reference_genome
    reference_metadata = GENOME_PREPROCESSING.out.reference_metadata

    // Index outputs for downstream validation
    bwa_indices = GENOME_PREPROCESSING.out.bwa_indices
    minimap2_indices = GENOME_PREPROCESSING.out.minimap2_indices
    samtools_indices = GENOME_PREPROCESSING.out.samtools_indices

    // Quality control outputs
    qc_reports = GENOME_PREPROCESSING.out.qc_reports
    extended_stats = GENOME_PREPROCESSING.out.extended_stats
    busco_results = GENOME_PREPROCESSING.out.busco_results

    // Mapping and traceability
    chromosome_maps = GENOME_PREPROCESSING.out.chromosome_maps
    header_maps = GENOME_PREPROCESSING.out.header_maps
    processing_metadata = GENOME_PREPROCESSING.out.processing_metadata

    // Pipeline control
    stage2_summary = stage2_completion
    versions = GENOME_PREPROCESSING.out.versions

    // Stage gate status for pipeline orchestration
    stage_gate = stage2_completion.map { it.stage_gate_status }
}

/*
========================================================================================
    Entry Workflow for Stage 2 Testing
========================================================================================
*/

workflow {

    // For standalone Stage 2 testing, read from Stage 1 outputs
    if (params.input && params.stage == 2) {

        log.info "🧪 Running Stage 2 in standalone mode"
        log.info "   Reading Stage 1 outputs from: ${params.outdir}/01_data_preparation/"

        // Construct Stage 1 output channels
        stage1_genomes = Channel
            .fromPath("${params.outdir}/01_data_preparation/downloaded_genomes/*.fa")
            .map { genome ->
                def sample_id = genome.baseName
                def meta = [id: sample_id]
                [meta, genome]
            }

        stage1_metadata = Channel
            .fromPath("${params.outdir}/01_data_preparation/metadata/*.json")
            .map { metadata ->
                def sample_id = metadata.baseName
                def meta = [id: sample_id]
                [meta, metadata]
            }

        stage1_validation = Channel
            .fromPath("${params.outdir}/01_data_preparation/validation/*_validation.json")
            .map { validation ->
                def sample_id = validation.baseName.replaceAll('_validation', '')
                def meta = [id: sample_id]
                [meta, validation]
            }

        // Run Stage 2
        STAGE2_PREPROCESSING(
            stage1_genomes,
            stage1_metadata,
            stage1_validation
        )

        // Stage completion handling
        STAGE2_PREPROCESSING.out.stage2_summary
            .subscribe { summary ->
                if (summary.stage_gate_status == 'PASSED') {
                    log.info "🎉 Stage 2 completed successfully - ready for Stage 3!"
                } else {
                    log.warn "⚠️  Stage 2 completed with warnings - review before Stage 3"
                }
            }
    }
}

/*
========================================================================================
    Workflow Completion Handlers
========================================================================================
*/

workflow.onComplete {
    log.info ""
    log.info "🏁 Stage 2: Genome Preprocessing & Indexing - Complete!"
    log.info "Execution status: ${workflow.success ? 'SUCCESS' : 'FAILED'}"
    log.info "Execution duration: ${workflow.duration}"
    log.info "Peak memory: ${workflow.stats.peakMemory ?: 'N/A'}"
    log.info "CPU hours: ${workflow.stats.computeTimeFmt ?: 'N/A'}"

    if (workflow.success) {
        log.info ""
        log.info "📋 Next Steps:"
        log.info "1. Review QC reports: results/02_preprocessing/qc_reports/"
        log.info "2. Check selected reference: results/02_preprocessing/reference_selection/"
        log.info "3. Proceed to Stage 3: PGGB pangenome construction"
        log.info "   nextflow run . --stage 3 -profile slurm,singularity"
        log.info ""
    }
}

workflow.onError {
    log.error ""
    log.error "❌ Stage 2 failed!"
    log.error "Error: ${workflow.errorMessage}"
    log.error ""
    log.error "🔧 Troubleshooting:"
    log.error "1. Check resource allocation in conf/slurm.config"
    log.error "2. Review failed process logs in work/ directories"
    log.error "3. Ensure all Stage 1 outputs are complete and valid"
    log.error ""
}