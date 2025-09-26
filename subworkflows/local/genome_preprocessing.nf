/*
========================================================================================
    Stage 2: Genome Preprocessing & Indexing Subworkflow
========================================================================================
    Description: Main orchestrator for genome standardization, QC, and indexing
    Input: Validated genomes from Stage 1
    Output: Analysis-ready genomes with indices for Stage 3
========================================================================================
*/

// Import subworkflows
include { STANDARDIZATION } from './standardization'
include { QUALITY_CONTROL } from './quality_control'
include { INDEX_GENERATION } from './index_generation'
include { REFERENCE_SELECTION } from './reference_selection'

workflow GENOME_PREPROCESSING {

    take:
    genomes              // channel: [meta, genome.fa] from Stage 1
    stage1_metadata      // channel: [meta, metadata.json] from Stage 1
    stage1_validation    // channel: [meta, validation.json] from Stage 1

    main:

    log.info "🔄 Starting Stage 2: Genome Preprocessing & Indexing"

    // ================================================================================
    // Input Validation & Contract Enforcement
    // ================================================================================

    // Validate Stage 1 outputs meet Stage 2 requirements
    validated_inputs = stage1_validation
        .join(genomes, by: 0)
        .join(stage1_metadata, by: 0)
        .map { meta, validation_file, genome_file, metadata_file ->

            // Parse validation results
            def validation = new groovy.json.JsonSlurper().parse(validation_file)
            def metadata = new groovy.json.JsonSlurper().parse(metadata_file)

            // Contract validation
            if (validation.validation_status != "PASS") {
                log.error "❌ Sample ${meta.id} failed Stage 1 validation - cannot proceed to Stage 2"
                return null
            }

            // Genome size validation (sheep-specific)
            if (validation.total_length < 2.4e9 || validation.total_length > 3.2e9) {
                log.warn "⚠️  Sample ${meta.id} genome size outside expected range (2.4-3.2Gb)"
            }

            // Quality score assessment
            def quality_tier = "Unknown"
            if (metadata.containsKey('quality_score')) {
                quality_tier = metadata.quality_score
            } else if (validation.total_sequences <= 1000 && validation.n_content < 1.0) {
                quality_tier = "A"
            } else if (validation.total_sequences <= 5000 && validation.n_content < 2.0) {
                quality_tier = "B"
            } else {
                quality_tier = "C"
            }

            // Enhanced metadata for Stage 2
            def enhanced_meta = meta + [
                genome_length: validation.total_length,
                sequence_count: validation.total_sequences,
                gc_content: validation.gc_content,
                n_content: validation.n_content,
                quality_tier: quality_tier,
                stage1_validated: true
            ]

            return [enhanced_meta, genome_file, metadata_file, validation_file]
        }
        .filter { it != null }

    // Report input validation results
    validated_inputs.count().subscribe { count ->
        log.info "✅ ${count} genomes passed Stage 1→2 validation"
    }

    // ================================================================================
    // Phase 1: Genome Standardization
    // ================================================================================

    log.info "🏷️  Phase 1: Genome Standardization"

    STANDARDIZATION(
        validated_inputs.map { enhanced_meta, genome, metadata, validation ->
            [enhanced_meta, genome]
        }
    )

    // ================================================================================
    // Phase 2: Extended Quality Control
    // ================================================================================

    log.info "🔍 Phase 2: Extended Quality Control"

    QUALITY_CONTROL(
        STANDARDIZATION.out.standardized_genomes,
        validated_inputs.map { enhanced_meta, genome, metadata, validation ->
            [enhanced_meta, metadata]
        }
    )

    // ================================================================================
    // Phase 3: Index Generation
    // ================================================================================

    log.info "🗂️  Phase 3: Multi-Platform Index Generation"

    INDEX_GENERATION(
        STANDARDIZATION.out.standardized_genomes
    )

    // ================================================================================
    // Phase 4: Reference Genome Selection
    // ================================================================================

    log.info "🎯 Phase 4: Reference Genome Selection"

    REFERENCE_SELECTION(
        STANDARDIZATION.out.standardized_genomes.collect(),
        QUALITY_CONTROL.out.extended_stats.collect(),
        validated_inputs.map { enhanced_meta, genome, metadata, validation ->
            [enhanced_meta, metadata]
        }.collect()
    )

    // ================================================================================
    // Final Validation & Reporting
    // ================================================================================

    // Collect all processed genomes with comprehensive metadata
    processed_genomes = STANDARDIZATION.out.standardized_genomes
        .join(QUALITY_CONTROL.out.extended_stats, by: 0)
        .join(INDEX_GENERATION.out.bwa_indices, by: 0)
        .join(INDEX_GENERATION.out.minimap2_indices, by: 0)
        .join(INDEX_GENERATION.out.samtools_indices, by: 0)
        .map { meta, genome, stats, bwa_idx, mm2_idx, sam_idx ->

            def processing_metadata = [
                sample_id: meta.id,
                stage2_processed: true,
                standardized: true,
                indexed: [bwa: true, minimap2: true, samtools: true],
                quality_tier: meta.quality_tier,
                genome_stats: stats,
                processing_timestamp: new Date().format("yyyy-MM-dd'T'HH:mm:ss.SSS")
            ]

            return [meta, genome, processing_metadata]
        }

    // Generate summary statistics
    stage2_summary = processed_genomes
        .collect()
        .map { genome_list ->
            def total = genome_list.size()
            def by_quality = genome_list.groupBy { it[0].quality_tier }
            def avg_size = genome_list.collect { it[0].genome_length }.sum() / total

            def summary = [
                total_genomes: total,
                avg_genome_size: avg_size,
                quality_distribution: by_quality.collectEntries { tier, genomes ->
                    [tier, genomes.size()]
                },
                standardized: total,
                indexed: total,
                ready_for_stage3: total
            ]

            log.info ""
            log.info "📊 Stage 2 Processing Summary:"
            log.info "  Total genomes processed: ${total}"
            log.info "  Average genome size: ${(avg_size/1e9).round(2)} Gb"
            log.info "  Quality distribution: ${summary.quality_distribution}"
            log.info "  All genomes standardized: ✅"
            log.info "  All indices generated: ✅"
            log.info "  Ready for Stage 3: ${total} genomes"
            log.info ""

            return summary
        }

    emit:
    // Primary outputs for Stage 3
    standardized_genomes = STANDARDIZATION.out.standardized_genomes
    reference_genome = REFERENCE_SELECTION.out.reference_genome
    reference_metadata = REFERENCE_SELECTION.out.reference_metadata

    // Index outputs
    bwa_indices = INDEX_GENERATION.out.bwa_indices
    minimap2_indices = INDEX_GENERATION.out.minimap2_indices
    samtools_indices = INDEX_GENERATION.out.samtools_indices

    // Quality control outputs
    qc_reports = QUALITY_CONTROL.out.qc_reports
    extended_stats = QUALITY_CONTROL.out.extended_stats
    busco_results = QUALITY_CONTROL.out.busco_results

    // Standardization outputs
    chromosome_maps = STANDARDIZATION.out.chromosome_maps
    header_maps = STANDARDIZATION.out.header_maps

    // Metadata and versioning
    processing_metadata = processed_genomes.map { meta, genome, proc_meta -> [meta, proc_meta] }
    versions = STANDARDIZATION.out.versions
        .mix(QUALITY_CONTROL.out.versions)
        .mix(INDEX_GENERATION.out.versions)
        .mix(REFERENCE_SELECTION.out.versions)

    // Summary for pipeline orchestration
    summary = stage2_summary
}