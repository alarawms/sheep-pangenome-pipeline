/*
========================================================================================
    Index Generation Subworkflow
========================================================================================
    Description: Generate BWA, minimap2, and samtools indices for standardized genomes
    Purpose: Enable efficient alignment and pangenome validation in Stage 3
========================================================================================
*/

include { CREATE_BWA_INDEX } from '../../modules/local/create_bwa_index'
include { CREATE_MINIMAP2_INDEX } from '../../modules/local/create_minimap2_index'
include { CREATE_SAMTOOLS_INDEX } from '../../modules/local/create_samtools_index'

workflow INDEX_GENERATION {

    take:
    standardized_genomes  // channel: [meta, standardized_genome.fa]

    main:

    log.info "🗂️  Starting multi-platform index generation"

    // ================================================================================
    // Step 1: BWA Index Generation
    // ================================================================================

    log.info "🧬 Step 1: Creating BWA-MEM2 indices for alignment validation"

    CREATE_BWA_INDEX(standardized_genomes)

    // ================================================================================
    // Step 2: Minimap2 Index Generation
    // ================================================================================

    log.info "🗺️  Step 2: Creating minimap2 indices for long-read alignment"

    CREATE_MINIMAP2_INDEX(standardized_genomes)

    // ================================================================================
    // Step 3: Samtools Index Generation
    // ================================================================================

    log.info "🔧 Step 3: Creating samtools FAIDX indices for random access"

    CREATE_SAMTOOLS_INDEX(standardized_genomes)

    // ================================================================================
    // Index Validation & Summary
    // ================================================================================

    // Validate all indices were created successfully
    all_indices = CREATE_BWA_INDEX.out.index
        .join(CREATE_MINIMAP2_INDEX.out.index, by: 0)
        .join(CREATE_SAMTOOLS_INDEX.out.index, by: 0)
        .map { meta, bwa_idx, mm2_idx, sam_idx ->

            // Verify index files exist
            def bwa_valid = bwa_idx && bwa_idx.exists()
            def mm2_valid = mm2_idx && mm2_idx.exists()
            def sam_valid = sam_idx && sam_idx.exists()

            def index_status = [
                sample_id: meta.id,
                bwa_index: bwa_valid,
                minimap2_index: mm2_valid,
                samtools_index: sam_valid,
                all_indices_created: bwa_valid && mm2_valid && sam_valid,
                index_timestamp: new Date().format("yyyy-MM-dd'T'HH:mm:ss.SSS")
            ]

            if (index_status.all_indices_created) {
                log.info "✅ All indices created for ${meta.id}"
            } else {
                log.warn "⚠️  Incomplete indexing for ${meta.id}: BWA:${bwa_valid}, MM2:${mm2_valid}, SAM:${sam_valid}"
            }

            return [meta, index_status]
        }

    // Generate indexing summary
    indexing_summary = all_indices
        .collect()
        .map { index_list ->
            def total = index_list.size()
            def successful = index_list.count { it[1].all_indices_created }
            def bwa_success = index_list.count { it[1].bwa_index }
            def mm2_success = index_list.count { it[1].minimap2_index }
            def sam_success = index_list.count { it[1].samtools_index }

            log.info ""
            log.info "🗂️  Index Generation Summary:"
            log.info "  Total genomes: ${total}"
            log.info "  Complete indexing: ${successful}/${total} (${(successful/total*100).round(1)}%)"
            log.info "  BWA indices: ${bwa_success}/${total}"
            log.info "  Minimap2 indices: ${mm2_success}/${total}"
            log.info "  Samtools indices: ${sam_success}/${total}"

            if (successful == total) {
                log.info "🎉 All genomes successfully indexed for downstream analysis"
            } else {
                log.warn "⚠️  ${total - successful} genomes have incomplete indexing"
            }

            return [
                total: total,
                successful: successful,
                success_rate: successful / total,
                bwa_success: bwa_success,
                minimap2_success: mm2_success,
                samtools_success: sam_success
            ]
        }

    emit:
    bwa_indices = CREATE_BWA_INDEX.out.index
    minimap2_indices = CREATE_MINIMAP2_INDEX.out.index
    samtools_indices = CREATE_SAMTOOLS_INDEX.out.index
    index_status = all_indices
    indexing_summary = indexing_summary
    versions = CREATE_BWA_INDEX.out.versions
        .mix(CREATE_MINIMAP2_INDEX.out.versions)
        .mix(CREATE_SAMTOOLS_INDEX.out.versions)
}