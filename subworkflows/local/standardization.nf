/*
========================================================================================
    Genome Standardization Subworkflow
========================================================================================
    Description: Normalize chromosome naming, sequence headers, and coordinate systems
    Focus: Sheep-specific chromosome naming (1-27, MT) and header standardization
========================================================================================
*/

include { NORMALIZE_HEADERS } from '../../modules/local/normalize_headers'
include { RENAME_CHROMOSOMES } from '../../modules/local/rename_chromosomes'

workflow STANDARDIZATION {

    take:
    genomes  // channel: [meta, genome.fa]

    main:

    log.info "🏷️  Starting genome standardization process"

    // ================================================================================
    // Step 1: Chromosome Naming Standardization
    // ================================================================================

    log.info "🧬 Step 1: Sheep chromosome naming standardization (1-27, MT)"

    RENAME_CHROMOSOMES(genomes)

    // ================================================================================
    // Step 2: FASTA Header Normalization
    // ================================================================================

    log.info "📝 Step 2: FASTA header normalization"

    NORMALIZE_HEADERS(RENAME_CHROMOSOMES.out.renamed_genomes)

    // ================================================================================
    // Step 3: Quality Assessment of Standardization
    // ================================================================================

    // Generate standardization report
    standardization_report = NORMALIZE_HEADERS.out.normalized_genomes
        .join(RENAME_CHROMOSOMES.out.chromosome_map, by: 0)
        .join(NORMALIZE_HEADERS.out.header_map, by: 0)
        .map { meta, genome, chr_map, hdr_map ->

            // Parse mapping files
            def chr_mappings = chr_map.readLines().size() - 1  // subtract header
            def hdr_mappings = hdr_map.readLines().size() - 1  // subtract header

            def report = [
                sample_id: meta.id,
                chromosomes_renamed: chr_mappings,
                headers_normalized: hdr_mappings,
                standardization_complete: true,
                sheep_chromosome_count: chr_map.readLines().findAll {
                    it.contains('chr1') || it.contains('chr2') || it.contains('chrMT')
                }.size() > 0
            ]

            log.info "✅ Standardized ${meta.id}: ${chr_mappings} chromosomes, ${hdr_mappings} headers"

            return [meta, report]
        }

    emit:
    standardized_genomes = NORMALIZE_HEADERS.out.normalized_genomes
    chromosome_maps = RENAME_CHROMOSOMES.out.chromosome_map
    header_maps = NORMALIZE_HEADERS.out.header_map
    standardization_reports = standardization_report
    versions = RENAME_CHROMOSOMES.out.versions.mix(NORMALIZE_HEADERS.out.versions)
}