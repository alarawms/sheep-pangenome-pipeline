/*
========================================================================================
    Stage 2: Genome Preprocessing & Indexing Subworkflow
========================================================================================
    Description: Orchestrates genome standardization, QC, indexing, and reference selection
    Input: Validated genomes from Stage 1
    Output: Standardized, indexed genomes with selected reference
========================================================================================
*/

include { STANDARDIZE_GENOME      } from '../../modules/local/standardize_genome'
include { GENOME_QC_EXTENDED      } from '../../modules/local/genome_qc_extended'
include { BWA_INDEX               } from '../../modules/local/bwa_index'
include { MINIMAP2_INDEX          } from '../../modules/local/minimap2_index'
include { SAMTOOLS_FAIDX          } from '../../modules/local/samtools_faidx'
include { SELECT_REFERENCE_MANUAL  } from '../../modules/local/select_reference_manual'

workflow PREPROCESSING {
    take:
    genome_input_ch    // Channel: [meta, genome.fa, metadata.json]

    main:
    // Initialize version collection
    ch_versions = Channel.empty()

    //
    // MODULE: Standardize genome chromosome naming and format
    //
    STANDARDIZE_GENOME(
        genome_input_ch
    )
    ch_versions = ch_versions.mix(STANDARDIZE_GENOME.out.versions)

    // Combine standardized genomes with original metadata for QC
    standardized_for_qc = STANDARDIZE_GENOME.out.fasta
        .join(genome_input_ch.map { meta, fasta, metadata -> [meta, metadata] })
        .map { meta, std_fasta, metadata -> [meta, std_fasta, metadata] }

    //
    // MODULE: Extended quality control and tier assessment
    //
    GENOME_QC_EXTENDED(
        standardized_for_qc
    )
    ch_versions = ch_versions.mix(GENOME_QC_EXTENDED.out.versions)

    // Create parallel indexing channels
    genomes_for_indexing = STANDARDIZE_GENOME.out.fasta

    //
    // MODULE: Create BWA-MEM2 index for short-read validation (multi-threaded)
    //
    BWA_INDEX(
        genomes_for_indexing
    )
    ch_versions = ch_versions.mix(BWA_INDEX.out.versions)

    //
    // MODULE: Create minimap2 index for long-read alignment
    //
    MINIMAP2_INDEX(
        genomes_for_indexing
    )
    ch_versions = ch_versions.mix(MINIMAP2_INDEX.out.versions)

    //
    // MODULE: Create samtools faidx index for sequence access
    //
    SAMTOOLS_FAIDX(
        genomes_for_indexing
    )
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)

    //
    // MODULE: Select reference genome manually (optional)
    //
    if (params.enable_reference_selection) {
        // Prepare data for manual reference selection
        // Collect all standardized genomes and QC stats for manual selection
        reference_selection_data = STANDARDIZE_GENOME.out.fasta
            .join(GENOME_QC_EXTENDED.out.stats, by: 0)
            .collect()
            .map { collected_data ->
                def genome_files = collected_data.collect { it[1] }
                def stats_files = collected_data.collect { it[2] }
                [collected_data[0][0], genome_files, stats_files]  // Use first meta as representative
            }

        SELECT_REFERENCE_MANUAL(
            reference_selection_data
        )
        ch_versions = ch_versions.mix(SELECT_REFERENCE_MANUAL.out.versions)

        // Set reference outputs when enabled
        reference_genome_out = SELECT_REFERENCE_MANUAL.out.reference_genome
        reference_metadata_out = SELECT_REFERENCE_MANUAL.out.reference_metadata
        reference_report_out = SELECT_REFERENCE_MANUAL.out.selection_report
    } else {
        // Create empty channels when reference selection is disabled
        reference_genome_out = Channel.empty()
        reference_metadata_out = Channel.empty()
        reference_report_out = Channel.empty()
    }

    // Create comprehensive output channels
    preprocessed_genomes = STANDARDIZE_GENOME.out.fasta
        .join(STANDARDIZE_GENOME.out.mapping, by: 0)
        .join(GENOME_QC_EXTENDED.out.report, by: 0)
        .join(GENOME_QC_EXTENDED.out.stats, by: 0)
        .join(BWA_INDEX.out.index, by: 0)
        .join(MINIMAP2_INDEX.out.index, by: 0)
        .join(SAMTOOLS_FAIDX.out.index, by: 0)

    emit:
    // Standardized genomes and mapping
    standardized_genomes = STANDARDIZE_GENOME.out.fasta              // [meta, standardized.fa]
    chromosome_mapping   = STANDARDIZE_GENOME.out.mapping            // [meta, mapping.tsv]
    standardization_logs = STANDARDIZE_GENOME.out.log               // [meta, log.txt]

    // Quality control results
    qc_html_reports      = GENOME_QC_EXTENDED.out.report           // [meta, report.html]
    qc_stats_json        = GENOME_QC_EXTENDED.out.stats            // [meta, stats.json]
    qc_logs              = GENOME_QC_EXTENDED.out.log              // [meta, log.txt]

    // Indexing results
    bwa_indices          = BWA_INDEX.out.index                      // [meta, bwa-mem2_files]
    bwa_logs             = BWA_INDEX.out.log                        // [meta, log.txt]
    minimap2_indices     = MINIMAP2_INDEX.out.index                 // [meta, .mmi]
    minimap2_logs        = MINIMAP2_INDEX.out.log                   // [meta, log.txt]
    samtools_indices     = SAMTOOLS_FAIDX.out.index                 // [meta, .fai]
    samtools_logs        = SAMTOOLS_FAIDX.out.log                   // [meta, log.txt]

    // Reference selection results (optional)
    reference_genome     = reference_genome_out     // [selected_meta, reference.fa] or empty
    reference_metadata   = reference_metadata_out   // reference_metadata.json or empty
    selection_report     = reference_report_out     // selection_report.json or empty

    // Complete preprocessed data
    preprocessed_genomes = preprocessed_genomes                     // All data combined

    // Version information
    versions             = ch_versions                              // versions.yml
}

// Completion handlers removed to prevent false completion messages