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
include { SELECT_REFERENCE        } from '../../modules/local/select_reference'

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
    // MODULE: Create BWA index for short-read validation
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

    // Combine all data for reference selection
    // Each genome needs: [meta, genome.fa, stats.json, metadata.json]
    combined_genome_data = STANDARDIZE_GENOME.out.fasta
        .join(GENOME_QC_EXTENDED.out.stats, by: 0)
        .join(genome_input_ch.map { meta, fasta, metadata -> [meta, metadata] }, by: 0)
        .map { meta, std_fasta, stats_json, metadata ->
            [meta, std_fasta, stats_json, metadata]
        }
        .collect()

    //
    // MODULE: Select optimal reference genome
    //
    SELECT_REFERENCE(
        combined_genome_data
    )
    ch_versions = ch_versions.mix(SELECT_REFERENCE.out.versions)

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
    bwa_indices          = BWA_INDEX.out.index                      // [meta, bwa_files]
    bwa_logs             = BWA_INDEX.out.log                        // [meta, log.txt]
    minimap2_indices     = MINIMAP2_INDEX.out.index                 // [meta, .mmi]
    minimap2_logs        = MINIMAP2_INDEX.out.log                   // [meta, log.txt]
    samtools_indices     = SAMTOOLS_FAIDX.out.index                 // [meta, .fai]
    samtools_logs        = SAMTOOLS_FAIDX.out.log                   // [meta, log.txt]

    // Reference selection results
    reference_genome     = SELECT_REFERENCE.out.reference_genome    // [selected_meta, reference.fa]
    reference_metadata   = SELECT_REFERENCE.out.reference_metadata  // reference_metadata.json
    selection_report     = SELECT_REFERENCE.out.selection_report    // selection_report.json

    // Complete preprocessed data
    preprocessed_genomes = preprocessed_genomes                     // All data combined

    // Version information
    versions             = ch_versions                              // versions.yml
}

workflow.onError {
    log.error "Stage 2 Preprocessing workflow failed: ${workflow.errorMessage}"
}

workflow.onComplete {
    log.info "Stage 2 Preprocessing completed successfully"
    log.info "Reference genome selected and all genomes standardized and indexed"
}