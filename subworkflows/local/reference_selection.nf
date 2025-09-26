/*
========================================================================================
    Reference Genome Selection Subworkflow
========================================================================================
    Description: Identify optimal reference genome for PGGB pangenome construction
    Algorithm: Multi-criteria scoring based on quality, completeness, and suitability
========================================================================================
*/

include { SELECT_REFERENCE } from '../../modules/local/select_reference'

workflow REFERENCE_SELECTION {

    take:
    standardized_genomes  // channel: [meta, genome.fa] - collected
    extended_stats        // channel: [meta, stats.json] - collected
    metadata_collection   // channel: [meta, metadata.json] - collected

    main:

    log.info "🎯 Starting reference genome selection process"

    // ================================================================================
    // Reference Selection Algorithm
    // ================================================================================

    log.info "🏆 Applying multi-criteria reference selection algorithm"

    // Combine all data for comprehensive evaluation
    combined_data = Channel.from(standardized_genomes)
        .flatten()
        .map { it -> [it[0].id, it] }
        .combine(
            Channel.from(extended_stats)
                .flatten()
                .map { it -> [it[0].id, it] },
            by: 0
        )
        .combine(
            Channel.from(metadata_collection)
                .flatten()
                .map { it -> [it[0].id, it] },
            by: 0
        )
        .map { id, genome_data, stats_data, meta_data ->
            [genome_data[0], genome_data[1], stats_data[1], meta_data[1]]
        }

    // Run reference selection analysis
    SELECT_REFERENCE(
        combined_data.collect()
    )

    // ================================================================================
    // Reference Validation & Reporting
    // ================================================================================

    // Process selection results
    reference_analysis = SELECT_REFERENCE.out.selection_report
        .map { report_file ->
            def selection = new groovy.json.JsonSlurper().parse(report_file)

            log.info ""
            log.info "🎯 Reference Genome Selection Results:"
            log.info "  Selected reference: ${selection.selected_reference.sample_id}"
            log.info "  Quality tier: ${selection.selected_reference.quality_tier}"
            log.info "  Selection score: ${selection.selected_reference.selection_score}/100"
            log.info "  Assembly level: ${selection.selected_reference.assembly_level}"
            log.info "  BUSCO completeness: ${selection.selected_reference.busco_completeness}%"
            log.info "  Scaffold N50: ${(selection.selected_reference.scaffold_n50/1e6).round(1)} Mb"
            log.info ""
            log.info "📊 Selection criteria breakdown:"
            selection.selected_reference.score_breakdown.each { criterion, score ->
                log.info "  ${criterion}: ${score}/100"
            }
            log.info ""
            log.info "🔍 Alternative candidates:"
            selection.alternative_candidates.each { candidate ->
                log.info "  ${candidate.sample_id} (Tier ${candidate.quality_tier}): ${candidate.selection_score}/100"
            }
            log.info ""

            // Validation checks for selected reference
            def validations = []
            if (selection.selected_reference.quality_tier in ['A+', 'A']) {
                validations.add("✅ High quality assembly (Tier ${selection.selected_reference.quality_tier})")
            } else {
                validations.add("⚠️  Lower quality reference (Tier ${selection.selected_reference.quality_tier})")
            }

            if (selection.selected_reference.busco_completeness >= 95) {
                validations.add("✅ Highly complete gene space (${selection.selected_reference.busco_completeness}%)")
            } else if (selection.selected_reference.busco_completeness >= 90) {
                validations.add("⚠️  Good but not optimal completeness (${selection.selected_reference.busco_completeness}%)")
            } else {
                validations.add("❌ Low gene space completeness (${selection.selected_reference.busco_completeness}%)")
            }

            if (selection.selected_reference.scaffold_n50 >= 50e6) {
                validations.add("✅ Excellent contiguity (N50 = ${(selection.selected_reference.scaffold_n50/1e6).round(1)} Mb)")
            } else if (selection.selected_reference.scaffold_n50 >= 10e6) {
                validations.add("⚠️  Good contiguity (N50 = ${(selection.selected_reference.scaffold_n50/1e6).round(1)} Mb)")
            } else {
                validations.add("❌ Poor contiguity (N50 = ${(selection.selected_reference.scaffold_n50/1e6).round(1)} Mb)")
            }

            log.info "🔬 Reference validation:"
            validations.each { log.info "  ${it}" }
            log.info ""

            // Determine suitability for pangenome construction
            def pangenome_suitability = "Unknown"
            if (selection.selected_reference.selection_score >= 85) {
                pangenome_suitability = "Excellent"
                log.info "🌟 Selected reference is EXCELLENT for pangenome construction"
            } else if (selection.selected_reference.selection_score >= 70) {
                pangenome_suitability = "Good"
                log.info "✅ Selected reference is GOOD for pangenome construction"
            } else if (selection.selected_reference.selection_score >= 55) {
                pangenome_suitability = "Adequate"
                log.info "⚠️  Selected reference is ADEQUATE for pangenome construction"
            } else {
                pangenome_suitability = "Suboptimal"
                log.info "❌ Selected reference is SUBOPTIMAL for pangenome construction"
                log.info "   Consider improving reference quality or selecting alternative"
            }

            return [
                selected_reference_id: selection.selected_reference.sample_id,
                selection_score: selection.selected_reference.selection_score,
                pangenome_suitability: pangenome_suitability,
                validation_results: validations,
                total_candidates: selection.total_candidates,
                alternative_count: selection.alternative_candidates.size()
            ]
        }

    emit:
    reference_genome = SELECT_REFERENCE.out.reference_genome
    reference_metadata = SELECT_REFERENCE.out.reference_metadata
    selection_report = SELECT_REFERENCE.out.selection_report
    reference_analysis = reference_analysis
    versions = SELECT_REFERENCE.out.versions
}