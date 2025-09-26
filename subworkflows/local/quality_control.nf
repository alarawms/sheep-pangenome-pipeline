/*
========================================================================================
    Extended Quality Control Subworkflow
========================================================================================
    Description: Comprehensive QC beyond Stage 1 basic validation
    Features: BUSCO completeness, detailed statistics, contamination detection
========================================================================================
*/

include { GENOME_STATS_EXTENDED } from '../../modules/local/genome_stats_extended'
include { BUSCO_ASSESSMENT } from '../../modules/local/busco_assessment'
include { GENERATE_QC_REPORT } from '../../modules/local/generate_qc_report'

workflow QUALITY_CONTROL {

    take:
    genomes      // channel: [meta, standardized_genome.fa]
    metadata     // channel: [meta, stage1_metadata.json]

    main:

    log.info "🔍 Starting extended quality control assessment"

    // ================================================================================
    // Step 1: Extended Genome Statistics
    // ================================================================================

    log.info "📊 Step 1: Computing extended genome statistics"

    GENOME_STATS_EXTENDED(genomes)

    // ================================================================================
    // Step 2: BUSCO Completeness Assessment
    // ================================================================================

    log.info "🧬 Step 2: BUSCO completeness assessment (Mammalia lineage)"

    BUSCO_ASSESSMENT(genomes)

    // ================================================================================
    // Step 3: Quality Tier Assignment
    // ================================================================================

    // Combine all QC metrics for comprehensive assessment
    combined_qc = genomes
        .join(GENOME_STATS_EXTENDED.out.stats, by: 0)
        .join(BUSCO_ASSESSMENT.out.summary, by: 0)
        .join(metadata, by: 0)
        .map { meta, genome, extended_stats, busco_summary, stage1_meta ->

            // Parse QC results
            def stats = new groovy.json.JsonSlurper().parse(extended_stats)
            def busco = new groovy.json.JsonSlurper().parse(busco_summary)
            def s1_meta = new groovy.json.JsonSlurper().parse(stage1_meta)

            // Quality scoring algorithm
            def quality_score = 0

            // Assembly contiguity (30%)
            if (stats.scaffold_n50 > 50e6) quality_score += 30      // Excellent contiguity
            else if (stats.scaffold_n50 > 10e6) quality_score += 20 // Good contiguity
            else if (stats.scaffold_n50 > 1e6) quality_score += 10  // Fair contiguity

            // BUSCO completeness (40%)
            def busco_complete = busco.complete_buscos_percentage ?: 0
            if (busco_complete >= 95) quality_score += 40      // Excellent completeness
            else if (busco_complete >= 90) quality_score += 30 // Good completeness
            else if (busco_complete >= 80) quality_score += 20 // Fair completeness
            else if (busco_complete >= 70) quality_score += 10 // Poor completeness

            // Sequence quality (20%)
            def n_content = stats.n_percentage ?: s1_meta.n_content ?: 5.0
            if (n_content <= 0.5) quality_score += 20      // Excellent sequence quality
            else if (n_content <= 1.0) quality_score += 15 // Good sequence quality
            else if (n_content <= 2.0) quality_score += 10 // Fair sequence quality
            else if (n_content <= 3.0) quality_score += 5  // Poor sequence quality

            // Assembly size accuracy (10%)
            def expected_size = 2.8e9  // Expected sheep genome size
            def size_diff = Math.abs(stats.total_length - expected_size) / expected_size
            if (size_diff <= 0.05) quality_score += 10     // Within 5%
            else if (size_diff <= 0.10) quality_score += 7 // Within 10%
            else if (size_diff <= 0.20) quality_score += 5 // Within 20%

            // Assign quality tier based on total score
            def quality_tier = "Unknown"
            def tier_rationale = []

            if (quality_score >= 85) {
                quality_tier = "A+"
                tier_rationale.add("Excellent assembly with >95% BUSCO completeness")
            } else if (quality_score >= 75) {
                quality_tier = "A"
                tier_rationale.add("High quality assembly suitable for reference use")
            } else if (quality_score >= 65) {
                quality_tier = "B+"
                tier_rationale.add("Good quality assembly with minor limitations")
            } else if (quality_score >= 50) {
                quality_tier = "B"
                tier_rationale.add("Adequate quality for pangenome inclusion")
            } else if (quality_score >= 35) {
                quality_tier = "C+"
                tier_rationale.add("Acceptable quality with noted limitations")
            } else if (quality_score >= 20) {
                quality_tier = "C"
                tier_rationale.add("Lower quality but potentially useful")
            } else {
                quality_tier = "D"
                tier_rationale.add("Poor quality - recommend exclusion")
            }

            // Add specific quality metrics to rationale
            if (stats.scaffold_n50 > 50e6) tier_rationale.add("Excellent contiguity (N50 > 50Mb)")
            if (busco_complete >= 95) tier_rationale.add("Highly complete gene space")
            if (n_content <= 1.0) tier_rationale.add("Low gap content")

            def enhanced_meta = meta + [
                quality_score: quality_score,
                quality_tier: quality_tier,
                quality_rationale: tier_rationale,
                busco_complete: busco_complete,
                scaffold_n50: stats.scaffold_n50,
                n_content: n_content
            ]

            log.info "🏆 ${meta.id}: Tier ${quality_tier} (score: ${quality_score}/100)"

            return [enhanced_meta, genome, extended_stats, busco_summary]
        }

    // ================================================================================
    // Step 4: Comprehensive QC Report Generation
    // ================================================================================

    log.info "📋 Step 4: Generating comprehensive QC reports"

    GENERATE_QC_REPORT(
        combined_qc.map { meta, genome, stats, busco -> [meta, stats] },
        combined_qc.map { meta, genome, stats, busco -> [meta, busco] }
    )

    // ================================================================================
    // Quality Control Summary
    // ================================================================================

    qc_summary = combined_qc
        .collect()
        .map { qc_list ->
            def total = qc_list.size()
            def by_tier = qc_list.groupBy { it[0].quality_tier }
            def avg_score = qc_list.collect { it[0].quality_score }.sum() / total
            def avg_busco = qc_list.collect { it[0].busco_complete }.sum() / total

            log.info ""
            log.info "🔍 Quality Control Summary:"
            log.info "  Total genomes assessed: ${total}"
            log.info "  Average quality score: ${avg_score.round(1)}/100"
            log.info "  Average BUSCO completeness: ${avg_busco.round(1)}%"
            log.info "  Quality distribution:"
            by_tier.each { tier, genomes ->
                log.info "    ${tier}: ${genomes.size()} genomes"
            }

            // Quality recommendations
            def tier_a_plus = by_tier['A+']?.size() ?: 0
            def tier_a = by_tier['A']?.size() ?: 0
            def tier_b_plus = by_tier['B+']?.size() ?: 0
            def tier_b = by_tier['B']?.size() ?: 0

            def recommended_for_reference = tier_a_plus + tier_a
            def suitable_for_pangenome = tier_a_plus + tier_a + tier_b_plus + tier_b

            log.info ""
            log.info "📈 Quality Recommendations:"
            log.info "  Recommended for reference: ${recommended_for_reference} genomes (Tier A+ or A)"
            log.info "  Suitable for pangenome: ${suitable_for_pangenome} genomes (Tier B or better)"

            return [
                total: total,
                average_quality_score: avg_score,
                average_busco_completeness: avg_busco,
                quality_distribution: by_tier.collectEntries { tier, genomes ->
                    [tier, genomes.size()]
                },
                recommended_for_reference: recommended_for_reference,
                suitable_for_pangenome: suitable_for_pangenome
            ]
        }

    emit:
    extended_stats = combined_qc.map { meta, genome, stats, busco -> [meta, stats] }
    busco_results = combined_qc.map { meta, genome, stats, busco -> [meta, busco] }
    qc_reports = GENERATE_QC_REPORT.out.qc_report
    quality_enhanced_genomes = combined_qc.map { meta, genome, stats, busco -> [meta, genome] }
    qc_summary = qc_summary
    versions = GENOME_STATS_EXTENDED.out.versions
        .mix(BUSCO_ASSESSMENT.out.versions)
        .mix(GENERATE_QC_REPORT.out.versions)
}