/*
========================================================================================
    Stage 3: Pangenome Graph Construction Subworkflow
========================================================================================
    Description: Orchestrates PGGB-based pangenome graph construction and analysis
    Input: Standardized genomes from Stage 2 with selected reference
    Output: Complete pangenome graph with validation, statistics, and visualizations
========================================================================================
*/

include { PGGB_CONSTRUCT    } from '../../modules/local/pggb_construct'
include { GRAPH_VALIDATE    } from '../../modules/local/graph_validate'
include { GRAPH_STATS       } from '../../modules/local/graph_stats'

workflow GRAPH_CONSTRUCTION {
    take:
    standardized_genomes_ch    // Channel: [meta, standardized.fa]
    reference_metadata_ch      // Channel: reference_metadata.json

    main:
    // Initialize version collection
    ch_versions = Channel.empty()

    //
    // Prepare genomes for PGGB construction
    //
    // Collect all standardized genomes into single input for PGGB
    genomes_for_pggb = standardized_genomes_ch
        .map { meta, fasta -> fasta }
        .collect()
        .map { fasta_files ->
            // Create meta for pangenome
            def pangenome_meta = [
                id: "sheep_pangenome",
                genome_count: fasta_files.size(),
                stage: "graph_construction"
            ]
            [pangenome_meta, fasta_files]
        }

    //
    // MODULE: PGGB pangenome graph construction
    //
    PGGB_CONSTRUCT(
        genomes_for_pggb,
        reference_metadata_ch
    )
    ch_versions = ch_versions.mix(PGGB_CONSTRUCT.out.versions)

    // Combine graph outputs for validation and statistics
    graph_data_for_validation = PGGB_CONSTRUCT.out.graph
        .join(PGGB_CONSTRUCT.out.odgi_graph, by: 0)

    //
    // MODULE: Graph validation and quality control
    //
    GRAPH_VALIDATE(
        graph_data_for_validation
    )
    ch_versions = ch_versions.mix(GRAPH_VALIDATE.out.versions)

    //
    // MODULE: Comprehensive graph statistics
    //
    GRAPH_STATS(
        graph_data_for_validation
    )
    ch_versions = ch_versions.mix(GRAPH_STATS.out.versions)

    emit:
    // Primary graph outputs
    pangenome_graph      = PGGB_CONSTRUCT.out.graph              // [meta, sheep_pangenome.gfa]
    odgi_graph           = PGGB_CONSTRUCT.out.odgi_graph         // [meta, sheep_pangenome.og]
    smoothed_graph       = PGGB_CONSTRUCT.out.smoothed_graph     // [meta, sheep_pangenome.smooth.gfa]
    graph_visualization  = PGGB_CONSTRUCT.out.visualization      // [meta, sheep_pangenome.viz.png]
    pggb_output_dir      = PGGB_CONSTRUCT.out.output_dir         // [meta, pggb_output/]

    // Graph validation results
    validation_report    = GRAPH_VALIDATE.out.report            // [meta, validation_report.html]
    validation_stats     = GRAPH_VALIDATE.out.stats             // [meta, validation.json]
    validation_log       = GRAPH_VALIDATE.out.log               // [meta, validation.log]

    // Graph statistics and analysis
    graph_stats_json     = GRAPH_STATS.out.stats_json           // [meta, graph_stats.json]
    graph_stats_html     = GRAPH_STATS.out.stats_html           // [meta, graph_report.html]
    node_statistics      = GRAPH_STATS.out.node_stats           // [meta, node_stats.tsv]
    path_statistics      = GRAPH_STATS.out.path_stats           // [meta, path_stats.tsv]
    complexity_analysis  = GRAPH_STATS.out.complexity           // [meta, complexity.tsv]

    // Basic construction outputs
    construction_stats   = PGGB_CONSTRUCT.out.stats             // [meta, stats.txt]

    // Version information
    versions             = ch_versions                           // versions.yml
}

workflow.onError {
    log.error "Stage 3 Graph Construction workflow failed: ${workflow.errorMessage}"
    log.error "Common issues:"
    log.error "  - Insufficient memory for PGGB (requires 64GB+ for large genomes)"
    log.error "  - PGGB container not available or network issues"
    log.error "  - Input genomes not properly standardized"
    log.error "  - Graph complexity too high for available resources"
}

workflow.onComplete {
    log.info "Stage 3 Graph Construction completed successfully"
    log.info "Pangenome graph constructed and validated"
    log.info "Next stage: Variant calling and comparative analysis"
}