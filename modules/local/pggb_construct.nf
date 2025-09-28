process PGGB_CONSTRUCT {
    tag "pangenome_graph"
    label 'process_high'

    container "pangenome/pggb:latest"

    input:
    tuple val(meta), path(genomes)

    output:
    tuple val(meta), path("*.gfa")           , emit: graph
    tuple val(meta), path("*.og")            , emit: odgi_graph
    tuple val(meta), path("*.smooth.gfa")    , emit: smoothed_graph
    tuple val(meta), path("*.viz.png")       , emit: visualization
    tuple val(meta), path("*.stats.txt")     , emit: stats
    tuple val(meta), path("pggb_output/")    , emit: output_dir
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "sheep_pangenome"

    // PGGB parameters optimized for sheep genomes
    def segment_length = params.pggb_segment_length ?: 10000
    def block_length = params.pggb_block_length ?: 5000
    def identity_threshold = params.pggb_identity_threshold ?: 95
    def map_pct_id = params.pggb_map_pct_id ?: 95
    def n_mappings = params.pggb_n_mappings ?: 10
    def wfmash_params = params.pggb_wfmash_params ?: '--hg-filter-ani-diff 30 --merge-segments'
    def smoothxg_ratio = params.smoothxg_ratio ?: 10
    def threads = task.cpus
    def memory = "${task.memory.toGiga()}g"

    // Count number of genomes to set haplotype count
    def genome_count = genomes.size()
    def n_haplotypes = genome_count  // One haplotype per genome

    """
    echo "🧬 Starting PGGB pangenome graph construction"
    echo "📊 Parameters: segment=${segment_length}, block=${block_length}, identity=${identity_threshold}%"

    # Using de novo construction
    echo "📋 Using de novo pangenome construction"
    REFERENCE_GUIDED=false

    # Create output directory
    mkdir -p pggb_output

    # Combine all genomes into single input file for PGGB
    echo "📋 Preparing genome input file"
    rm -f combined_genomes.fa genomes.fa.gz
    for genome in ${genomes}; do
        echo "Processing \$genome"
        # Extract sample name from filename and add to headers
        sample_name=\$(basename \$genome .fa | sed 's/_standardized//')

        # Use awk for faster header modification on large files
        awk -v prefix="\$sample_name#" '/^>/ {print ">" prefix substr(\$0,2); next} {print}' \$genome >> combined_genomes.fa
    done

    # Compress input for PGGB (use gzip if bgzip not available)
    if command -v bgzip >/dev/null 2>&1; then
        bgzip -c combined_genomes.fa > genomes.fa.gz
        if command -v samtools >/dev/null 2>&1; then
            samtools faidx genomes.fa.gz
        fi
    else
        gzip -c combined_genomes.fa > genomes.fa.gz
    fi

    # Run PGGB with sheep-optimized parameters
    echo "🔧 Running PGGB graph construction with sheep-optimized parameters"
    echo "   Segment length: ${segment_length}bp, Block length: ${block_length}bp"
    echo "   Identity threshold: ${identity_threshold}%, Mappings: ${n_mappings}"

    # Run PGGB with basic parameters
    pggb \\
        -i genomes.fa.gz \\
        -o pggb_output \\
        -t ${threads} \\
        -p ${identity_threshold} \\
        -s ${segment_length} \\
        -l ${block_length} \\
        -n ${n_haplotypes}

    # Move and rename key outputs
    echo "📦 Organizing outputs"
    find pggb_output -name "*.final.gfa" -exec cp {} ${prefix}.gfa \\;
    find pggb_output -name "*.og" -exec cp {} ${prefix}.og \\;
    find pggb_output -name "*.smooth.final.gfa" -exec cp {} ${prefix}.smooth.gfa \\;
    find pggb_output -name "*.viz.png" -exec cp {} ${prefix}.viz.png \\;

    # Generate comprehensive statistics
    echo "📊 Generating graph statistics"
    if [ -f "${prefix}.gfa" ]; then
        echo "Graph Statistics for ${prefix}" > ${prefix}.stats.txt
        echo "================================" >> ${prefix}.stats.txt
        echo "Generated: \$(date)" >> ${prefix}.stats.txt
        echo "" >> ${prefix}.stats.txt

        # Basic GFA stats
        echo "GFA File Statistics:" >> ${prefix}.stats.txt
        echo "  File size: \$(du -h ${prefix}.gfa | cut -f1)" >> ${prefix}.stats.txt
        echo "  Sequences: \$(grep '^S' ${prefix}.gfa | wc -l)" >> ${prefix}.stats.txt
        echo "  Links: \$(grep '^L' ${prefix}.gfa | wc -l)" >> ${prefix}.stats.txt
        echo "  Paths: \$(grep '^P' ${prefix}.gfa | wc -l)" >> ${prefix}.stats.txt
        echo "" >> ${prefix}.stats.txt

        # ODGI stats if available
        if [ -f "${prefix}.og" ]; then
            echo "ODGI Graph Statistics:" >> ${prefix}.stats.txt
            odgi stats -i ${prefix}.og -S >> ${prefix}.stats.txt 2>/dev/null || echo "  ODGI stats not available" >> ${prefix}.stats.txt
        fi

        echo "" >> ${prefix}.stats.txt
        echo "PGGB Parameters Used:" >> ${prefix}.stats.txt
        echo "  Segment length: ${segment_length}" >> ${prefix}.stats.txt
        echo "  Block length: ${block_length}" >> ${prefix}.stats.txt
        echo "  Identity threshold: ${identity_threshold}%" >> ${prefix}.stats.txt
        echo "  Threads: ${threads}" >> ${prefix}.stats.txt
        echo "  Memory: ${memory}" >> ${prefix}.stats.txt

    else
        echo "ERROR: PGGB failed to generate graph output" > ${prefix}.stats.txt
        echo "Check pggb_output directory for detailed logs" >> ${prefix}.stats.txt
    fi

    echo "✅ PGGB pangenome graph construction completed"
    echo "📈 Graph file: ${prefix}.gfa"
    echo "🎨 Visualization: ${prefix}.viz.png"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pggb: \$(pggb --version 2>&1 | head -1 | sed 's/pggb //')
        odgi: \$(odgi version 2>&1 | head -1 || echo "unknown")
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
        bgzip: \$(bgzip --version 2>&1 | head -1 || echo "unknown")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "sheep_pangenome"
    """
    touch ${prefix}.gfa
    touch ${prefix}.og
    touch ${prefix}.smooth.gfa
    touch ${prefix}.viz.png
    touch ${prefix}.stats.txt
    mkdir -p pggb_output
    touch pggb_output/stub_output.txt
    touch versions.yml
    """
}