process GRAPH_STATS {
    tag "$meta.id"
    label 'process_medium'

    container "pangenome/odgi:latest"

    input:
    tuple val(meta), path(gfa_graph), path(odgi_graph)

    output:
    tuple val(meta), path("*_graph_stats.json")  , emit: stats_json
    tuple val(meta), path("*_graph_report.html") , emit: stats_html
    tuple val(meta), path("*_node_stats.tsv")    , emit: node_stats
    tuple val(meta), path("*_path_stats.tsv")    , emit: path_stats
    tuple val(meta), path("*_complexity.tsv")    , emit: complexity
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo "📊 Starting comprehensive graph statistics for ${meta.id}"

    # Initialize outputs
    stats_json="${prefix}_graph_stats.json"
    stats_html="${prefix}_graph_report.html"
    node_stats="${prefix}_node_stats.tsv"
    path_stats="${prefix}_path_stats.tsv"
    complexity_stats="${prefix}_complexity.tsv"

    # Basic ODGI statistics
    echo "🔧 Computing basic ODGI statistics..."
    odgi stats -i ${odgi_graph} -S > basic_stats.txt 2>&1

    # Extract key metrics
    nodes=\$(grep "nodes" basic_stats.txt | awk '{print \$2}' || echo "0")
    edges=\$(grep "edges" basic_stats.txt | awk '{print \$2}' || echo "0")
    paths=\$(grep "paths" basic_stats.txt | awk '{print \$2}' || echo "0")
    length=\$(grep "length" basic_stats.txt | awk '{print \$2}' || echo "0")

    echo "📈 Basic metrics: \$nodes nodes, \$edges edges, \$paths paths"

    # Node degree distribution
    echo "🕸️ Computing node degree distribution..."
    if odgi degree -i ${odgi_graph} -d > node_degrees.txt 2>/dev/null; then
        # Calculate degree statistics
        sort -n node_degrees.txt > sorted_degrees.txt

        min_degree=\$(head -1 sorted_degrees.txt)
        max_degree=\$(tail -1 sorted_degrees.txt)
        avg_degree=\$(awk '{sum+=\$1} END {printf "%.2f", sum/NR}' sorted_degrees.txt)
        median_degree=\$(awk '{a[NR]=\$1} END {print (NR%2==1) ? a[int(NR/2)+1] : (a[NR/2]+a[NR/2+1])/2}' sorted_degrees.txt)

        echo "Node degree stats: min=\$min_degree, max=\$max_degree, avg=\$avg_degree, median=\$median_degree"
    else
        min_degree="N/A"
        max_degree="N/A"
        avg_degree="N/A"
        median_degree="N/A"
        echo "⚠️ Could not compute node degree distribution"
    fi

    # Path statistics
    echo "🛤️ Computing path statistics..."
    if odgi paths -i ${odgi_graph} -L > path_lengths.txt 2>/dev/null; then
        # Create detailed path stats
        echo -e "path_name\\tpath_length\\tpath_nodes\\tstatus" > \$path_stats

        while read line; do
            path_name=\$(echo \$line | awk '{print \$1}')
            path_length=\$(echo \$line | awk '{print \$2}')

            # Get path node count
            if odgi paths -i ${odgi_graph} -p \$path_name -v > path_nodes_temp.txt 2>/dev/null; then
                path_nodes=\$(wc -l < path_nodes_temp.txt)
            else
                path_nodes="N/A"
            fi

            # Classify path status
            if [ \$path_length -gt 1000000 ]; then
                status="large_genome"
            elif [ \$path_length -gt 100000 ]; then
                status="medium_genome"
            else
                status="small_fragment"
            fi

            echo -e "\$path_name\\t\$path_length\\t\$path_nodes\\t\$status" >> \$path_stats
        done < path_lengths.txt

        # Path summary statistics
        min_path_length=\$(tail -n +2 \$path_stats | awk -F'\\t' '{print \$2}' | sort -n | head -1)
        max_path_length=\$(tail -n +2 \$path_stats | awk -F'\\t' '{print \$2}' | sort -n | tail -1)
        avg_path_length=\$(tail -n +2 \$path_stats | awk -F'\\t' '{sum+=\$2} END {printf "%.0f", sum/NR}')

    else
        echo -e "path_name\\tpath_length\\tpath_nodes\\tstatus\\nN/A\\t0\\t0\\terror" > \$path_stats
        min_path_length="0"
        max_path_length="0"
        avg_path_length="0"
        echo "⚠️ Could not compute path statistics"
    fi

    # Graph complexity analysis
    echo "🧮 Analyzing graph complexity..."
    echo -e "metric\\tvalue\\tinterpretation" > \$complexity_stats

    # Node density
    if [ \$nodes -gt 0 ] && [ \$length -gt 0 ]; then
        node_density=\$(echo "scale=6; \$nodes / \$length * 1000000" | bc -l 2>/dev/null || echo "N/A")
        echo -e "node_density_per_mb\\t\$node_density\\tnodes per megabase" >> \$complexity_stats
    else
        echo -e "node_density_per_mb\\tN/A\\tcannot calculate" >> \$complexity_stats
    fi

    # Edge/node ratio
    if [ \$nodes -gt 0 ] && [ \$edges -gt 0 ]; then
        edge_node_ratio=\$(echo "scale=3; \$edges / \$nodes" | bc -l 2>/dev/null || echo "N/A")
        echo -e "edge_node_ratio\\t\$edge_node_ratio\\taverage edges per node" >> \$complexity_stats

        # Complexity classification
        if [ "\$edge_node_ratio" != "N/A" ]; then
            complexity_class=\$(echo "\$edge_node_ratio > 2.0" | bc -l 2>/dev/null)
            if [ "\$complexity_class" = "1" ]; then
                complexity_level="high_complexity"
            else
                complexity_level="low_complexity"
            fi
        else
            complexity_level="unknown"
        fi
    else
        echo -e "edge_node_ratio\\tN/A\\tcannot calculate" >> \$complexity_stats
        complexity_level="unknown"
    fi

    echo -e "complexity_classification\\t\$complexity_level\\toverall graph complexity" >> \$complexity_stats

    # Bubbles and variants (if available)
    echo "🫧 Detecting graph bubbles..."
    if command -v vg &> /dev/null && [ -f "${gfa_graph}" ]; then
        # Convert GFA to vg format and detect bubbles
        if vg convert -g ${gfa_graph} > graph.vg 2>/dev/null; then
            bubble_count=\$(vg snarls graph.vg 2>/dev/null | wc -l || echo "0")
            echo "Detected \$bubble_count potential bubbles/variants"
        else
            bubble_count="N/A"
            echo "⚠️ Could not detect bubbles (vg conversion failed)"
        fi
    else
        bubble_count="N/A"
        echo "⚠️ Bubble detection not available (vg not found)"
    fi

    # Node statistics table
    echo "📋 Creating node statistics table..."
    echo -e "metric\\tvalue\\tdescription" > \$node_stats
    echo -e "total_nodes\\t\$nodes\\ttotal number of graph nodes" >> \$node_stats
    echo -e "total_edges\\t\$edges\\ttotal number of graph edges" >> \$node_stats
    echo -e "graph_length\\t\$length\\ttotal sequence length in graph" >> \$node_stats
    echo -e "min_node_degree\\t\$min_degree\\tminimum node connectivity" >> \$node_stats
    echo -e "max_node_degree\\t\$max_degree\\tmaximum node connectivity" >> \$node_stats
    echo -e "avg_node_degree\\t\$avg_degree\\taverage node connectivity" >> \$node_stats
    echo -e "median_node_degree\\t\$median_degree\\tmedian node connectivity" >> \$node_stats

    # Generate comprehensive JSON statistics
    echo "📊 Generating comprehensive statistics JSON..."
    cat > \$stats_json << EOF
{
    "graph_id": "${meta.id}",
    "analysis_timestamp": "\$(date -Iseconds)",
    "basic_statistics": {
        "nodes": \$nodes,
        "edges": \$edges,
        "paths": \$paths,
        "total_length": \$length,
        "edge_node_ratio": "\$edge_node_ratio"
    },
    "node_statistics": {
        "min_degree": "\$min_degree",
        "max_degree": "\$max_degree",
        "avg_degree": "\$avg_degree",
        "median_degree": "\$median_degree",
        "node_density_per_mb": "\$node_density"
    },
    "path_statistics": {
        "total_paths": \$paths,
        "min_path_length": \$min_path_length,
        "max_path_length": \$max_path_length,
        "avg_path_length": \$avg_path_length
    },
    "complexity_analysis": {
        "classification": "\$complexity_level",
        "edge_node_ratio": "\$edge_node_ratio",
        "estimated_bubbles": "\$bubble_count"
    },
    "file_information": {
        "gfa_size_bytes": "\$(stat -f%z "${gfa_graph}" 2>/dev/null || stat -c%s "${gfa_graph}" 2>/dev/null)",
        "odgi_size_bytes": "\$(stat -f%z "${odgi_graph}" 2>/dev/null || stat -c%s "${odgi_graph}" 2>/dev/null)",
        "gfa_size_human": "\$(du -h ${gfa_graph} | cut -f1)",
        "odgi_size_human": "\$(du -h ${odgi_graph} | cut -f1)"
    },
    "analysis_parameters": {
        "tools_used": ["odgi", "basic statistics"],
        "computation_host": "\$(hostname)",
        "analysis_duration_seconds": \$SECONDS
    }
}
EOF

    # Generate comprehensive HTML report
    echo "📄 Generating HTML statistics report..."
    cat > \$stats_html << 'EOF'
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Pangenome Graph Statistics Report</title>
    <style>
        body { font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; margin: 0; padding: 20px; background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); }
        .container { max-width: 1200px; margin: 0 auto; background: white; border-radius: 15px; overflow: hidden; box-shadow: 0 20px 40px rgba(0,0,0,0.1); }
        .header { background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; padding: 40px; text-align: center; }
        .header h1 { margin: 0; font-size: 2.5em; font-weight: 300; }
        .header p { margin: 10px 0 0 0; opacity: 0.9; }
        .content { padding: 40px; }
        .stats-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(300px, 1fr)); gap: 25px; margin: 30px 0; }
        .stat-card { background: #f8f9fa; border: 1px solid #e9ecef; border-radius: 10px; padding: 25px; transition: transform 0.2s; }
        .stat-card:hover { transform: translateY(-5px); box-shadow: 0 10px 25px rgba(0,0,0,0.1); }
        .stat-card h3 { margin: 0 0 15px 0; color: #2c3e50; font-size: 1.3em; }
        .stat-number { font-size: 2.5em; font-weight: bold; color: #3498db; margin: 10px 0; }
        .stat-label { color: #7f8c8d; font-size: 0.9em; text-transform: uppercase; letter-spacing: 1px; }
        .complexity-high { border-left: 5px solid #e74c3c; }
        .complexity-medium { border-left: 5px solid #f39c12; }
        .complexity-low { border-left: 5px solid #27ae60; }
        .table-container { background: white; border-radius: 10px; overflow: hidden; margin: 30px 0; box-shadow: 0 2px 10px rgba(0,0,0,0.1); }
        table { width: 100%; border-collapse: collapse; }
        th { background: #34495e; color: white; padding: 15px; text-align: left; font-weight: 600; }
        td { padding: 12px 15px; border-bottom: 1px solid #ecf0f1; }
        tr:hover { background: #f8f9fa; }
        .section { margin: 40px 0; }
        .section h2 { color: #2c3e50; font-size: 1.8em; margin-bottom: 20px; padding-bottom: 10px; border-bottom: 3px solid #3498db; }
        .status-indicator { display: inline-block; padding: 5px 15px; border-radius: 20px; font-size: 0.8em; font-weight: bold; text-transform: uppercase; }
        .status-high { background: #e74c3c; color: white; }
        .status-medium { background: #f39c12; color: white; }
        .status-low { background: #27ae60; color: white; }
        .metric-description { font-size: 0.9em; color: #7f8c8d; margin-top: 5px; }
        .chart-placeholder { background: #f8f9fa; border: 2px dashed #bdc3c7; border-radius: 10px; padding: 40px; text-align: center; color: #7f8c8d; margin: 20px 0; }
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>🧬 Pangenome Graph Statistics</h1>
            <p><strong>Graph ID:</strong> ${meta.id} | <strong>Generated:</strong> \$(date)</p>
        </div>

        <div class="content">
            <div class="section">
                <h2>📊 Overview Statistics</h2>
                <div class="stats-grid">
                    <div class="stat-card">
                        <h3>🔵 Graph Nodes</h3>
                        <div class="stat-number">\$nodes</div>
                        <div class="stat-label">Total Sequence Nodes</div>
                        <div class="metric-description">Each node represents a unique sequence segment in the pangenome</div>
                    </div>
                    <div class="stat-card">
                        <h3>🔗 Graph Edges</h3>
                        <div class="stat-number">\$edges</div>
                        <div class="stat-label">Node Connections</div>
                        <div class="metric-description">Edges connect nodes that appear sequentially in at least one genome</div>
                    </div>
                    <div class="stat-card">
                        <h3>🛤️ Graph Paths</h3>
                        <div class="stat-number">\$paths</div>
                        <div class="stat-label">Genome Sequences</div>
                        <div class="metric-description">Each path represents a complete genome sequence through the graph</div>
                    </div>
                    <div class="stat-card complexity-\$\$(echo \$complexity_level | sed 's/.*_//')">
                        <h3>🧮 Complexity</h3>
                        <div class="stat-number">\$edge_node_ratio</div>
                        <div class="stat-label">Edge/Node Ratio</div>
                        <div class="status-indicator status-\$\$(echo \$complexity_level | sed 's/.*_//')">\$\$(echo \$complexity_level | tr '_' ' ')</div>
                        <div class="metric-description">Higher ratios indicate more complex genomic variation</div>
                    </div>
                </div>
            </div>

            <div class="section">
                <h2>📏 Size and Scale Metrics</h2>
                <div class="table-container">
                    <table>
                        <thead>
                            <tr>
                                <th>Metric</th>
                                <th>Value</th>
                                <th>Description</th>
                            </tr>
                        </thead>
                        <tbody>
                            <tr>
                                <td><strong>Total Graph Length</strong></td>
                                <td>\$\$(echo \$length | awk '{printf "%.2f MB", \$1/1000000}')</td>
                                <td>Total sequence content across all paths</td>
                            </tr>
                            <tr>
                                <td><strong>Node Density</strong></td>
                                <td>\$node_density nodes/MB</td>
                                <td>Graph fragmentation level - higher values indicate more variation</td>
                            </tr>
                            <tr>
                                <td><strong>Average Path Length</strong></td>
                                <td>\$\$(echo \$avg_path_length | awk '{printf "%.2f MB", \$1/1000000}')</td>
                                <td>Mean genome size in the pangenome</td>
                            </tr>
                            <tr>
                                <td><strong>GFA File Size</strong></td>
                                <td>\$(du -h ${gfa_graph} | cut -f1)</td>
                                <td>Disk space used by graph representation</td>
                            </tr>
                            <tr>
                                <td><strong>ODGI File Size</strong></td>
                                <td>\$(du -h ${odgi_graph} | cut -f1)</td>
                                <td>Optimized graph format for analysis</td>
                            </tr>
                        </tbody>
                    </table>
                </div>
            </div>

            <div class="section">
                <h2>🕸️ Connectivity Analysis</h2>
                <div class="stats-grid">
                    <div class="stat-card">
                        <h3>📊 Min Node Degree</h3>
                        <div class="stat-number">\$min_degree</div>
                        <div class="stat-label">Minimum Connections</div>
                    </div>
                    <div class="stat-card">
                        <h3>📈 Max Node Degree</h3>
                        <div class="stat-number">\$max_degree</div>
                        <div class="stat-label">Maximum Connections</div>
                    </div>
                    <div class="stat-card">
                        <h3>📊 Average Degree</h3>
                        <div class="stat-number">\$avg_degree</div>
                        <div class="stat-label">Mean Connections</div>
                    </div>
                    <div class="stat-card">
                        <h3>🫧 Estimated Bubbles</h3>
                        <div class="stat-number">\$bubble_count</div>
                        <div class="stat-label">Variation Structures</div>
                    </div>
                </div>
            </div>

            <div class="section">
                <h2>🔍 Analysis Summary</h2>
                <div class="stat-card">
                    <h3>Graph Quality Assessment</h3>
EOF

    # Add quality assessment based on metrics
    if [ \$nodes -gt 1000 ] && [ \$edges -gt 1000 ] && [ \$paths -ge 3 ]; then
        echo '                    <div class="status-indicator status-low">HIGH QUALITY GRAPH</div>' >> \$stats_html
        echo '                    <div class="metric-description">Graph has sufficient complexity and represents multiple genomes well</div>' >> \$stats_html
    elif [ \$nodes -gt 100 ] && [ \$edges -gt 100 ] && [ \$paths -ge 2 ]; then
        echo '                    <div class="status-indicator status-medium">MODERATE QUALITY GRAPH</div>' >> \$stats_html
        echo '                    <div class="metric-description">Graph is functional but may benefit from parameter optimization</div>' >> \$stats_html
    else
        echo '                    <div class="status-indicator status-high">LOW COMPLEXITY GRAPH</div>' >> \$stats_html
        echo '                    <div class="metric-description">Graph may be over-simplified or construction parameters need adjustment</div>' >> \$stats_html
    fi

    cat >> \$stats_html << 'EOF'
                </div>
            </div>

            <div class="section">
                <p style="text-align: center; color: #7f8c8d; font-style: italic;">
                    Generated by PGGB Graph Statistics Pipeline • Analysis completed in \$SECONDS seconds
                </p>
            </div>
        </div>
    </div>
</body>
</html>
EOF

    echo "✅ Graph statistics analysis completed"
    echo "📊 Generated comprehensive statistics for \$nodes nodes, \$edges edges, \$paths paths"
    echo "🔧 Graph complexity: \$complexity_level (edge/node ratio: \$edge_node_ratio)"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        odgi: \$(odgi version 2>&1 | head -1 || echo "unknown")
        bc: \$(bc --version 2>&1 | head -1 | awk '{print \$2}' || echo "unknown")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_graph_stats.json
    touch ${prefix}_graph_report.html
    touch ${prefix}_node_stats.tsv
    touch ${prefix}_path_stats.tsv
    touch ${prefix}_complexity.tsv
    touch versions.yml
    """
}