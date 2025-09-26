process GRAPH_VALIDATE {
    tag "$meta.id"
    label 'process_medium'

    container "pangenome/odgi:latest"

    input:
    tuple val(meta), path(gfa_graph), path(odgi_graph)

    output:
    tuple val(meta), path("*_validation_report.html"), emit: report
    tuple val(meta), path("*_validation.json")     , emit: stats
    tuple val(meta), path("*_validation.log")      , emit: log
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo "🔍 Starting graph validation for ${meta.id}"

    # Initialize log
    echo "Graph Validation Log for ${meta.id}" > ${prefix}_validation.log
    echo "Started: \$(date)" >> ${prefix}_validation.log
    echo "=====================================" >> ${prefix}_validation.log

    # Basic GFA validation
    echo "📋 Validating GFA structure..." | tee -a ${prefix}_validation.log

    validation_status="PASS"
    error_count=0
    warning_count=0

    # Check GFA file integrity
    if [ ! -s "${gfa_graph}" ]; then
        echo "❌ ERROR: GFA file is empty or missing" | tee -a ${prefix}_validation.log
        validation_status="FAIL"
        error_count=\$((error_count + 1))
    else
        echo "✅ GFA file exists and is non-empty" | tee -a ${prefix}_validation.log

        # Check GFA format
        gfa_lines=\$(wc -l < ${gfa_graph})
        s_lines=\$(grep -c "^S" ${gfa_graph} || echo 0)
        l_lines=\$(grep -c "^L" ${gfa_graph} || echo 0)
        p_lines=\$(grep -c "^P" ${gfa_graph} || echo 0)

        echo "  Total lines: \$gfa_lines" | tee -a ${prefix}_validation.log
        echo "  Sequence lines (S): \$s_lines" | tee -a ${prefix}_validation.log
        echo "  Link lines (L): \$l_lines" | tee -a ${prefix}_validation.log
        echo "  Path lines (P): \$p_lines" | tee -a ${prefix}_validation.log

        if [ \$s_lines -eq 0 ]; then
            echo "⚠️ WARNING: No sequence lines found in GFA" | tee -a ${prefix}_validation.log
            warning_count=\$((warning_count + 1))
        fi

        if [ \$l_lines -eq 0 ]; then
            echo "⚠️ WARNING: No link lines found in GFA" | tee -a ${prefix}_validation.log
            warning_count=\$((warning_count + 1))
        fi

        if [ \$p_lines -eq 0 ]; then
            echo "⚠️ WARNING: No path lines found in GFA" | tee -a ${prefix}_validation.log
            warning_count=\$((warning_count + 1))
        fi
    fi

    # ODGI validation
    echo "" | tee -a ${prefix}_validation.log
    echo "🔧 Validating ODGI graph..." | tee -a ${prefix}_validation.log

    if [ ! -s "${odgi_graph}" ]; then
        echo "❌ ERROR: ODGI file is empty or missing" | tee -a ${prefix}_validation.log
        validation_status="FAIL"
        error_count=\$((error_count + 1))
    else
        echo "✅ ODGI file exists and is non-empty" | tee -a ${prefix}_validation.log

        # Test ODGI operations
        if odgi stats -i ${odgi_graph} -S > odgi_stats_temp.txt 2>&1; then
            echo "✅ ODGI graph is valid and readable" | tee -a ${prefix}_validation.log

            # Extract key stats
            nodes=\$(grep "nodes" odgi_stats_temp.txt | awk '{print \$2}' || echo "unknown")
            edges=\$(grep "edges" odgi_stats_temp.txt | awk '{print \$2}' || echo "unknown")
            paths=\$(grep "paths" odgi_stats_temp.txt | awk '{print \$2}' || echo "unknown")

            echo "  Nodes: \$nodes" | tee -a ${prefix}_validation.log
            echo "  Edges: \$edges" | tee -a ${prefix}_validation.log
            echo "  Paths: \$paths" | tee -a ${prefix}_validation.log

            # Basic sanity checks
            if [ "\$nodes" != "unknown" ] && [ \$nodes -lt 100 ]; then
                echo "⚠️ WARNING: Very few nodes (\$nodes) - graph may be too simplified" | tee -a ${prefix}_validation.log
                warning_count=\$((warning_count + 1))
            fi

        else
            echo "❌ ERROR: ODGI graph validation failed" | tee -a ${prefix}_validation.log
            validation_status="FAIL"
            error_count=\$((error_count + 1))
        fi
    fi

    # Graph connectivity analysis
    echo "" | tee -a ${prefix}_validation.log
    echo "🕸️ Analyzing graph connectivity..." | tee -a ${prefix}_validation.log

    if [ -s "${odgi_graph}" ] && odgi stats -i ${odgi_graph} -S > connectivity_stats.txt 2>&1; then
        # Check for graph components
        if odgi components -i ${odgi_graph} -C > components.txt 2>/dev/null; then
            component_count=\$(wc -l < components.txt)
            echo "  Graph components: \$component_count" | tee -a ${prefix}_validation.log

            if [ \$component_count -gt 1 ]; then
                echo "⚠️ WARNING: Graph has multiple components (\$component_count) - may indicate fragmentation" | tee -a ${prefix}_validation.log
                warning_count=\$((warning_count + 1))
            else
                echo "✅ Graph is well-connected (single component)" | tee -a ${prefix}_validation.log
            fi
        else
            echo "⚠️ WARNING: Could not analyze graph components" | tee -a ${prefix}_validation.log
            warning_count=\$((warning_count + 1))
        fi
    else
        echo "❌ ERROR: Cannot analyze graph connectivity" | tee -a ${prefix}_validation.log
        error_count=\$((error_count + 1))
    fi

    # Generate validation JSON
    echo "📊 Generating validation report..." | tee -a ${prefix}_validation.log

    cat > ${prefix}_validation.json << EOF
{
    "graph_id": "${meta.id}",
    "validation_timestamp": "\$(date -Iseconds)",
    "overall_status": "\$validation_status",
    "summary": {
        "errors": \$error_count,
        "warnings": \$warning_count,
        "status": "\$validation_status"
    },
    "gfa_validation": {
        "file_exists": \$([ -s "${gfa_graph}" ] && echo true || echo false),
        "file_size_bytes": \$(stat -f%z "${gfa_graph}" 2>/dev/null || stat -c%s "${gfa_graph}" 2>/dev/null || echo 0),
        "total_lines": \$gfa_lines,
        "sequence_lines": \$s_lines,
        "link_lines": \$l_lines,
        "path_lines": \$p_lines
    },
    "odgi_validation": {
        "file_exists": \$([ -s "${odgi_graph}" ] && echo true || echo false),
        "file_size_bytes": \$(stat -f%z "${odgi_graph}" 2>/dev/null || stat -c%s "${odgi_graph}" 2>/dev/null || echo 0),
        "graph_readable": \$(odgi stats -i ${odgi_graph} -S >/dev/null 2>&1 && echo true || echo false),
        "nodes": "\$nodes",
        "edges": "\$edges",
        "paths": "\$paths"
    },
    "validation_criteria": {
        "min_nodes": 100,
        "min_edges": 100,
        "min_paths": 2,
        "max_components": 1
    }
}
EOF

    # Generate HTML report
    cat > ${prefix}_validation_report.html << EOF
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Pangenome Graph Validation Report</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 40px; background-color: #f5f5f5; }
        .container { max-width: 1000px; margin: 0 auto; background: white; padding: 30px; border-radius: 10px; box-shadow: 0 2px 10px rgba(0,0,0,0.1); }
        .header { text-align: center; margin-bottom: 30px; padding-bottom: 20px; border-bottom: 2px solid #eee; }
        .status-pass { color: #27ae60; font-weight: bold; }
        .status-fail { color: #e74c3c; font-weight: bold; }
        .status-warning { color: #f39c12; font-weight: bold; }
        .metric-card { background: #f8f9fa; padding: 15px; margin: 10px 0; border-radius: 5px; border-left: 4px solid #007bff; }
        .error { border-left-color: #e74c3c; }
        .warning { border-left-color: #f39c12; }
        .success { border-left-color: #27ae60; }
        table { width: 100%; border-collapse: collapse; margin: 20px 0; }
        th, td { padding: 12px; text-align: left; border-bottom: 1px solid #ddd; }
        th { background-color: #f8f9fa; font-weight: bold; }
        .section { margin: 30px 0; }
        .section h3 { color: #2c3e50; border-bottom: 2px solid #3498db; padding-bottom: 10px; }
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>🧬 Pangenome Graph Validation Report</h1>
            <p><strong>Graph ID:</strong> ${meta.id}</p>
            <p><strong>Validation Date:</strong> \$(date)</p>
            <div class="status-\$(echo \$validation_status | tr '[:upper:]' '[:lower:]')">
                <h2>Overall Status: \$validation_status</h2>
            </div>
        </div>

        <div class="section">
            <h3>📊 Summary</h3>
            <div class="metric-card \$([ \$error_count -gt 0 ] && echo 'error' || echo 'success')">
                <strong>Errors:</strong> \$error_count
            </div>
            <div class="metric-card \$([ \$warning_count -gt 0 ] && echo 'warning' || echo 'success')">
                <strong>Warnings:</strong> \$warning_count
            </div>
        </div>

        <div class="section">
            <h3>📋 GFA File Validation</h3>
            <table>
                <tr><th>Metric</th><th>Value</th><th>Status</th></tr>
                <tr><td>File Size</td><td>\$(du -h ${gfa_graph} | cut -f1)</td><td>\$([ -s "${gfa_graph}" ] && echo '✅ OK' || echo '❌ FAIL')</td></tr>
                <tr><td>Total Lines</td><td>\$gfa_lines</td><td>✅ OK</td></tr>
                <tr><td>Sequence Lines (S)</td><td>\$s_lines</td><td>\$([ \$s_lines -gt 0 ] && echo '✅ OK' || echo '⚠️ WARNING')</td></tr>
                <tr><td>Link Lines (L)</td><td>\$l_lines</td><td>\$([ \$l_lines -gt 0 ] && echo '✅ OK' || echo '⚠️ WARNING')</td></tr>
                <tr><td>Path Lines (P)</td><td>\$p_lines</td><td>\$([ \$p_lines -gt 0 ] && echo '✅ OK' || echo '⚠️ WARNING')</td></tr>
            </table>
        </div>

        <div class="section">
            <h3>🔧 ODGI Graph Validation</h3>
            <table>
                <tr><th>Metric</th><th>Value</th><th>Status</th></tr>
                <tr><td>File Size</td><td>\$(du -h ${odgi_graph} 2>/dev/null | cut -f1 || echo 'N/A')</td><td>\$([ -s "${odgi_graph}" ] && echo '✅ OK' || echo '❌ FAIL')</td></tr>
                <tr><td>Graph Readable</td><td>\$(odgi stats -i ${odgi_graph} -S >/dev/null 2>&1 && echo 'Yes' || echo 'No')</td><td>\$(odgi stats -i ${odgi_graph} -S >/dev/null 2>&1 && echo '✅ OK' || echo '❌ FAIL')</td></tr>
                <tr><td>Nodes</td><td>\$nodes</td><td>\$([ "\$nodes" != "unknown" ] && [ \$nodes -gt 100 ] && echo '✅ OK' || echo '⚠️ WARNING')</td></tr>
                <tr><td>Edges</td><td>\$edges</td><td>\$([ "\$edges" != "unknown" ] && [ \$edges -gt 100 ] && echo '✅ OK' || echo '⚠️ WARNING')</td></tr>
                <tr><td>Paths</td><td>\$paths</td><td>\$([ "\$paths" != "unknown" ] && [ \$paths -gt 1 ] && echo '✅ OK' || echo '⚠️ WARNING')</td></tr>
            </table>
        </div>

        <div class="section">
            <h3>🕸️ Connectivity Analysis</h3>
            <div class="metric-card">
                <strong>Graph Components:</strong> \$([ -f "components.txt" ] && wc -l < components.txt || echo "Unable to determine")
                <br><em>Single component indicates well-connected graph</em>
            </div>
        </div>

        <div class="section">
            <h3>📝 Recommendations</h3>
            <ul>
EOF

    # Add recommendations based on validation results
    if [ \$error_count -gt 0 ]; then
        echo "                <li class='error'>🚨 Address \$error_count critical errors before proceeding</li>" >> ${prefix}_validation_report.html
    fi

    if [ \$warning_count -gt 0 ]; then
        echo "                <li class='warning'>⚠️ Review \$warning_count warnings for potential improvements</li>" >> ${prefix}_validation_report.html
    fi

    if [ \$s_lines -lt 100 ]; then
        echo "                <li class='warning'>Consider increasing graph complexity - very few sequence nodes detected</li>" >> ${prefix}_validation_report.html
    fi

    if [ "\$validation_status" = "PASS" ]; then
        echo "                <li class='success'>✅ Graph validation passed - ready for downstream analysis</li>" >> ${prefix}_validation_report.html
    fi

    cat >> ${prefix}_validation_report.html << 'EOF'
            </ul>
        </div>

        <div class="section">
            <p><em>Generated by PGGB Graph Validation Pipeline</em></p>
        </div>
    </div>
</body>
</html>
EOF

    echo "" | tee -a ${prefix}_validation.log
    echo "✅ Graph validation completed" | tee -a ${prefix}_validation.log
    echo "📊 Status: \$validation_status" | tee -a ${prefix}_validation.log
    echo "❌ Errors: \$error_count" | tee -a ${prefix}_validation.log
    echo "⚠️ Warnings: \$warning_count" | tee -a ${prefix}_validation.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        odgi: \$(odgi version 2>&1 | head -1 || echo "unknown")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_validation_report.html
    touch ${prefix}_validation.json
    touch ${prefix}_validation.log
    touch versions.yml
    """
}