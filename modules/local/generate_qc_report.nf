/*
========================================================================================
    Generate QC Report Module
========================================================================================
    Description: Create comprehensive HTML QC reports for genome preprocessing
    Features: Interactive plots, quality metrics, and recommendations
========================================================================================
*/

process GENERATE_QC_REPORT {
    tag "$meta.id"
    label 'process_low'

    container 'biocontainers/python:3.9--1'

    publishDir "${params.outdir}/02_preprocessing/qc_reports", mode: params.publish_dir_mode,
        saveAs: { filename ->
            if (filename.endsWith('.html')) "${meta.id}_qc_report.html"
            else null
        }

    input:
    tuple val(meta), path(extended_stats)
    tuple val(meta), path(busco_results)

    output:
    tuple val(meta), path("${meta.id}_qc_report.html")  , emit: qc_report
    path "versions.yml"                                 , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    #!/usr/bin/env python3

    import json
    import sys
    from pathlib import Path
    from datetime import datetime

    def load_json_safely(json_file):
        \"\"\"Safely load JSON file with error handling\"\"\"
        try:
            if Path(json_file).exists():
                with open(json_file, 'r') as f:
                    return json.load(f)
            else:
                print(f"Warning: File not found: {json_file}")
                return {}
        except Exception as e:
            print(f"Error loading {json_file}: {e}")
            return {}

    def format_number(num, unit=''):
        \"\"\"Format numbers with appropriate units\"\"\"
        if num >= 1e9:
            return f"{num/1e9:.2f} G{unit}"
        elif num >= 1e6:
            return f"{num/1e6:.2f} M{unit}"
        elif num >= 1e3:
            return f"{num/1e3:.2f} K{unit}"
        else:
            return f"{num:.0f} {unit}"

    def get_quality_color(score, thresholds=[60, 80, 95]):
        \"\"\"Get color based on quality score\"\"\"
        if score >= thresholds[2]:
            return '#28a745'  # Green
        elif score >= thresholds[1]:
            return '#17a2b8'  # Blue
        elif score >= thresholds[0]:
            return '#ffc107'  # Yellow
        else:
            return '#dc3545'  # Red

    def create_progress_bar(value, max_value=100, label='', color='#007bff'):
        \"\"\"Create HTML progress bar\"\"\"
        percentage = min(100, (value / max_value) * 100) if max_value > 0 else 0
        return f'''
        <div class="progress-container">
            <div class="progress-label">{label}: {value:.1f}%</div>
            <div class="progress">
                <div class="progress-bar" style="width: {percentage}%; background-color: {color}"></div>
            </div>
        </div>
        '''

    def create_metric_card(title, value, description, color='#007bff'):
        \"\"\"Create metric card HTML\"\"\"
        return f'''
        <div class="metric-card">
            <div class="metric-header" style="border-left: 4px solid {color}">
                <h4>{title}</h4>
            </div>
            <div class="metric-value">{value}</div>
            <div class="metric-description">{description}</div>
        </div>
        '''

    def generate_qc_html_report(sample_id, extended_stats, busco_results):
        \"\"\"Generate comprehensive HTML QC report\"\"\"

        # Extract key metrics
        sample_info = extended_stats.get('sample_info', {})
        assembly_metrics = extended_stats.get('assembly_metrics', {})
        composition = extended_stats.get('nucleotide_composition', {})
        quality_assessment = extended_stats.get('quality_assessment', {})
        chromosome_candidates = extended_stats.get('chromosome_candidates', [])

        # BUSCO metrics
        busco_metrics = busco_results.get('busco_metrics', {})
        busco_quality = busco_results.get('quality_assessment', {})

        # Calculate overall quality score
        overall_score = 0
        if quality_assessment:
            scores = [
                quality_assessment.get('chromosome_completeness', 0),
                quality_assessment.get('contiguity_score', 0),
                quality_assessment.get('size_accuracy', 0),
                quality_assessment.get('gap_quality', 0)
            ]
            overall_score = sum(scores) / len(scores) if scores else 0

        # Add BUSCO score if available
        if busco_metrics.get('complete_buscos_percentage'):
            overall_score = (overall_score + busco_metrics['complete_buscos_percentage']) / 2

        html_content = f'''
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Genome QC Report - {sample_id}</title>
    <style>
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, 'Helvetica Neue', Arial, sans-serif;
            margin: 0;
            padding: 20px;
            background-color: #f8f9fa;
            color: #333;
        }}
        .container {{
            max-width: 1200px;
            margin: 0 auto;
            background: white;
            padding: 30px;
            border-radius: 8px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
        }}
        .header {{
            text-align: center;
            margin-bottom: 40px;
            padding-bottom: 20px;
            border-bottom: 2px solid #e9ecef;
        }}
        .header h1 {{
            color: #2c3e50;
            margin: 0;
            font-size: 2.5em;
        }}
        .header .sample-id {{
            color: #6c757d;
            font-size: 1.2em;
            margin-top: 10px;
        }}
        .section {{
            margin: 30px 0;
        }}
        .section h2 {{
            color: #2c3e50;
            border-bottom: 2px solid #007bff;
            padding-bottom: 10px;
            margin-bottom: 20px;
        }}
        .metrics-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(280px, 1fr));
            gap: 20px;
            margin: 20px 0;
        }}
        .metric-card {{
            background: #f8f9fa;
            padding: 20px;
            border-radius: 8px;
            border-left: 4px solid #007bff;
        }}
        .metric-header h4 {{
            margin: 0 0 10px 0;
            color: #2c3e50;
        }}
        .metric-value {{
            font-size: 1.8em;
            font-weight: bold;
            color: #007bff;
            margin: 10px 0;
        }}
        .metric-description {{
            color: #6c757d;
            font-size: 0.9em;
        }}
        .progress-container {{
            margin: 15px 0;
        }}
        .progress-label {{
            margin-bottom: 8px;
            font-weight: 500;
            color: #495057;
        }}
        .progress {{
            background-color: #e9ecef;
            border-radius: 4px;
            height: 20px;
            overflow: hidden;
        }}
        .progress-bar {{
            height: 100%;
            transition: width 0.3s ease;
            border-radius: 4px;
        }}
        .quality-score {{
            text-align: center;
            padding: 30px;
            margin: 20px 0;
            border-radius: 8px;
            color: white;
        }}
        .table {{
            width: 100%;
            border-collapse: collapse;
            margin: 20px 0;
        }}
        .table th, .table td {{
            padding: 12px;
            text-align: left;
            border-bottom: 1px solid #dee2e6;
        }}
        .table th {{
            background-color: #f8f9fa;
            font-weight: 600;
            color: #495057;
        }}
        .recommendations {{
            background-color: #d4edda;
            border: 1px solid #c3e6cb;
            border-radius: 5px;
            padding: 20px;
            margin: 20px 0;
        }}
        .recommendations h4 {{
            color: #155724;
            margin-top: 0;
        }}
        .recommendations ul {{
            margin: 10px 0;
            padding-left: 20px;
        }}
        .recommendations li {{
            margin: 5px 0;
            color: #155724;
        }}
        .warning {{
            background-color: #fff3cd;
            border: 1px solid #ffeaa7;
            color: #856404;
        }}
        .warning h4 {{
            color: #856404;
        }}
        .warning li {{
            color: #856404;
        }}
        .footer {{
            margin-top: 40px;
            padding-top: 20px;
            border-top: 1px solid #e9ecef;
            text-align: center;
            color: #6c757d;
            font-size: 0.9em;
        }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>🧬 Genome Quality Control Report</h1>
            <div class="sample-id">Sample: <strong>{sample_id}</strong></div>
            <div class="sample-id">Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</div>
        </div>

        <!-- Overall Quality Score -->
        <div class="section">
            <div class="quality-score" style="background-color: {get_quality_color(overall_score)}">
                <h2 style="margin: 0; color: white;">Overall Quality Score</h2>
                <div style="font-size: 3em; margin: 10px 0;">{overall_score:.1f}/100</div>
                <div style="font-size: 1.2em;">
                    {'Excellent' if overall_score >= 90 else 'Good' if overall_score >= 80 else 'Fair' if overall_score >= 70 else 'Poor' if overall_score >= 60 else 'Very Poor'}
                </div>
            </div>
        </div>

        <!-- Basic Assembly Metrics -->
        <div class="section">
            <h2>📊 Assembly Statistics</h2>
            <div class="metrics-grid">
                {create_metric_card(
                    "Total Length",
                    format_number(sample_info.get('total_length', 0), 'bp'),
                    "Total genome size",
                    get_quality_color(quality_assessment.get('size_accuracy', 0))
                )}
                {create_metric_card(
                    "Number of Sequences",
                    f"{sample_info.get('total_sequences', 0):,}",
                    "Total scaffolds/contigs",
                    '#6f42c1'
                )}
                {create_metric_card(
                    "Scaffold N50",
                    format_number(assembly_metrics.get('scaffold_n50', 0), 'bp'),
                    "50% of assembly in scaffolds this size or larger",
                    get_quality_color(quality_assessment.get('contiguity_score', 0))
                )}
                {create_metric_card(
                    "Largest Scaffold",
                    format_number(sample_info.get('max_sequence_length', 0), 'bp'),
                    "Size of longest scaffold",
                    '#20c997'
                )}
            </div>
        </div>

        <!-- Quality Metrics -->
        <div class="section">
            <h2>🔬 Quality Assessment</h2>
            <div class="metrics-grid">
                <div>
                    {create_progress_bar(
                        quality_assessment.get('chromosome_completeness', 0),
                        100,
                        'Chromosome Completeness',
                        get_quality_color(quality_assessment.get('chromosome_completeness', 0))
                    )}
                    {create_progress_bar(
                        quality_assessment.get('contiguity_score', 0),
                        100,
                        'Contiguity Score',
                        get_quality_color(quality_assessment.get('contiguity_score', 0))
                    )}
                </div>
                <div>
                    {create_progress_bar(
                        quality_assessment.get('size_accuracy', 0),
                        100,
                        'Size Accuracy',
                        get_quality_color(quality_assessment.get('size_accuracy', 0))
                    )}
                    {create_progress_bar(
                        quality_assessment.get('gap_quality', 0),
                        100,
                        'Gap Quality',
                        get_quality_color(quality_assessment.get('gap_quality', 0))
                    )}
                </div>
            </div>
        </div>

        <!-- Nucleotide Composition -->
        <div class="section">
            <h2>🧪 Nucleotide Composition</h2>
            <div class="metrics-grid">
                {create_metric_card(
                    "GC Content",
                    f"{composition.get('gc_percentage', 0):.2f}%",
                    "Guanine-Cytosine content",
                    '#17a2b8'
                )}
                {create_metric_card(
                    "AT Content",
                    f"{composition.get('at_percentage', 0):.2f}%",
                    "Adenine-Thymine content",
                    '#28a745'
                )}
                {create_metric_card(
                    "Gap Content",
                    f"{composition.get('n_percentage', 0):.2f}%",
                    "Unknown bases (N's)",
                    get_quality_color(100 - composition.get('n_percentage', 5) * 20)
                )}
            </div>
        </div>

        <!-- BUSCO Results (if available) -->
        '''

        if busco_metrics:
            html_content += f'''
        <div class="section">
            <h2>🎯 BUSCO Completeness Assessment</h2>
            <div class="metrics-grid">
                <div>
                    {create_progress_bar(
                        busco_metrics.get('complete_buscos_percentage', 0),
                        100,
                        'Complete BUSCOs',
                        get_quality_color(busco_metrics.get('complete_buscos_percentage', 0), [70, 90, 95])
                    )}
                    {create_progress_bar(
                        busco_metrics.get('fragmented_buscos_percentage', 0),
                        100,
                        'Fragmented BUSCOs',
                        '#ffc107'
                    )}
                </div>
                <div>
                    {create_progress_bar(
                        busco_metrics.get('missing_buscos_percentage', 0),
                        100,
                        'Missing BUSCOs',
                        '#dc3545'
                    )}
                    {create_metric_card(
                        "Total BUSCO Groups",
                        f"{busco_metrics.get('total_buscos', 0)}",
                        "From Mammalia lineage dataset",
                        '#6f42c1'
                    )}
                </div>
            </div>
        </div>
            '''

        # Chromosome candidates table
        if chromosome_candidates:
            html_content += '''
        <div class="section">
            <h2>🧬 Top Chromosome Candidates</h2>
            <table class="table">
                <thead>
                    <tr>
                        <th>Rank</th>
                        <th>Sequence Name</th>
                        <th>Length</th>
                        <th>Predicted Chromosome</th>
                    </tr>
                </thead>
                <tbody>
            '''

            for i, candidate in enumerate(chromosome_candidates[:10], 1):
                seq_name = candidate.get('header', 'Unknown')
                if len(seq_name) > 50:
                    seq_name = seq_name[:47] + '...'

                html_content += f'''
                    <tr>
                        <td>{i}</td>
                        <td>{seq_name}</td>
                        <td>{format_number(candidate.get('length', 0), 'bp')}</td>
                        <td>{candidate.get('predicted_chromosome', 'Unknown')}</td>
                    </tr>
                '''

            html_content += '''
                </tbody>
            </table>
        </div>
            '''

        # Recommendations
        recommendations = []
        warnings = []

        if overall_score >= 90:
            recommendations.append("Excellent quality genome suitable for reference use")
            recommendations.append("Ready for high-quality pangenome construction")
        elif overall_score >= 80:
            recommendations.append("Good quality genome suitable for pangenome inclusion")
            recommendations.append("Consider for population-level analyses")
        elif overall_score >= 70:
            recommendations.append("Adequate quality for basic genomic analyses")
            warnings.append("Consider quality limitations for reference applications")
        else:
            warnings.append("Quality concerns identified - review before inclusion")
            warnings.append("May require additional processing or exclusion")

        if composition.get('n_percentage', 0) > 3:
            warnings.append(f"High gap content ({composition.get('n_percentage', 0):.1f}%) may affect analysis quality")

        if assembly_metrics.get('scaffold_n50', 0) < 1e6:
            warnings.append("Low contiguity (N50 < 1Mb) may limit analysis capabilities")

        if busco_metrics.get('complete_buscos_percentage', 100) < 90:
            warnings.append(f"Lower BUSCO completeness ({busco_metrics.get('complete_buscos_percentage', 0):.1f}%) indicates potential gene space gaps")

        # Add recommendations section
        if recommendations:
            html_content += f'''
        <div class="section">
            <div class="recommendations">
                <h4>✅ Recommendations</h4>
                <ul>
                    {''.join(f'<li>{rec}</li>' for rec in recommendations)}
                </ul>
            </div>
        </div>
            '''

        if warnings:
            html_content += f'''
        <div class="section">
            <div class="recommendations warning">
                <h4>⚠️ Quality Considerations</h4>
                <ul>
                    {''.join(f'<li>{warn}</li>' for warn in warnings)}
                </ul>
            </div>
        </div>
            '''

        html_content += f'''
        <div class="footer">
            <p>Generated by Sheep Pangenome Pipeline Stage 2 - Genome Preprocessing & Indexing</p>
            <p>Report generated on {datetime.now().strftime('%Y-%m-%d at %H:%M:%S')}</p>
        </div>
    </div>
</body>
</html>
        '''

        return html_content

    # Main processing
    try:
        sample_id = "${meta.id}"
        stats_file = "${extended_stats}"
        busco_file = "${busco_results}"
        output_file = f"{sample_id}_qc_report.html"

        print(f"Generating QC report for {sample_id}")

        # Load data
        extended_stats = load_json_safely(stats_file)
        busco_results = load_json_safely(busco_file)

        # Generate HTML report
        html_content = generate_qc_html_report(sample_id, extended_stats, busco_results)

        # Write report
        with open(output_file, 'w') as f:
            f.write(html_content)

        print(f"✅ QC report generated: {output_file}")

    except Exception as e:
        print(f"❌ Error generating QC report for {sample_id}: {str(e)}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """
}