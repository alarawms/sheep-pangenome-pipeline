/*
========================================================================================
    BUSCO Assessment Module
========================================================================================
    Description: Assess genome completeness using BUSCO (Mammalia lineage)
    Features: Complete gene space analysis for quality tier assignment
========================================================================================
*/

process BUSCO_ASSESSMENT {
    tag "$meta.id"
    label 'process_high'
    label 'process_long'

    container 'ezlabgva/busco:v5.4.7_cv1'

    publishDir "${params.outdir}/02_preprocessing/busco", mode: params.publish_dir_mode,
        saveAs: { filename ->
            if (filename.endsWith('_summary.json')) "${meta.id}_busco_summary.json"
            else if (filename.endsWith('_full_table.tsv')) "${meta.id}_busco_full_table.tsv"
            else if (filename.endsWith('_short_summary.txt')) "${meta.id}_busco_short_summary.txt"
            else null
        }

    input:
    tuple val(meta), path(genome)

    output:
    tuple val(meta), path("${meta.id}_busco_summary.json")      , emit: summary
    tuple val(meta), path("${meta.id}_busco_full_table.tsv")    , emit: full_table
    tuple val(meta), path("${meta.id}_busco_short_summary.txt") , emit: short_summary
    path "versions.yml"                                         , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''

    """
    # Set up BUSCO working directory
    mkdir -p busco_work
    cd busco_work

    # Download BUSCO Mammalia lineage dataset if not present
    if [ ! -d "mammalia_odb10" ]; then
        echo "Downloading BUSCO Mammalia lineage dataset..."
        busco --download mammalia
    fi

    # Run BUSCO analysis
    echo "Running BUSCO assessment for ${meta.id}..."
    busco \\
        --in ../${genome} \\
        --out ${meta.id}_busco \\
        --mode genome \\
        --lineage_dataset mammalia_odb10 \\
        --cpu ${task.cpus} \\
        --force \\
        --quiet \\
        ${args}

    # Move results to expected locations
    if [ -d "${meta.id}_busco" ]; then
        # Copy short summary
        if [ -f "${meta.id}_busco/short_summary.*.${meta.id}_busco.txt" ]; then
            cp ${meta.id}_busco/short_summary.*.${meta.id}_busco.txt ../${meta.id}_busco_short_summary.txt
        fi

        # Copy full table
        if [ -f "${meta.id}_busco/run_mammalia_odb10/full_table.tsv" ]; then
            cp ${meta.id}_busco/run_mammalia_odb10/full_table.tsv ../${meta.id}_busco_full_table.tsv
        fi
    else
        echo "Warning: BUSCO output directory not found"
        # Create empty files to prevent pipeline failure
        touch ../${meta.id}_busco_short_summary.txt
        touch ../${meta.id}_busco_full_table.tsv
    fi

    # Return to parent directory
    cd ..

    # Parse BUSCO results and create JSON summary
    python3 << 'EOF'
import json
import re
import os
import sys
from pathlib import Path

def parse_busco_short_summary(summary_file):
    """Parse BUSCO short summary file"""

    if not os.path.exists(summary_file):
        return {
            'error': 'BUSCO summary file not found',
            'complete_buscos': 0,
            'complete_buscos_percentage': 0.0,
            'complete_single_copy': 0,
            'complete_duplicated': 0,
            'fragmented_buscos': 0,
            'fragmented_buscos_percentage': 0.0,
            'missing_buscos': 0,
            'missing_buscos_percentage': 0.0,
            'total_buscos': 0,
            'lineage': 'unknown'
        }

    results = {}

    try:
        with open(summary_file, 'r') as f:
            content = f.read()

        # Extract key metrics using regex patterns
        patterns = {
            'complete_buscos': r'Complete BUSCOs.*?(\\d+)',
            'complete_buscos_percentage': r'Complete BUSCOs.*?\\((\\d+\\.\\d+)%\\)',
            'complete_single_copy': r'Complete and single-copy BUSCOs.*?(\\d+)',
            'complete_duplicated': r'Complete and duplicated BUSCOs.*?(\\d+)',
            'fragmented_buscos': r'Fragmented BUSCOs.*?(\\d+)',
            'fragmented_buscos_percentage': r'Fragmented BUSCOs.*?\\((\\d+\\.\\d+)%\\)',
            'missing_buscos': r'Missing BUSCOs.*?(\\d+)',
            'missing_buscos_percentage': r'Missing BUSCOs.*?\\((\\d+\\.\\d+)%\\)',
            'total_buscos': r'Total BUSCO groups searched.*?(\\d+)',
            'lineage': r'The lineage dataset is: (\\w+)'
        }

        for key, pattern in patterns.items():
            match = re.search(pattern, content, re.IGNORECASE)
            if match:
                value = match.group(1)
                if key.endswith('_percentage'):
                    results[key] = float(value)
                elif key == 'lineage':
                    results[key] = value
                else:
                    results[key] = int(value)
            else:
                # Set defaults for missing values
                if key.endswith('_percentage'):
                    results[key] = 0.0
                elif key == 'lineage':
                    results[key] = 'mammalia_odb10'
                else:
                    results[key] = 0

    except Exception as e:
        print(f"Error parsing BUSCO summary: {e}")
        results = {
            'error': f'Parsing error: {str(e)}',
            'complete_buscos': 0,
            'complete_buscos_percentage': 0.0,
            'complete_single_copy': 0,
            'complete_duplicated': 0,
            'fragmented_buscos': 0,
            'fragmented_buscos_percentage': 0.0,
            'missing_buscos': 0,
            'missing_buscos_percentage': 0.0,
            'total_buscos': 0,
            'lineage': 'mammalia_odb10'
        }

    return results

def parse_busco_full_table(table_file):
    """Parse BUSCO full table for additional insights"""

    if not os.path.exists(table_file):
        return {'error': 'BUSCO full table not found'}

    try:
        status_counts = {'Complete': 0, 'Duplicated': 0, 'Fragmented': 0, 'Missing': 0}
        total_genes = 0

        with open(table_file, 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue

                fields = line.strip().split('\\t')
                if len(fields) >= 2:
                    status = fields[1]
                    status_counts[status] = status_counts.get(status, 0) + 1
                    total_genes += 1

        return {
            'detailed_counts': status_counts,
            'total_genes_analyzed': total_genes,
            'analysis_complete': True
        }

    except Exception as e:
        return {'error': f'Full table parsing error: {str(e)}'}

def assess_genome_quality(busco_results):
    """Assess genome quality based on BUSCO results"""

    complete_percentage = busco_results.get('complete_buscos_percentage', 0)
    missing_percentage = busco_results.get('missing_buscos_percentage', 100)
    duplicated_count = busco_results.get('complete_duplicated', 0)
    total_buscos = busco_results.get('total_buscos', 1)

    # Quality tier assignment based on completeness
    if complete_percentage >= 95:
        quality_tier = 'Excellent'
        quality_note = 'Highly complete genome suitable for reference use'
    elif complete_percentage >= 90:
        quality_tier = 'Good'
        quality_note = 'Good quality genome with minor gaps'
    elif complete_percentage >= 80:
        quality_tier = 'Fair'
        quality_note = 'Adequate quality with some gene space gaps'
    elif complete_percentage >= 70:
        quality_tier = 'Poor'
        quality_note = 'Lower quality with significant gene space gaps'
    else:
        quality_tier = 'Very Poor'
        quality_note = 'Very incomplete genome, not recommended for analysis'

    # Check for potential issues
    issues = []
    if missing_percentage > 20:
        issues.append(f'High missing gene content ({missing_percentage:.1f}%)')

    if total_buscos > 0 and (duplicated_count / total_buscos) > 0.1:
        issues.append(f'High duplication rate ({duplicated_count/total_buscos*100:.1f}%)')

    if complete_percentage < 80:
        issues.append('Below standard completeness threshold')

    return {
        'quality_tier': quality_tier,
        'quality_note': quality_note,
        'potential_issues': issues,
        'recommended_for_pangenome': complete_percentage >= 80 and missing_percentage <= 15
    }

# Main processing
try:
    sample_id = "${meta.id}"
    summary_file = f"{sample_id}_busco_short_summary.txt"
    table_file = f"{sample_id}_busco_full_table.tsv"
    json_output = f"{sample_id}_busco_summary.json"

    print(f"Processing BUSCO results for {sample_id}")

    # Parse BUSCO results
    short_summary = parse_busco_short_summary(summary_file)
    full_table = parse_busco_full_table(table_file)

    # Assess quality
    quality_assessment = assess_genome_quality(short_summary)

    # Compile comprehensive results
    busco_summary = {
        'sample_id': sample_id,
        'analysis_timestamp': os.popen('date -Iseconds').read().strip(),
        'busco_version': 'v5.4.7',
        'lineage_dataset': short_summary.get('lineage', 'mammalia_odb10'),

        # Core BUSCO metrics
        'busco_metrics': short_summary,

        # Additional analysis
        'detailed_analysis': full_table,

        # Quality assessment
        'quality_assessment': quality_assessment,

        # Summary statistics
        'summary_statistics': {
            'completeness_score': short_summary.get('complete_buscos_percentage', 0),
            'fragmentation_rate': short_summary.get('fragmented_buscos_percentage', 0),
            'missing_rate': short_summary.get('missing_buscos_percentage', 0),
            'duplication_rate': (short_summary.get('complete_duplicated', 0) /
                               max(1, short_summary.get('total_buscos', 1))) * 100
        }
    }

    # Write JSON summary
    with open(json_output, 'w') as f:
        json.dump(busco_summary, f, indent=2)

    # Print summary
    print(f"")
    print(f"BUSCO Assessment Results for {sample_id}:")
    print(f"  Completeness: {short_summary.get('complete_buscos_percentage', 0):.1f}%")
    print(f"  Missing: {short_summary.get('missing_buscos_percentage', 0):.1f}%")
    print(f"  Fragmented: {short_summary.get('fragmented_buscos_percentage', 0):.1f}%")
    print(f"  Quality tier: {quality_assessment.get('quality_tier', 'Unknown')}")
    print(f"  Recommended for pangenome: {quality_assessment.get('recommended_for_pangenome', False)}")

    if quality_assessment.get('potential_issues'):
        print(f"  Issues identified:")
        for issue in quality_assessment['potential_issues']:
            print(f"    - {issue}")

    print(f"  Results: {json_output}")

except Exception as e:
    print(f"❌ Error processing BUSCO results for {sample_id}: {str(e)}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

EOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        busco: \$(busco --version | sed 's/BUSCO //')
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """
}