/*
========================================================================================
    Manual Reference Genome Selection Module
========================================================================================
    Description: Select user-specified reference genome for PGGB pangenome construction
    Input: User specifies reference sample via --reference_sample parameter
========================================================================================
*/

process SELECT_REFERENCE_MANUAL {
    label 'process_low'

    publishDir "${params.outdir}/02_preprocessing/reference_selection", mode: params.publish_dir_mode

    input:
    tuple val(meta), path(standardized_genomes), path(qc_stats)

    output:
    tuple val(selected_meta), path("selected_reference.fa"), emit: reference_genome
    path "reference_metadata.json",                          emit: reference_metadata
    path "reference_selection_report.json",                  emit: selection_report
    path "versions.yml",                                     emit: versions

    when:
    params.enable_reference_selection && (task.ext.when == null || task.ext.when)

    script:
    def reference_sample = params.reference_sample
    """
    #!/bin/bash

    echo "🔍 Manual Reference Selection Process"
    echo "===================================="

    # Check if reference sample is specified
    if [ -z "${reference_sample}" ]; then
        echo "❌ ERROR: Reference selection is enabled but no reference sample specified!"
        echo "Please provide --reference_sample parameter"
        echo ""
        echo "Available samples:"
        for genome in ${standardized_genomes}; do
            sample_name=\$(basename \$genome _standardized.fa)
            echo "  - \$sample_name"
        done
        echo ""
        echo "Example: --enable_reference_selection true --reference_sample rambouillet"
        echo ""
        echo "To disable reference selection entirely, use: --enable_reference_selection false"
        exit 1
    fi

    echo "🎯 Looking for reference sample: ${reference_sample}"

    # Find the specified reference genome
    reference_found=false
    selected_genome=""
    selected_stats=""

    for genome in ${standardized_genomes}; do
        sample_name=\$(basename \$genome _standardized.fa)
        if [ "\$sample_name" = "${reference_sample}" ]; then
            echo "✅ Found reference genome: \$genome"
            selected_genome=\$genome
            reference_found=true
            break
        fi
    done

    # Find corresponding stats file
    for stats in ${qc_stats}; do
        sample_name=\$(basename \$stats _extended_stats.json)
        if [ "\$sample_name" = "${reference_sample}" ]; then
            echo "✅ Found stats file: \$stats"
            selected_stats=\$stats
            break
        fi
    done

    if [ "\$reference_found" = false ]; then
        echo "❌ ERROR: Reference sample '${reference_sample}' not found!"
        echo ""
        echo "Available samples:"
        for genome in ${standardized_genomes}; do
            sample_name=\$(basename \$genome _standardized.fa)
            echo "  - \$sample_name"
        done
        exit 1
    fi

    # Copy the selected reference
    cp "\$selected_genome" selected_reference.fa
    echo "📋 Reference genome copied to: selected_reference.fa"

    # Create reference metadata
    echo "📊 Creating reference metadata..."

    # Get basic stats if available
    genome_size=0
    if [ -f "\$selected_stats" ]; then
        genome_size=\$(python3 -c "
import json
try:
    with open('\$selected_stats', 'r') as f:
        data = json.load(f)
    size = data.get('sample_info', {}).get('total_length', 0)
    print(int(size))
except:
    print(0)
")
    fi

    cat > reference_metadata.json << EOF
{
  "selected_reference": {
    "sample_id": "${reference_sample}",
    "selection_method": "manual",
    "selection_timestamp": "\$(date -Iseconds)",
    "genome_size": \$genome_size,
    "file_path": "selected_reference.fa"
  },
  "selection_criteria": {
    "method": "manual_user_selection",
    "user_specified": true,
    "automatic_scoring": false
  },
  "pangenome_suitability": "User_Selected"
}
EOF

    # Create selection report
    cat > reference_selection_report.json << EOF
{
  "selection_summary": {
    "method": "manual",
    "selected_reference": "${reference_sample}",
    "selection_timestamp": "\$(date -Iseconds)",
    "user_specified": true
  },
  "selected_reference": {
    "sample_id": "${reference_sample}",
    "selection_method": "manual_user_selection",
    "genome_size": \$genome_size,
    "notes": "Reference manually selected by user via --reference_sample parameter"
  },
  "available_samples": [
\$(for genome in ${standardized_genomes}; do
    sample_name=\$(basename \$genome _standardized.fa)
    echo "    \"\$sample_name\","
done | sed '\$s/,\$//')
  ],
  "selection_criteria": {
    "method": "manual",
    "rationale": "User-specified reference for pangenome construction"
  }
}
EOF

    echo "✅ Manual reference selection completed"
    echo "   Selected: ${reference_sample}"
    echo "   Reference file: selected_reference.fa"
    echo "   Metadata: reference_metadata.json"
    echo "   Report: reference_selection_report.json"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(bash --version | head -1 | sed 's/.*version //' | sed 's/ .*//')
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    touch selected_reference.fa
    touch reference_metadata.json
    touch reference_selection_report.json
    touch versions.yml
    """
}