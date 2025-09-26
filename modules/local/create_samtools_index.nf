/*
========================================================================================
    Create Samtools Index Module
========================================================================================
    Description: Generate samtools FAIDX indices for random sequence access
    Purpose: Enable efficient sequence retrieval for pangenome analysis tools
========================================================================================
*/

process CREATE_SAMTOOLS_INDEX {
    tag "$meta.id"
    label 'process_low'

    container 'biocontainers/samtools:1.17--h00cdaf9_0'

    publishDir "${params.outdir}/02_preprocessing/indices/samtools", mode: params.publish_dir_mode,
        saveAs: { filename ->
            if (filename.startsWith("${meta.id}")) filename
            else null
        }

    input:
    tuple val(meta), path(genome)

    output:
    tuple val(meta), path("${meta.id}_samtools_index")  , emit: index
    path "versions.yml"                                 , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''

    """
    # Create output directory
    mkdir -p ${prefix}_samtools_index

    # Copy genome file to index directory with consistent naming
    cp ${genome} ${prefix}_samtools_index/${prefix}.fa

    # Create samtools FAIDX index
    echo "Creating samtools FAIDX index for ${meta.id}..."
    samtools faidx ${prefix}_samtools_index/${prefix}.fa ${args}

    # Verify index creation
    if [ ! -f "${prefix}_samtools_index/${prefix}.fa.fai" ]; then
        echo "❌ Samtools index creation failed for ${meta.id}"
        exit 1
    fi

    # Generate additional statistics using samtools
    echo "Generating sequence statistics..."

    # Create sequence statistics from the FAIDX index
    python3 << 'EOF'
import json
import os

def parse_faidx_index(fai_file):
    """Parse samtools FAIDX index file to extract sequence statistics"""

    sequences = []
    total_length = 0

    try:
        with open(fai_file, 'r') as f:
            for line in f:
                fields = line.strip().split('\\t')
                if len(fields) >= 5:
                    seq_name = fields[0]
                    seq_length = int(fields[1])
                    offset = int(fields[2])
                    line_bases = int(fields[3])
                    line_width = int(fields[4])

                    sequences.append({
                        'name': seq_name,
                        'length': seq_length,
                        'offset': offset,
                        'line_bases': line_bases,
                        'line_width': line_width
                    })

                    total_length += seq_length

        return sequences, total_length

    except Exception as e:
        print(f"Error parsing FAIDX file: {e}")
        return [], 0

def generate_sequence_summary(sequences, total_length, sample_id):
    """Generate summary statistics from sequence information"""

    if not sequences:
        return {'error': 'No sequences found'}

    # Sort sequences by length
    sorted_seqs = sorted(sequences, key=lambda x: x['length'], reverse=True)

    # Calculate N50/N90 statistics
    target_50 = total_length * 0.5
    target_90 = total_length * 0.9

    cumulative = 0
    n50 = n90 = l50 = l90 = 0

    for i, seq in enumerate(sorted_seqs):
        cumulative += seq['length']

        if cumulative >= target_50 and n50 == 0:
            n50 = seq['length']
            l50 = i + 1

        if cumulative >= target_90 and n90 == 0:
            n90 = seq['length']
            l90 = i + 1
            break

    # Size categories
    very_large = sum(1 for s in sequences if s['length'] >= 50e6)
    large = sum(1 for s in sequences if 1e6 <= s['length'] < 50e6)
    medium = sum(1 for s in sequences if 10e3 <= s['length'] < 1e6)
    small = sum(1 for s in sequences if s['length'] < 10e3)

    summary = {
        'sample_id': sample_id,
        'total_sequences': len(sequences),
        'total_length': total_length,
        'longest_sequence': sorted_seqs[0]['length'] if sorted_seqs else 0,
        'shortest_sequence': sorted_seqs[-1]['length'] if sorted_seqs else 0,
        'mean_length': total_length / len(sequences) if sequences else 0,
        'n50': n50,
        'l50': l50,
        'n90': n90,
        'l90': l90,
        'size_distribution': {
            'very_large_scaffolds': very_large,
            'large_scaffolds': large,
            'medium_scaffolds': medium,
            'small_scaffolds': small
        },
        'top_10_sequences': [
            {'name': seq['name'], 'length': seq['length']}
            for seq in sorted_seqs[:10]
        ]
    }

    return summary

# Process the FAIDX index
try:
    sample_id = "${meta.id}"
    fai_file = f"${prefix}_samtools_index/${prefix}.fa.fai"

    print(f"Processing FAIDX index for {sample_id}")

    sequences, total_length = parse_faidx_index(fai_file)
    summary = generate_sequence_summary(sequences, total_length, sample_id)

    # Write summary to JSON
    with open(f"${prefix}_samtools_index/sequence_summary.json", 'w') as f:
        json.dump(summary, f, indent=2)

    print(f"✅ Sequence summary generated:")
    print(f"   Total sequences: {len(sequences):,}")
    print(f"   Total length: {total_length:,} bp ({total_length/1e9:.2f} Gb)")
    print(f"   N50: {summary.get('n50', 0):,} bp ({summary.get('n50', 0)/1e6:.1f} Mb)")
    print(f"   Longest sequence: {summary.get('longest_sequence', 0):,} bp")

except Exception as e:
    print(f"❌ Error generating sequence summary: {e}")

EOF

    # Create index metadata
    cat > ${prefix}_samtools_index/index_info.json << EOF
{
    "sample_id": "${meta.id}",
    "index_type": "samtools_faidx",
    "reference_file": "${prefix}.fa",
    "index_file": "${prefix}.fa.fai",
    "creation_timestamp": "\$(date -Iseconds)",
    "samtools_version": "\$(samtools --version | head -1 | cut -d' ' -f2)",
    "genome_size": \$(wc -c < ${prefix}_samtools_index/${prefix}.fa),
    "sequence_count": \$(wc -l < ${prefix}_samtools_index/${prefix}.fa.fai),
    "ready_for_access": true
}
EOF

    # Print summary
    echo "✅ Samtools FAIDX index created successfully for ${meta.id}"
    echo "   Index directory: ${prefix}_samtools_index/"
    echo "   Reference file: ${prefix}.fa"
    echo "   Index file: ${prefix}.fa.fai"
    echo "   Sequences indexed: \$(wc -l < ${prefix}_samtools_index/${prefix}.fa.fai)"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -1 | cut -d' ' -f2)
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """
}