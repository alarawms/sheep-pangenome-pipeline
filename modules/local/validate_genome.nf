process VALIDATE_GENOME {
    tag "$meta.id"
    label 'process_low'

    container "python:3.9-slim"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*_validation.json"), emit: validation
    tuple val(meta), path("*_stats.txt")     , emit: stats
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    # Create validation script
    cat > validate_genome.py << 'PYTHON_EOF'
import os
import json
from datetime import datetime
from collections import defaultdict

def analyze_fasta(filename):
    stats = {
        "total_sequences": 0,
        "total_length": 0,
        "gc_content": 0,
        "n_content": 0,
        "sequences": []
    }

    gc_count = 0
    n_count = 0
    current_seq = {"name": "", "length": 0}

    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_seq["name"]:
                    stats["sequences"].append(current_seq)
                stats["total_sequences"] += 1
                current_seq = {"name": line[1:], "length": 0}
            else:
                seq_len = len(line)
                current_seq["length"] += seq_len
                stats["total_length"] += seq_len

                # Count GC and N content
                gc_count += line.upper().count('G') + line.upper().count('C')
                n_count += line.upper().count('N')

    # Add last sequence
    if current_seq["name"]:
        stats["sequences"].append(current_seq)

    # Calculate percentages
    if stats["total_length"] > 0:
        stats["gc_content"] = (gc_count / stats["total_length"]) * 100
        stats["n_content"] = (n_count / stats["total_length"]) * 100

    return stats

# Analyze the genome
fasta_stats = analyze_fasta("${fasta}")

# Validation criteria for sheep genomes
validation = {
    "sample_id": "${prefix}",
    "validation_timestamp": datetime.now().isoformat(),
    "total_length": fasta_stats["total_length"],
    "total_sequences": fasta_stats["total_sequences"],
    "gc_content": fasta_stats["gc_content"],
    "n_content": fasta_stats["n_content"],
    "validation_status": "PENDING"
}

# Apply validation rules
validation_passed = True
validation_messages = []

# Size validation (sheep genomes ~2.6-3.0 Gb)
if fasta_stats["total_length"] < 2.4e9:
    validation_passed = False
    validation_messages.append("Genome too small (<2.4Gb)")
elif fasta_stats["total_length"] > 3.2e9:
    validation_passed = False
    validation_messages.append("Genome too large (>3.2Gb)")

# GC content validation (mammals ~40-45%)
if fasta_stats["gc_content"] < 35 or fasta_stats["gc_content"] > 50:
    validation_passed = False
    validation_messages.append(f"Unusual GC content: {fasta_stats['gc_content']:.1f}%")

# N content validation (should be <5%)
if fasta_stats["n_content"] > 5.0:
    validation_passed = False
    validation_messages.append(f"High N content: {fasta_stats['n_content']:.1f}%")

# Sequence count validation (reasonable number of contigs)
if fasta_stats["total_sequences"] > 50000:
    validation_passed = False
    validation_messages.append(f"Too many sequences: {fasta_stats['total_sequences']}")

validation["validation_status"] = "PASS" if validation_passed else "FAIL"
validation["validation_messages"] = validation_messages

# Save validation results
with open("${prefix}_validation.json", 'w') as f:
    json.dump(validation, f, indent=2)

# Save detailed statistics
with open("${prefix}_stats.txt", 'w') as f:
    f.write(f"Genome Statistics for ${prefix}\\n")
    f.write(f"{'='*50}\\n")
    f.write(f"Total sequences: {fasta_stats['total_sequences']:,}\\n")
    f.write(f"Total length: {fasta_stats['total_length']:,} bp\\n")
    f.write(f"GC content: {fasta_stats['gc_content']:.2f}%\\n")
    f.write(f"N content: {fasta_stats['n_content']:.2f}%\\n")
    f.write(f"Validation: {validation['validation_status']}\\n")

    if validation_messages:
        f.write(f"\\nValidation messages:\\n")
        for msg in validation_messages:
            f.write(f"  - {msg}\\n")

    f.write(f"\\nLargest sequences:\\n")
    sorted_seqs = sorted(fasta_stats['sequences'], key=lambda x: x['length'], reverse=True)
    for i, seq in enumerate(sorted_seqs[:10]):
        f.write(f"  {i+1:2d}. {seq['name'][:50]:50s} {seq['length']:>12,} bp\\n")

PYTHON_EOF

    # Run the validation script
    python3 validate_genome.py

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | cut -d' ' -f2)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_validation.json
    touch ${prefix}_stats.txt
    touch versions.yml
    """
}