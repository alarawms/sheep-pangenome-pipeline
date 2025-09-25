process DOWNLOAD_GENOME {
    tag "$meta.id"
    label 'process_medium'

    container "quay.io/biocontainers/ncbi-datasets-cli:16.9.0--h2d05d04_0"

    input:
    tuple val(meta), val(accession)

    output:
    tuple val(meta), path("*.fa")          , emit: fasta
    tuple val(meta), path("*.json")        , emit: metadata
    tuple val(meta), path("download_log.txt"), emit: log
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def accession_clean = accession.trim()

    """
    # Download genome using NCBI datasets
    datasets download genome accession ${accession_clean} \\
        --include genome \\
        --filename ${accession_clean}.zip \\
        ${args}

    # Extract and process files
    unzip -q ${accession_clean}.zip

    # Find and rename genome file
    find . -name "*.fna" -o -name "*.fa" | head -1 | xargs -I {} cp {} ${prefix}.fa

    # Extract metadata
    find . -name "*.jsonl" | head -1 | xargs -I {} cp {} ${prefix}_temp.json || echo '{}' > ${prefix}_temp.json

    # Create structured metadata
    python3 << 'EOF'
import json
import os
from pathlib import Path

# Basic genome validation
fasta_file = "${prefix}.fa"
if os.path.exists(fasta_file):
    seq_count = 0
    seq_length = 0
    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                seq_count += 1
            else:
                seq_length += len(line.strip())

    # Create metadata
    metadata = {
        "sample_id": "${prefix}",
        "accession": "${accession_clean}",
        "download_success": True,
        "sequence_count": seq_count,
        "genome_length": seq_length,
        "file_path": f"{fasta_file}"
    }

    # Try to merge with NCBI metadata if available
    try:
        with open("${prefix}_temp.json", 'r') as f:
            ncbi_meta = json.load(f)
        metadata.update(ncbi_meta)
    except:
        pass

    with open("${prefix}.json", 'w') as f:
        json.dump(metadata, f, indent=2)
else:
    metadata = {
        "sample_id": "${prefix}",
        "accession": "${accession_clean}",
        "download_success": False,
        "error": "Genome file not found after download"
    }
    with open("${prefix}.json", 'w') as f:
        json.dump(metadata, f, indent=2)
EOF

    # Create download log
    echo "Downloaded: ${accession_clean}" > download_log.txt
    echo "Sample ID: ${prefix}" >> download_log.txt
    echo "Timestamp: \$(date)" >> download_log.txt

    if [ -f "${prefix}.fa" ]; then
        seq_count=\$(grep -c "^>" ${prefix}.fa || echo "0")
        seq_length=\$(grep -v "^>" ${prefix}.fa | tr -d '\\n' | wc -c || echo "0")
        echo "Sequences: \$seq_count" >> download_log.txt
        echo "Length: \$seq_length bp" >> download_log.txt
        echo "Status: SUCCESS" >> download_log.txt
    else
        echo "Status: FAILED" >> download_log.txt
    fi

    # Clean up temporary files
    rm -rf ncbi_dataset/ ${accession_clean}.zip ${prefix}_temp.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ncbi-datasets: \$(datasets --version 2>&1 | head -1 | cut -d' ' -f3 | cut -d'v' -f2)
        python: \$(python3 --version | cut -d' ' -f2)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.fa
    touch ${prefix}.json
    touch download_log.txt
    touch versions.yml
    """
}