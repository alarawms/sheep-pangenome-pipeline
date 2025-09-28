/*
========================================================================================
    Create BWA Index Module
========================================================================================
    Description: Generate BWA-MEM2 indices for efficient short-read alignment
    Purpose: Enable alignment validation and quality control for pangenome construction
========================================================================================
*/

process CREATE_BWA_INDEX {
    tag "$meta.id"
    label 'process_high'
    label 'process_long'

    container 'quay.io/biocontainers/bwa-mem2:2.2.1--he513fc3_0'

    publishDir "${params.outdir}/02_preprocessing/indices/bwa", mode: params.publish_dir_mode,
        saveAs: { filename ->
            if (filename.startsWith("${meta.id}")) filename
            else null
        }

    input:
    tuple val(meta), path(genome)

    output:
    tuple val(meta), path("${meta.id}_bwa_index")   , emit: index
    path "versions.yml"                             , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''

    """
    # Create output directory
    mkdir -p ${prefix}_bwa_index

    # Copy genome file to index directory with consistent naming
    cp ${genome} ${prefix}_bwa_index/${prefix}.fa

    # Create BWA-MEM2 index
    echo "Creating BWA-MEM2 index for ${meta.id}..."
    bwa-mem2 index \\
        -p ${prefix}_bwa_index/${prefix} \\
        ${prefix}_bwa_index/${prefix}.fa \\
        ${args}

    # Verify index creation
    if [ ! -f "${prefix}_bwa_index/${prefix}.bwt.2bit.64" ]; then
        echo "❌ BWA index creation failed for ${meta.id}"
        exit 1
    fi

    # Create index metadata
    cat > ${prefix}_bwa_index/index_info.json << EOF
{
    "sample_id": "${meta.id}",
    "index_type": "bwa-mem2",
    "reference_file": "${prefix}.fa",
    "creation_timestamp": "\$(date -Iseconds)",
    "bwa_version": "\$(bwa-mem2 version | head -1 | cut -d' ' -f2)",
    "index_files": [
        "${prefix}.0123",
        "${prefix}.amb",
        "${prefix}.ann",
        "${prefix}.bwt.2bit.64",
        "${prefix}.pac"
    ],
    "genome_size": \$(wc -c < ${prefix}_bwa_index/${prefix}.fa),
    "ready_for_alignment": true
}
EOF

    # Print summary
    echo "✅ BWA-MEM2 index created successfully for ${meta.id}"
    echo "   Index directory: ${prefix}_bwa_index/"
    echo "   Reference file: ${prefix}.fa"
    echo "   Index files: \$(ls -la ${prefix}_bwa_index/${prefix}.* | wc -l) files"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa-mem2: \$(bwa-mem2 version | head -1 | cut -d' ' -f2)
    END_VERSIONS
    """
}