/*
========================================================================================
    Create Minimap2 Index Module
========================================================================================
    Description: Generate minimap2 indices for efficient long-read and genome alignment
    Purpose: Enable efficient pangenome graph validation and sequence mapping
========================================================================================
*/

process CREATE_MINIMAP2_INDEX {
    tag "$meta.id"
    label 'process_medium'
    label 'process_long'

    container 'quay.io/biocontainers/minimap2:2.30--h577a1d6_0'

    publishDir "${params.outdir}/02_preprocessing/indices/minimap2", mode: params.publish_dir_mode,
        saveAs: { filename ->
            if (filename.startsWith("${meta.id}")) filename
            else null
        }

    input:
    tuple val(meta), path(genome)

    output:
    tuple val(meta), path("${meta.id}_minimap2_index")  , emit: index
    path "versions.yml"                                 , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''

    """
    # Create output directory
    mkdir -p ${prefix}_minimap2_index

    # Copy genome file to index directory with consistent naming
    cp ${genome} ${prefix}_minimap2_index/${prefix}.fa

    # Create minimap2 index for different alignment modes
    echo "Creating minimap2 indices for ${meta.id}..."

    # Genome-to-genome alignment index (asm5/asm10/asm20 presets)
    minimap2 -x asm5 -t ${task.cpus} -d ${prefix}_minimap2_index/${prefix}_asm5.mmi ${prefix}_minimap2_index/${prefix}.fa ${args}

    # Long-read alignment index (map-ont/map-pb presets)
    minimap2 -x map-ont -t ${task.cpus} -d ${prefix}_minimap2_index/${prefix}_ont.mmi ${prefix}_minimap2_index/${prefix}.fa ${args}

    # PacBio long-read alignment index
    minimap2 -x map-pb -t ${task.cpus} -d ${prefix}_minimap2_index/${prefix}_pb.mmi ${prefix}_minimap2_index/${prefix}.fa ${args}

    # Verify index creation
    if [ ! -f "${prefix}_minimap2_index/${prefix}_asm5.mmi" ]; then
        echo "❌ Minimap2 index creation failed for ${meta.id}"
        exit 1
    fi

    # Create index metadata
    cat > ${prefix}_minimap2_index/index_info.json << EOF
{
    "sample_id": "${meta.id}",
    "index_type": "minimap2",
    "reference_file": "${prefix}.fa",
    "creation_timestamp": "\$(date -Iseconds)",
    "minimap2_version": "\$(minimap2 --version)",
    "index_files": [
        "${prefix}_asm5.mmi",
        "${prefix}_ont.mmi",
        "${prefix}_pb.mmi"
    ],
    "preset_modes": {
        "asm5": "Accurate genome-to-genome alignment",
        "map-ont": "Oxford Nanopore long reads",
        "map-pb": "PacBio long reads"
    },
    "genome_size": \$(wc -c < ${prefix}_minimap2_index/${prefix}.fa),
    "ready_for_alignment": true
}
EOF

    # Print index file sizes and summary
    echo "✅ Minimap2 indices created successfully for ${meta.id}"
    echo "   Index directory: ${prefix}_minimap2_index/"
    echo "   Reference file: ${prefix}.fa"
    echo "   Index files created:"
    ls -lh ${prefix}_minimap2_index/*.mmi | while read line; do
        echo "     \$line"
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version)
    END_VERSIONS
    """
}