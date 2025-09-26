process MINIMAP2_INDEX {
    tag "$meta.id"
    label 'process_medium'

    container "quay.io/biocontainers/minimap2:2.28--he4a0461_2"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.mmi")                     , emit: index
    tuple val(meta), path("minimap2_index_log.txt")   , emit: log
    path "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo "🔧 Creating minimap2 index for ${meta.id}" | tee minimap2_index_log.txt
    echo "Input: ${fasta}" >> minimap2_index_log.txt
    echo "Started: \$(date)" >> minimap2_index_log.txt
    echo "" >> minimap2_index_log.txt

    # Create minimap2 index
    minimap2 \\
        ${args} \\
        -t ${task.cpus} \\
        -d ${prefix}.mmi \\
        ${fasta} \\
        2>&1 | tee -a minimap2_index_log.txt

    # Check if index was created successfully
    if [ -f "${prefix}.mmi" ]; then
        echo "✅ Minimap2 index created successfully" >> minimap2_index_log.txt
        echo "Index file: ${prefix}.mmi" >> minimap2_index_log.txt
        ls -lh ${prefix}.mmi >> minimap2_index_log.txt

        # Get index statistics
        echo "" >> minimap2_index_log.txt
        echo "Index statistics:" >> minimap2_index_log.txt
        file_size=\$(stat -c%s "${prefix}.mmi")
        echo "  File size: \$(echo \$file_size | numfmt --to=iec)" >> minimap2_index_log.txt
    else
        echo "❌ Minimap2 index creation failed" >> minimap2_index_log.txt
        exit 1
    fi

    echo "Completed: \$(date)" >> minimap2_index_log.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.mmi
    touch minimap2_index_log.txt
    touch versions.yml
    """
}