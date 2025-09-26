process BWA_INDEX {
    tag "$meta.id"
    label 'process_single'

    container "quay.io/biocontainers/bwa:0.7.18--he4a0461_0"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.{amb,ann,bwt,pac,sa}")    , emit: index
    tuple val(meta), path("bwa_index_log.txt")         , emit: log
    path "versions.yml"                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo "🔧 Creating BWA index for ${meta.id}" | tee bwa_index_log.txt
    echo "Input: ${fasta}" >> bwa_index_log.txt
    echo "Started: \$(date)" >> bwa_index_log.txt
    echo "" >> bwa_index_log.txt

    # Create BWA index
    bwa index ${args} ${fasta} 2>&1 | tee -a bwa_index_log.txt

    # Check if index was created successfully
    if [ -f "${fasta}.bwt" ]; then
        echo "✅ BWA index created successfully" >> bwa_index_log.txt
        echo "Index files:" >> bwa_index_log.txt
        ls -lh ${fasta}.* | grep -E "\\.(amb|ann|bwt|pac|sa)\$" >> bwa_index_log.txt
    else
        echo "❌ BWA index creation failed" >> bwa_index_log.txt
        exit 1
    fi

    echo "Completed: \$(date)" >> bwa_index_log.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(bwa 2>&1 | grep -E '^Version' | cut -d' ' -f2)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${fasta}.amb
    touch ${fasta}.ann
    touch ${fasta}.bwt
    touch ${fasta}.pac
    touch ${fasta}.sa
    touch bwa_index_log.txt
    touch versions.yml
    """
}