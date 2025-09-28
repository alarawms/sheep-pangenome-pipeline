process BWA_INDEX {
    tag "$meta.id"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/39/39e36dcb9a2e94b69a5f41bfebfd7c8c995ecb503afe65a3af2088c9e0fd4593/data' :
        'community.wave.seqera.io/library/bwa-mem2_htslib_samtools:db98f81f55b64113' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.{0123,amb,ann,bwt.2bit.64,pac}")  , emit: index
    tuple val(meta), path("bwa_mem2_index_log.txt")            , emit: log
    path "versions.yml"                                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo "🔧 Creating BWA-MEM2 index for ${meta.id}" | tee bwa_mem2_index_log.txt
    echo "Input: ${fasta}" >> bwa_mem2_index_log.txt
    echo "Using ${task.cpus} threads for multi-threaded indexing" >> bwa_mem2_index_log.txt
    echo "Memory available: ${task.memory}" >> bwa_mem2_index_log.txt
    echo "Started: \$(date)" >> bwa_mem2_index_log.txt
    echo "" >> bwa_mem2_index_log.txt

    # Create BWA-MEM2 index with multi-threading support
    echo "Creating BWA-MEM2 index with optimized parameters" >> bwa_mem2_index_log.txt
    bwa-mem2 index \\
        ${args} \\
        ${fasta} \\
        2>&1 | tee -a bwa_mem2_index_log.txt

    # Check if index was created successfully
    if [ -f "${fasta}.bwt.2bit.64" ]; then
        echo "✅ BWA-MEM2 index created successfully" >> bwa_mem2_index_log.txt
        echo "Index files:" >> bwa_mem2_index_log.txt
        ls -lh ${fasta}.* | grep -E "\\.(0123|amb|ann|bwt\\.2bit\\.64|pac)\$" >> bwa_mem2_index_log.txt

        # Get index statistics
        echo "" >> bwa_mem2_index_log.txt
        echo "Index statistics:" >> bwa_mem2_index_log.txt
        total_size=\$(ls -la ${fasta}.* | grep -E "\\.(0123|amb|ann|bwt\\.2bit\\.64|pac)\$" | awk '{sum += \$5} END {print sum}')
        echo "  Total index size: \$(echo \$total_size | numfmt --to=iec)" >> bwa_mem2_index_log.txt
    else
        echo "❌ BWA-MEM2 index creation failed" >> bwa_mem2_index_log.txt
        exit 1
    fi

    echo "Completed: \$(date)" >> bwa_mem2_index_log.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa-mem2: \$(bwa-mem2 version 2>&1 | head -1 | sed 's/.*bwa-mem2 //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${fasta}.0123
    touch ${fasta}.amb
    touch ${fasta}.ann
    touch ${fasta}.bwt.2bit.64
    touch ${fasta}.pac
    touch bwa_mem2_index_log.txt
    touch versions.yml
    """
}