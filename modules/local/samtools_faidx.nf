process SAMTOOLS_FAIDX {
    tag "$meta.id"
    label 'process_single'

    container "quay.io/biocontainers/samtools:1.19.2--h50ea8bc_1"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.fai")                    , emit: index
    tuple val(meta), path("samtools_index_log.txt")  , emit: log
    path "versions.yml"                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo "🔧 Creating samtools faidx index for ${meta.id}" | tee samtools_index_log.txt
    echo "Input: ${fasta}" >> samtools_index_log.txt
    echo "Started: \$(date)" >> samtools_index_log.txt
    echo "" >> samtools_index_log.txt

    # Create samtools faidx index
    samtools faidx ${args} ${fasta} 2>&1 | tee -a samtools_index_log.txt

    # Check if index was created successfully
    if [ -f "${fasta}.fai" ]; then
        echo "✅ Samtools faidx index created successfully" >> samtools_index_log.txt
        echo "Index file: ${fasta}.fai" >> samtools_index_log.txt

        # Show index statistics
        echo "" >> samtools_index_log.txt
        echo "Index contents (first 10 sequences):" >> samtools_index_log.txt
        head -10 ${fasta}.fai >> samtools_index_log.txt

        echo "" >> samtools_index_log.txt
        echo "Index summary:" >> samtools_index_log.txt
        total_seqs=\$(wc -l < ${fasta}.fai)
        total_length=\$(awk '{sum+=\$2} END {print sum}' ${fasta}.fai)
        echo "  Total sequences: \$total_seqs" >> samtools_index_log.txt
        echo "  Total length: \$total_length bp" >> samtools_index_log.txt
        echo "  Index file size: \$(stat -c%s ${fasta}.fai) bytes" >> samtools_index_log.txt

        # Check for standard chromosomes
        echo "" >> samtools_index_log.txt
        echo "Chromosome detection:" >> samtools_index_log.txt
        for chr in {1..26} X MT; do
            if grep -q "^chr\$chr" ${fasta}.fai; then
                length=\$(grep "^chr\$chr" ${fasta}.fai | cut -f2)
                echo "  chr\$chr: \$length bp ✅" >> samtools_index_log.txt
            else
                echo "  chr\$chr: Not found ❌" >> samtools_index_log.txt
            fi
        done

    else
        echo "❌ Samtools faidx index creation failed" >> samtools_index_log.txt
        exit 1
    fi

    echo "Completed: \$(date)" >> samtools_index_log.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -1 | cut -d' ' -f2)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${fasta}.fai
    touch samtools_index_log.txt
    touch versions.yml
    """
}