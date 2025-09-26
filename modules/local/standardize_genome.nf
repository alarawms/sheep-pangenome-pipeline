process STANDARDIZE_GENOME {
    tag "$meta.id"
    label 'process_medium'

    container "python:3.9-slim"

    input:
    tuple val(meta), path(fasta), path(metadata)

    output:
    tuple val(meta), path("*_standardized.fa")         , emit: fasta
    tuple val(meta), path("*_mapping.tsv")             , emit: mapping
    tuple val(meta), path("*_standardization_log.txt") , emit: log
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    python3 << 'EOF'
import re
import json
import datetime
from collections import defaultdict

def standardize_chromosome_name(seq_id):
    \"\"\"Standardize chromosome names for sheep genomes\"\"\"
    seq_id_clean = seq_id.strip()

    # Standard autosomes (1-26)
    match = re.match(r'^(chr)?([1-9]|1[0-9]|2[0-6])(\\.1)?\$', seq_id_clean, re.IGNORECASE)
    if match:
        return f"chr{match.group(2)}"

    match = re.match(r'^chromosome_?([1-9]|1[0-9]|2[0-6])', seq_id_clean, re.IGNORECASE)
    if match:
        return f"chr{match.group(1)}"

    # RefSeq autosome patterns for sheep (NC_019458.2 to NC_019483.2 = chr1-26)
    match = re.match(r'^NC_01945([8-9])|NC_0194[6-7][0-9]|NC_019483\\.[0-9]+', seq_id_clean)
    if match:
        if 'NC_019458' in seq_id_clean: return 'chr1'
        elif 'NC_019459' in seq_id_clean: return 'chr2'
        elif 'NC_019460' in seq_id_clean: return 'chr3'
        elif 'NC_019461' in seq_id_clean: return 'chr4'
        elif 'NC_019462' in seq_id_clean: return 'chr5'
        elif 'NC_019463' in seq_id_clean: return 'chr6'
        elif 'NC_019464' in seq_id_clean: return 'chr7'
        elif 'NC_019465' in seq_id_clean: return 'chr8'
        elif 'NC_019466' in seq_id_clean: return 'chr9'
        elif 'NC_019467' in seq_id_clean: return 'chr10'
        elif 'NC_019468' in seq_id_clean: return 'chr11'
        elif 'NC_019469' in seq_id_clean: return 'chr12'
        elif 'NC_019470' in seq_id_clean: return 'chr13'
        elif 'NC_019471' in seq_id_clean: return 'chr14'
        elif 'NC_019472' in seq_id_clean: return 'chr15'
        elif 'NC_019473' in seq_id_clean: return 'chr16'
        elif 'NC_019474' in seq_id_clean: return 'chr17'
        elif 'NC_019475' in seq_id_clean: return 'chr18'
        elif 'NC_019476' in seq_id_clean: return 'chr19'
        elif 'NC_019477' in seq_id_clean: return 'chr20'
        elif 'NC_019478' in seq_id_clean: return 'chr21'
        elif 'NC_019479' in seq_id_clean: return 'chr22'
        elif 'NC_019480' in seq_id_clean: return 'chr23'
        elif 'NC_019481' in seq_id_clean: return 'chr24'
        elif 'NC_019482' in seq_id_clean: return 'chr25'
        elif 'NC_019483' in seq_id_clean: return 'chr26'

    # X chromosome
    if re.match(r'^(chr)?(x|X)(\\.1)?\$', seq_id_clean):
        return 'chrX'
    if re.match(r'^chromosome_?[xX]', seq_id_clean, re.IGNORECASE):
        return 'chrX'
    if re.match(r'^NC_019484\\.[0-9]+', seq_id_clean):  # Sheep X chromosome
        return 'chrX'

    # Mitochondrial
    if re.match(r'^(chr)?(mt|MT|m|M)(\\.1)?\$', seq_id_clean):
        return 'chrMT'
    if re.match(r'^mitochondrion', seq_id_clean, re.IGNORECASE):
        return 'chrMT'
    if re.match(r'^NC_001941\\.[0-9]+', seq_id_clean):  # Sheep mitochondrial
        return 'chrMT'

    # Unplaced/unlocalized scaffolds
    if re.match(r'^(chr)?un|unplaced|scaffold', seq_id_clean, re.IGNORECASE):
        return f"scaffold_{seq_id_clean}"
    if re.match(r'^SUPER_', seq_id_clean):
        return f"scaffold_{seq_id_clean}"
    if re.match(r'^AMGL[0-9]+', seq_id_clean):
        return f"scaffold_{seq_id_clean}"

    # Default: keep as scaffold
    return f"scaffold_{seq_id_clean}"

# Load metadata
try:
    with open("${metadata}", 'r') as f:
        meta_data = json.load(f)
except:
    meta_data = {}

# Process FASTA file
input_file = "${fasta}"
output_file = "${prefix}_standardized.fa"
mapping_file = "${prefix}_mapping.tsv"
log_file = "${prefix}_standardization_log.txt"

mapping = []
chromosome_counts = defaultdict(int)
total_sequences = 0
standardized_sequences = 0

print(f"🔧 Standardizing chromosome names for ${meta.id}")

with open(input_file, 'r') as infile, \\
     open(output_file, 'w') as outfile, \\
     open(mapping_file, 'w') as mapfile, \\
     open(log_file, 'w') as logfile:

    # Write mapping header
    mapfile.write("original_id\\tstandardized_id\\tsequence_length\\tstatus\\n")

    # Write log header
    logfile.write(f"Standardization log for ${meta.id}\\n")
    logfile.write(f"Started: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\\n\\n")

    current_seq = []
    current_header = None

    for line in infile:
        line = line.strip()
        if line.startswith('>'):
            # Process previous sequence
            if current_header:
                seq_len = len(''.join(current_seq))
                seq_id = current_header[1:].split()[0]
                std_name = standardize_chromosome_name(seq_id)

                # Write standardized sequence
                outfile.write(f">{std_name}\\n")
                seq_str = ''.join(current_seq)
                for i in range(0, len(seq_str), 80):
                    outfile.write(seq_str[i:i+80] + '\\n')

                # Record mapping
                status = "STANDARD" if not std_name.startswith("scaffold_") else "SCAFFOLD"
                mapfile.write(f"{seq_id}\\t{std_name}\\t{seq_len}\\t{status}\\n")

                logfile.write(f"{seq_id} -> {std_name} ({seq_len} bp, {status})\\n")

                chromosome_counts[std_name] += 1
                if status == "STANDARD":
                    standardized_sequences += 1

                total_sequences += 1

            current_header = line
            current_seq = []
        else:
            current_seq.append(line)

    # Process last sequence
    if current_header:
        seq_len = len(''.join(current_seq))
        seq_id = current_header[1:].split()[0]
        std_name = standardize_chromosome_name(seq_id)

        outfile.write(f">{std_name}\\n")
        seq_str = ''.join(current_seq)
        for i in range(0, len(seq_str), 80):
            outfile.write(seq_str[i:i+80] + '\\n')

        status = "STANDARD" if not std_name.startswith("scaffold_") else "SCAFFOLD"
        mapfile.write(f"{seq_id}\\t{std_name}\\t{seq_len}\\t{status}\\n")

        logfile.write(f"{seq_id} -> {std_name} ({seq_len} bp, {status})\\n")

        chromosome_counts[std_name] += 1
        if status == "STANDARD":
            standardized_sequences += 1

        total_sequences += 1

    # Write summary
    logfile.write(f"\\nStandardization Summary:\\n")
    logfile.write(f"  Total sequences: {total_sequences}\\n")
    logfile.write(f"  Standard chromosomes: {standardized_sequences}\\n")
    logfile.write(f"  Scaffolds: {total_sequences - standardized_sequences}\\n")
    logfile.write(f"\\nChromosome counts:\\n")

    # Count standard chromosomes
    standard_chrs = [f"chr{i}" for i in range(1, 27)] + ['chrX', 'chrMT']
    found_chrs = []

    for chr_name in standard_chrs:
        if chr_name in chromosome_counts:
            found_chrs.append(chr_name)
            logfile.write(f"  {chr_name}: {chromosome_counts[chr_name]}\\n")

    missing_chrs = set(standard_chrs) - set(found_chrs)
    if missing_chrs:
        logfile.write(f"\\nMissing standard chromosomes: {', '.join(sorted(missing_chrs))}\\n")

    logfile.write(f"\\nCompleted: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\\n")

print(f"✅ Standardization completed for ${meta.id}")
print(f"   Standard chromosomes: {standardized_sequences}/{total_sequences}")
print(f"   Output: {output_file}")

EOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_standardized.fa
    touch ${prefix}_mapping.tsv
    touch ${prefix}_standardization_log.txt
    touch versions.yml
    """
}