/*
========================================================================================
    Rename Chromosomes Module
========================================================================================
    Description: Standardize sheep chromosome naming to conventional format (1-27, MT)
    Algorithm: Pattern matching and systematic renaming based on sheep karyotype
========================================================================================
*/

process RENAME_CHROMOSOMES {
    tag "$meta.id"
    label 'process_medium'
    label 'process_long'

    container 'biocontainers/biopython:1.79'

    publishDir "${params.outdir}/02_preprocessing/chromosome_mapping", mode: params.publish_dir_mode,
        saveAs: { filename ->
            if (filename.endsWith('.tsv')) "${meta.id}_chromosome_mapping.tsv"
            else null
        }

    input:
    tuple val(meta), path(genome)

    output:
    tuple val(meta), path("${meta.id}_renamed.fa")          , emit: renamed_genomes
    tuple val(meta), path("${meta.id}_chr_mapping.tsv")     , emit: chromosome_map
    path "versions.yml"                                     , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    #!/usr/bin/env python3

    import re
    import sys
    from pathlib import Path
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    def create_sheep_chromosome_mapping():
        \"\"\"Create comprehensive sheep chromosome naming mapping\"\"\"

        # Standard sheep chromosomes (1-26 autosomes + X chromosome = 27 + MT)
        sheep_chromosomes = {
            # Standard numeric chromosomes
            **{f'chr{i}': f'chr{i}' for i in range(1, 27)},  # chr1-chr26
            **{f'{i}': f'chr{i}' for i in range(1, 27)},     # 1-26 -> chr1-chr26

            # X chromosome variations
            'chrX': 'chr27',
            'X': 'chr27',
            'chr27': 'chr27',
            '27': 'chr27',

            # Mitochondrial variations
            'chrMT': 'chrMT',
            'chrM': 'chrMT',
            'MT': 'chrMT',
            'M': 'chrMT',
            'chrMito': 'chrMT',
            'chrMitochondrion': 'chrMT',
            'mitochondrion': 'chrMT',
            'mito': 'chrMT',

            # Assembly-specific naming patterns (common in sheep assemblies)
            **{f'chromosome_{i}': f'chr{i}' for i in range(1, 27)},
            **{f'Chromosome_{i}': f'chr{i}' for i in range(1, 27)},
            **{f'CHROMOSOME_{i}': f'chr{i}' for i in range(1, 27)},

            # Scaffold/Contig patterns for major chromosomes
            **{f'scaffold_{i}': f'chr{i}' for i in range(1, 27)},
            **{f'Scaffold_{i}': f'chr{i}' for i in range(1, 27)},

            # NCBI/RefSeq patterns
            **{f'NC_0{40000+i-1:05d}.1': f'chr{i}' for i in range(1, 27)},
            **{f'NC_0{40000+i-1:05d}': f'chr{i}' for i in range(1, 27)},
        }

        # Add common assembly accession patterns
        sheep_chromosomes.update({
            'NC_040252.1': 'chr1', 'NC_040253.1': 'chr2', 'NC_040254.1': 'chr3',
            'NC_040255.1': 'chr4', 'NC_040256.1': 'chr5', 'NC_040257.1': 'chr6',
            'NC_040258.1': 'chr7', 'NC_040259.1': 'chr8', 'NC_040260.1': 'chr9',
            'NC_040261.1': 'chr10', 'NC_040262.1': 'chr11', 'NC_040263.1': 'chr12',
            'NC_040264.1': 'chr13', 'NC_040265.1': 'chr14', 'NC_040266.1': 'chr15',
            'NC_040267.1': 'chr16', 'NC_040268.1': 'chr17', 'NC_040269.1': 'chr18',
            'NC_040270.1': 'chr19', 'NC_040271.1': 'chr20', 'NC_040272.1': 'chr21',
            'NC_040273.1': 'chr22', 'NC_040274.1': 'chr23', 'NC_040275.1': 'chr24',
            'NC_040276.1': 'chr25', 'NC_040277.1': 'chr26', 'NC_040278.1': 'chr27',
            'NC_001941.1': 'chrMT'
        })

        return sheep_chromosomes

    def intelligent_chromosome_detection(seq_name, seq_length, sheep_mapping):
        \"\"\"Intelligent chromosome detection using multiple heuristics\"\"\"

        # Direct mapping check
        clean_name = seq_name.split()[0]  # Take first part before space
        if clean_name in sheep_mapping:
            return sheep_mapping[clean_name], 'direct_mapping', 1.0

        # Pattern-based detection
        patterns = [
            # Standard chromosome patterns
            (r'^chr(\d+)$', lambda m: f'chr{int(m.group(1))}' if 1 <= int(m.group(1)) <= 27 else None),
            (r'^(\d+)$', lambda m: f'chr{int(m.group(1))}' if 1 <= int(m.group(1)) <= 27 else None),
            (r'chromosome[_\\s]+(\\d+)', lambda m: f'chr{int(m.group(1))}' if 1 <= int(m.group(1)) <= 27 else None),

            # X chromosome patterns
            (r'^(chr)?X$', lambda m: 'chr27'),
            (r'chromosome[_\\s]+X', lambda m: 'chr27'),

            # Mitochondrial patterns
            (r'^(chr)?(MT|M|mito)', lambda m: 'chrMT'),
            (r'mitochondr', lambda m: 'chrMT'),

            # Assembly-specific patterns
            (r'scaffold[_\\s]+(\\d+)', lambda m: f'chr{int(m.group(1))}' if 1 <= int(m.group(1)) <= 27 else None),
            (r'contig[_\\s]+(\\d+)', lambda m: f'chr{int(m.group(1))}' if 1 <= int(m.group(1)) <= 27 else None),
        ]

        for pattern, transformer in patterns:
            match = re.search(pattern, seq_name, re.IGNORECASE)
            if match:
                result = transformer(match)
                if result:
                    return result, 'pattern_matching', 0.9

        # Size-based heuristics for major chromosomes
        # Sheep chromosome sizes (approximate, in bp)
        sheep_chr_sizes = {
            1: (275e6, 285e6), 2: (248e6, 252e6), 3: (222e6, 228e6),
            4: (117e6, 123e6), 5: (107e6, 113e6), 6: (115e6, 121e6),
            7: (98e6, 104e6),  8: (89e6, 95e6),   9: (93e6, 99e6),
            10: (86e6, 92e6),  11: (59e6, 65e6),  12: (76e6, 82e6),
            13: (81e6, 87e6),  14: (60e6, 66e6),  15: (74e6, 80e6),
            16: (68e6, 74e6),  17: (66e6, 72e6),  18: (62e6, 68e6),
            19: (58e6, 64e6),  20: (48e6, 54e6),  21: (46e6, 52e6),
            22: (48e6, 54e6),  23: (45e6, 51e6),  24: (40e6, 46e6),
            25: (42e6, 48e6),  26: (43e6, 49e6),  27: (135e6, 145e6),  # X chromosome
        }

        if seq_length > 30e6:  # Only for substantial sequences
            for chr_num, (min_size, max_size) in sheep_chr_sizes.items():
                if min_size <= seq_length <= max_size:
                    return f'chr{chr_num}', 'size_heuristic', 0.7

        # Mitochondrial size check
        if 15000 <= seq_length <= 20000:
            if any(term in seq_name.lower() for term in ['mt', 'mito', 'mitochondr']):
                return 'chrMT', 'size_mito_heuristic', 0.8

        return None, 'unassigned', 0.0

    def process_genome_file(input_file, output_file, mapping_file):
        \"\"\"Process genome file and rename chromosomes\"\"\"

        print(f"Processing genome file: {input_file}")

        sheep_mapping = create_sheep_chromosome_mapping()
        renamed_sequences = []
        mapping_records = []

        # Statistics
        total_sequences = 0
        renamed_sequences_count = 0
        unassigned_sequences = 0

        # Process sequences
        for record in SeqIO.parse(input_file, "fasta"):
            total_sequences += 1
            original_name = record.id
            original_description = record.description
            sequence_length = len(record.seq)

            # Attempt intelligent chromosome detection
            new_name, method, confidence = intelligent_chromosome_detection(
                original_name, sequence_length, sheep_mapping
            )

            if new_name:
                # Create renamed record
                new_record = SeqRecord(
                    record.seq,
                    id=new_name,
                    description=f"{new_name} [original: {original_name}] [method: {method}] [confidence: {confidence:.2f}]"
                )
                renamed_sequences.append(new_record)
                renamed_sequences_count += 1

                # Record mapping
                mapping_records.append({
                    'original_name': original_name,
                    'new_name': new_name,
                    'sequence_length': sequence_length,
                    'method': method,
                    'confidence': confidence,
                    'original_description': original_description
                })

                print(f"Renamed: {original_name} -> {new_name} (method: {method}, confidence: {confidence:.2f})")

            else:
                # Keep original name for unassigned sequences
                renamed_sequences.append(record)
                unassigned_sequences += 1

                mapping_records.append({
                    'original_name': original_name,
                    'new_name': original_name,
                    'sequence_length': sequence_length,
                    'method': 'unchanged',
                    'confidence': 1.0,
                    'original_description': original_description
                })

        # Write renamed genome
        with open(output_file, 'w') as out_handle:
            SeqIO.write(renamed_sequences, out_handle, "fasta")

        # Write mapping file
        with open(mapping_file, 'w') as map_handle:
            map_handle.write("original_name\\tnew_name\\tsequence_length\\tmethod\\tconfidence\\toriginal_description\\n")
            for record in mapping_records:
                map_handle.write(f"{record['original_name']}\\t{record['new_name']}\\t{record['sequence_length']}\\t{record['method']}\\t{record['confidence']:.2f}\\t{record['original_description']}\\n")

        # Print summary
        print(f"")
        print(f"Chromosome renaming summary for {meta.id}:")
        print(f"  Total sequences: {total_sequences}")
        print(f"  Renamed sequences: {renamed_sequences_count}")
        print(f"  Unassigned sequences: {unassigned_sequences}")
        print(f"  Success rate: {(renamed_sequences_count/total_sequences*100):.1f}%")

        # Count sheep chromosomes
        sheep_chrs = set()
        for record in mapping_records:
            if record['new_name'].startswith('chr') and record['new_name'] != record['original_name']:
                sheep_chrs.add(record['new_name'])

        print(f"  Identified sheep chromosomes: {len(sheep_chrs)}")
        print(f"  Chromosomes found: {sorted(sheep_chrs)}")

        return {
            'total_sequences': total_sequences,
            'renamed_sequences': renamed_sequences_count,
            'unassigned_sequences': unassigned_sequences,
            'success_rate': renamed_sequences_count/total_sequences,
            'sheep_chromosomes_found': len(sheep_chrs),
            'chromosomes_identified': sorted(list(sheep_chrs))
        }

    # Main processing
    try:
        input_file = "${genome}"
        output_file = "${meta.id}_renamed.fa"
        mapping_file = "${meta.id}_chr_mapping.tsv"

        result = process_genome_file(input_file, output_file, mapping_file)

        print(f"")
        print(f"✅ Successfully processed {meta.id}")
        print(f"   Output: {output_file}")
        print(f"   Mapping: {mapping_file}")

    except Exception as e:
        print(f"❌ Error processing {meta.id}: {str(e)}")
        sys.exit(1)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        biopython: \$(python -c "import Bio; print(Bio.__version__)")
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """
}