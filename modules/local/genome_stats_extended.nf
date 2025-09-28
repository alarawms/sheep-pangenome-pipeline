/*
========================================================================================
    Extended Genome Statistics Module
========================================================================================
    Description: Comprehensive genome statistics beyond basic validation
    Features: N50/L50, GC distribution, repeat content estimation, assembly metrics
========================================================================================
*/

process GENOME_STATS_EXTENDED {
    tag "$meta.id"
    label 'process_medium'

    container 'quay.io/biocontainers/seqkit:2.3.1--h9ee0642_0'

    publishDir "${params.outdir}/02_preprocessing/statistics", mode: params.publish_dir_mode,
        saveAs: { filename ->
            if (filename.endsWith('.json')) "${meta.id}_extended_stats.json"
            else if (filename.endsWith('.txt')) "${meta.id}_extended_stats.txt"
            else null
        }

    input:
    tuple val(meta), path(genome)

    output:
    tuple val(meta), path("${meta.id}_extended_stats.json")    , emit: stats
    tuple val(meta), path("${meta.id}_extended_stats.txt")     , emit: stats_txt
    path "versions.yml"                                        , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    #!/usr/bin/env python3

    import json
    import sys
    import subprocess
    from pathlib import Path
    from collections import defaultdict, Counter
    import re

    def run_command(cmd):
        \"\"\"Run shell command and return output\"\"\"
        try:
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True)
            return result.stdout.strip()
        except subprocess.CalledProcessError as e:
            print(f"Command failed: {cmd}")
            print(f"Error: {e.stderr}")
            return None

    def parse_fasta_sequences(fasta_file):
        \"\"\"Parse FASTA file and extract detailed statistics\"\"\"

        sequences = []
        current_seq = ""
        current_header = ""

        with open(fasta_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if current_seq:
                        sequences.append({
                            'header': current_header,
                            'sequence': current_seq,
                            'length': len(current_seq)
                        })
                    current_header = line[1:]
                    current_seq = ""
                else:
                    current_seq += line.upper()

            # Add last sequence
            if current_seq:
                sequences.append({
                    'header': current_header,
                    'sequence': current_seq,
                    'length': len(current_seq)
                })

        return sequences

    def calculate_nx_lx_stats(lengths, x_values=[10, 20, 25, 30, 40, 50, 60, 70, 75, 80, 90]):
        \"\"\"Calculate NX and LX statistics for multiple percentiles\"\"\"

        sorted_lengths = sorted(lengths, reverse=True)
        total_length = sum(lengths)
        nx_lx_stats = {}

        for x in x_values:
            target_length = total_length * (x / 100)
            cumulative_length = 0
            lx_count = 0

            for length in sorted_lengths:
                cumulative_length += length
                lx_count += 1
                if cumulative_length >= target_length:
                    nx_lx_stats[f'N{x}'] = length
                    nx_lx_stats[f'L{x}'] = lx_count
                    break

        return nx_lx_stats

    def analyze_gc_distribution(sequences):
        \"\"\"Analyze GC content distribution across sequences\"\"\"

        gc_contents = []
        gc_by_length_bins = defaultdict(list)

        for seq_data in sequences:
            sequence = seq_data['sequence']
            length = seq_data['length']

            if length > 0:
                gc_count = sequence.count('G') + sequence.count('C')
                gc_content = (gc_count / length) * 100
                gc_contents.append(gc_content)

                # Bin by sequence length
                if length >= 1e6:
                    gc_by_length_bins['large_scaffolds'].append(gc_content)
                elif length >= 100e3:
                    gc_by_length_bins['medium_scaffolds'].append(gc_content)
                elif length >= 10e3:
                    gc_by_length_bins['small_scaffolds'].append(gc_content)
                else:
                    gc_by_length_bins['contigs'].append(gc_content)

        # Calculate statistics
        gc_stats = {}
        if gc_contents:
            gc_contents_sorted = sorted(gc_contents)
            n = len(gc_contents_sorted)

            gc_stats = {
                'mean': sum(gc_contents) / len(gc_contents),
                'median': gc_contents_sorted[n//2] if n % 2 == 1 else (gc_contents_sorted[n//2-1] + gc_contents_sorted[n//2]) / 2,
                'min': min(gc_contents),
                'max': max(gc_contents),
                'std': (sum((x - gc_stats.get('mean', 0))**2 for x in gc_contents) / len(gc_contents))**0.5 if len(gc_contents) > 1 else 0
            }

        # GC by length bins
        gc_by_bins = {}
        for bin_name, gc_values in gc_by_length_bins.items():
            if gc_values:
                gc_by_bins[bin_name] = {
                    'count': len(gc_values),
                    'mean_gc': sum(gc_values) / len(gc_values),
                    'min_gc': min(gc_values),
                    'max_gc': max(gc_values)
                }

        return gc_stats, gc_by_bins

    def estimate_repeat_content(sequences, sample_size=10):
        \"\"\"Estimate repeat content using k-mer analysis\"\"\"

        # Sample largest sequences for efficiency
        large_sequences = sorted([s for s in sequences if s['length'] >= 1000],
                                key=lambda x: x['length'], reverse=True)[:sample_size]

        if not large_sequences:
            return {'estimated_repeat_percentage': 0, 'analysis_note': 'No sequences >= 1kb for repeat analysis'}

        total_analyzed_length = 0
        repeat_regions = 0
        kmer_size = 21

        for seq_data in large_sequences:
            sequence = seq_data['sequence']
            seq_length = len(sequence)

            if seq_length < kmer_size:
                continue

            total_analyzed_length += seq_length

            # Count k-mer frequencies
            kmer_counts = Counter()
            for i in range(seq_length - kmer_size + 1):
                kmer = sequence[i:i + kmer_size]
                if 'N' not in kmer:  # Skip k-mers with gaps
                    kmer_counts[kmer] += 1

            # Estimate repeats (k-mers appearing more than once)
            repetitive_kmers = sum(count - 1 for count in kmer_counts.values() if count > 1)
            repeat_regions += repetitive_kmers

        if total_analyzed_length > 0:
            repeat_percentage = (repeat_regions / total_analyzed_length) * 100
        else:
            repeat_percentage = 0

        return {
            'estimated_repeat_percentage': repeat_percentage,
            'analyzed_sequences': len(large_sequences),
            'analyzed_length': total_analyzed_length,
            'kmer_size': kmer_size,
            'analysis_note': f'Estimated from largest {len(large_sequences)} sequences using {kmer_size}-mer analysis'
        }

    def analyze_sequence_composition(sequences):
        \"\"\"Analyze nucleotide composition and quality\"\"\"

        total_bases = {'A': 0, 'T': 0, 'G': 0, 'C': 0, 'N': 0, 'OTHER': 0}
        total_length = 0

        for seq_data in sequences:
            sequence = seq_data['sequence']
            total_length += len(sequence)

            for base in sequence:
                if base in total_bases:
                    total_bases[base] += 1
                else:
                    total_bases['OTHER'] += 1

        # Calculate percentages
        composition = {}
        if total_length > 0:
            for base, count in total_bases.items():
                composition[f'{base.lower()}_percentage'] = (count / total_length) * 100

        # Additional metrics
        composition['gc_percentage'] = composition.get('g_percentage', 0) + composition.get('c_percentage', 0)
        composition['at_percentage'] = composition.get('a_percentage', 0) + composition.get('t_percentage', 0)

        return composition

    def identify_chromosome_sequences(sequences):
        \"\"\"Identify likely chromosome-level sequences\"\"\"

        # Sheep chromosome size expectations (approximate, in bp)
        sheep_chr_sizes = {
            'chr1': (270e6, 290e6), 'chr2': (240e6, 260e6), 'chr3': (210e6, 240e6),
            'chr4': (110e6, 130e6), 'chr5': (100e6, 120e6), 'chr6': (110e6, 130e6),
            'chr7': (90e6, 110e6), 'chr8': (80e6, 100e6), 'chr9': (85e6, 105e6),
            'chr10': (80e6, 100e6), 'chr11': (55e6, 70e6), 'chr12': (70e6, 90e6),
            'chr13': (75e6, 95e6), 'chr14': (55e6, 75e6), 'chr15': (70e6, 90e6),
            'chr16': (60e6, 80e6), 'chr17': (60e6, 80e6), 'chr18': (55e6, 75e6),
            'chr19': (50e6, 70e6), 'chr20': (40e6, 60e6), 'chr21': (40e6, 60e6),
            'chr22': (40e6, 60e6), 'chr23': (40e6, 60e6), 'chr24': (35e6, 55e6),
            'chr25': (35e6, 55e6), 'chr26': (35e6, 55e6), 'chr27': (130e6, 150e6),  # X chromosome
            'chrMT': (15000, 20000)  # Mitochondrial
        }

        chromosome_candidates = []
        large_sequences = [s for s in sequences if s['length'] >= 30e6]  # Potential chromosomes

        for seq_data in large_sequences:
            header = seq_data['header'].lower()
            length = seq_data['length']

            # Check for chromosome indicators in header
            chr_match = None
            for chr_name in sheep_chr_sizes:
                if chr_name.lower() in header or chr_name.replace('chr', '') in header:
                    min_size, max_size = sheep_chr_sizes[chr_name]
                    if min_size <= length <= max_size:
                        chr_match = chr_name
                        break

            chromosome_candidates.append({
                'header': seq_data['header'],
                'length': length,
                'predicted_chromosome': chr_match,
                'size_rank': len([s for s in sequences if s['length'] > length]) + 1
            })

        return sorted(chromosome_candidates, key=lambda x: x['length'], reverse=True)

    def generate_extended_statistics(fasta_file, sample_id):
        \"\"\"Generate comprehensive genome statistics\"\"\"

        print(f"Generating extended statistics for {sample_id}")

        # Parse sequences
        sequences = parse_fasta_sequences(fasta_file)

        if not sequences:
            return {'error': 'No sequences found in FASTA file'}

        # Basic statistics
        lengths = [s['length'] for s in sequences]
        total_length = sum(lengths)
        num_sequences = len(sequences)

        basic_stats = {
            'sample_id': sample_id,
            'total_sequences': num_sequences,
            'total_length': total_length,
            'mean_sequence_length': total_length / num_sequences if num_sequences > 0 else 0,
            'median_sequence_length': sorted(lengths)[num_sequences // 2] if num_sequences > 0 else 0,
            'min_sequence_length': min(lengths) if lengths else 0,
            'max_sequence_length': max(lengths) if lengths else 0
        }

        # N50/L50 statistics
        nx_lx_stats = calculate_nx_lx_stats(lengths)

        # Sequence size categories
        size_categories = {
            'very_large': len([l for l in lengths if l >= 50e6]),    # >= 50 Mb (likely chromosomes)
            'large': len([l for l in lengths if 1e6 <= l < 50e6]),   # 1-50 Mb (large scaffolds)
            'medium': len([l for l in lengths if 10e3 <= l < 1e6]),  # 10kb-1Mb (small scaffolds)
            'small': len([l for l in lengths if 1e3 <= l < 10e3]),   # 1-10kb (contigs)
            'tiny': len([l for l in lengths if l < 1e3])             # < 1kb (fragments)
        }

        # Composition analysis
        composition = analyze_sequence_composition(sequences)

        # GC distribution analysis
        gc_stats, gc_by_bins = analyze_gc_distribution(sequences)

        # Repeat content estimation
        repeat_analysis = estimate_repeat_content(sequences)

        # Chromosome identification
        chromosome_candidates = identify_chromosome_sequences(sequences)

        # Assembly quality metrics
        assembly_metrics = {
            'scaffold_n50': nx_lx_stats.get('N50', 0),
            'scaffold_l50': nx_lx_stats.get('L50', 0),
            'contig_n50': nx_lx_stats.get('N50', 0),  # Approximation - would need gap info for true contig N50
            'assembly_span': total_length,
            'num_scaffolds': num_sequences,
            'largest_scaffold': max(lengths) if lengths else 0,
            'chromosome_level_scaffolds': size_categories['very_large'],
            'gaps_estimated': composition.get('n_percentage', 0)
        }

        # Quality assessment
        quality_assessment = {
            'chromosome_completeness': min(27, size_categories['very_large']) / 27 * 100,  # 27 sheep chromosomes
            'contiguity_score': min(100, nx_lx_stats.get('N50', 0) / 1e6 * 2),  # N50 in Mb * 2, capped at 100
            'size_accuracy': max(0, 100 - abs(total_length - 2.8e9) / 2.8e9 * 100),  # Expected sheep genome ~2.8Gb
            'gap_quality': max(0, 100 - composition.get('n_percentage', 0) * 20)  # Penalize high N content
        }

        # Compile final statistics
        extended_stats = {
            'analysis_timestamp': subprocess.run(['date', '-Iseconds'], capture_output=True, text=True).stdout.strip(),
            'sample_info': basic_stats,
            'assembly_metrics': assembly_metrics,
            'nx_lx_statistics': nx_lx_stats,
            'sequence_size_distribution': size_categories,
            'nucleotide_composition': composition,
            'gc_content_analysis': {
                'overall_statistics': gc_stats,
                'by_sequence_size': gc_by_bins
            },
            'repeat_content_estimation': repeat_analysis,
            'chromosome_candidates': chromosome_candidates[:10],  # Top 10 largest
            'quality_assessment': quality_assessment
        }

        return extended_stats

    # Main processing
    try:
        fasta_file = "${genome}"
        sample_id = "${meta.id}"
        json_output = "${meta.id}_extended_stats.json"
        txt_output = "${meta.id}_extended_stats.txt"

        # Generate statistics
        stats = generate_extended_statistics(fasta_file, sample_id)

        # Write JSON output
        with open(json_output, 'w') as f:
            json.dump(stats, f, indent=2)

        # Generate human-readable report
        with open(txt_output, 'w') as f:
            f.write(f"Extended Genome Statistics Report\\n")
            f.write(f"Sample: {sample_id}\\n")
            f.write(f"Analysis: {stats.get('analysis_timestamp', 'Unknown')}\\n")
            f.write(f"{'='*60}\\n\\n")

            # Basic information
            basic = stats.get('sample_info', {})
            f.write(f"BASIC STATISTICS\\n")
            f.write(f"  Total sequences: {basic.get('total_sequences', 0):,}\\n")
            f.write(f"  Total length: {basic.get('total_length', 0):,} bp ({basic.get('total_length', 0)/1e9:.2f} Gb)\\n")
            f.write(f"  Mean length: {basic.get('mean_sequence_length', 0):,.0f} bp\\n")
            f.write(f"  Largest sequence: {basic.get('max_sequence_length', 0):,} bp\\n\\n")

            # Assembly metrics
            assembly = stats.get('assembly_metrics', {})
            f.write(f"ASSEMBLY METRICS\\n")
            f.write(f"  Scaffold N50: {assembly.get('scaffold_n50', 0):,} bp ({assembly.get('scaffold_n50', 0)/1e6:.1f} Mb)\\n")
            f.write(f"  Scaffold L50: {assembly.get('scaffold_l50', 0)}\\n")
            f.write(f"  Chromosome-level scaffolds: {assembly.get('chromosome_level_scaffolds', 0)}\\n")
            f.write(f"  Estimated gaps: {assembly.get('gaps_estimated', 0):.2f}%\\n\\n")

            # Quality assessment
            quality = stats.get('quality_assessment', {})
            f.write(f"QUALITY ASSESSMENT\\n")
            f.write(f"  Chromosome completeness: {quality.get('chromosome_completeness', 0):.1f}%\\n")
            f.write(f"  Contiguity score: {quality.get('contiguity_score', 0):.1f}/100\\n")
            f.write(f"  Size accuracy: {quality.get('size_accuracy', 0):.1f}%\\n")
            f.write(f"  Gap quality: {quality.get('gap_quality', 0):.1f}/100\\n\\n")

            # Composition
            comp = stats.get('nucleotide_composition', {})
            f.write(f"NUCLEOTIDE COMPOSITION\\n")
            f.write(f"  GC content: {comp.get('gc_percentage', 0):.2f}%\\n")
            f.write(f"  AT content: {comp.get('at_percentage', 0):.2f}%\\n")
            f.write(f"  N content: {comp.get('n_percentage', 0):.2f}%\\n\\n")

            # Chromosome candidates
            chr_candidates = stats.get('chromosome_candidates', [])
            if chr_candidates:
                f.write(f"TOP CHROMOSOME CANDIDATES\\n")
                for i, candidate in enumerate(chr_candidates[:5], 1):
                    f.write(f"  {i}. {candidate.get('header', 'Unknown')[:50]}{'...' if len(candidate.get('header', '')) > 50 else ''}\\n")
                    f.write(f"     Length: {candidate.get('length', 0):,} bp ({candidate.get('length', 0)/1e6:.1f} Mb)\\n")
                    if candidate.get('predicted_chromosome'):
                        f.write(f"     Predicted: {candidate['predicted_chromosome']}\\n")
                    f.write(f"\\n")

        print(f"✅ Extended statistics generated for {sample_id}")
        print(f"   JSON: {json_output}")
        print(f"   Report: {txt_output}")

        # Print key metrics
        if 'assembly_metrics' in stats:
            metrics = stats['assembly_metrics']
            print(f"   N50: {metrics.get('scaffold_n50', 0)/1e6:.1f} Mb")
            print(f"   Chromosome-level scaffolds: {metrics.get('chromosome_level_scaffolds', 0)}")

    except Exception as e:
        print(f"❌ Error generating statistics for {sample_id}: {str(e)}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
        seqkit: \$(seqkit version | cut -d' ' -f2)
    END_VERSIONS
    """
}