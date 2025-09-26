process GENOME_QC_EXTENDED {
    tag "$meta.id"
    label 'process_medium'

    container "quay.io/biocontainers/python:3.9--1"

    input:
    tuple val(meta), path(fasta), path(standardization_log)

    output:
    tuple val(meta), path("*_extended_stats.json")    , emit: stats
    tuple val(meta), path("*_quality_report.html")    , emit: report
    tuple val(meta), path("*_qc_log.txt")            , emit: log
    path "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    # Extended Quality Control Analysis for Sheep Genomes
    python3 << 'EOF'
import json
import re
import math
import datetime
from collections import Counter, defaultdict
from pathlib import Path

def calculate_n50(lengths):
    \"\"\"Calculate N50 from list of sequence lengths\"\"\"
    if not lengths:
        return 0

    sorted_lengths = sorted(lengths, reverse=True)
    total_length = sum(sorted_lengths)
    target = total_length / 2

    cumulative = 0
    for length in sorted_lengths:
        cumulative += length
        if cumulative >= target:
            return length
    return sorted_lengths[-1]

def analyze_sequence_composition(seq):
    \"\"\"Analyze nucleotide composition and patterns\"\"\"
    seq = seq.upper()
    counts = Counter(seq)
    total = len(seq)

    if total == 0:
        return {}

    # Basic composition
    composition = {
        'A': counts.get('A', 0),
        'T': counts.get('T', 0),
        'G': counts.get('G', 0),
        'C': counts.get('C', 0),
        'N': counts.get('N', 0),
        'other': total - sum(counts.get(base, 0) for base in 'ATGCN')
    }

    # Percentages
    percentages = {f"{base}_percent": (count/total)*100 for base, count in composition.items()}

    # GC content
    gc_count = composition['G'] + composition['C']
    at_count = composition['A'] + composition['T']
    gc_percent = (gc_count / total) * 100 if total > 0 else 0
    at_percent = (at_count / total) * 100 if total > 0 else 0

    # N content (gaps)
    n_percent = (composition['N'] / total) * 100 if total > 0 else 0

    return {
        **composition,
        **percentages,
        'total_length': total,
        'gc_content': gc_percent,
        'at_content': at_percent,
        'n_content': n_percent,
        'gc_count': gc_count,
        'at_count': at_count
    }

def detect_repeats(seq, min_repeat_length=3, max_check_length=10000):
    \"\"\"Detect simple tandem repeats\"\"\"
    seq = seq.upper()[:max_check_length]  # Limit for performance
    repeats = defaultdict(int)

    for repeat_len in range(min_repeat_length, min(20, len(seq)//4)):
        for i in range(len(seq) - repeat_len * 2):
            motif = seq[i:i+repeat_len]
            if 'N' not in motif:  # Skip motifs with N
                # Check for tandem repeats
                count = 1
                pos = i + repeat_len
                while pos + repeat_len <= len(seq) and seq[pos:pos+repeat_len] == motif:
                    count += 1
                    pos += repeat_len

                if count >= 3:  # At least 3 repeats
                    repeats[motif] += count

    return dict(repeats)

# Process FASTA file for extended statistics
input_file = "${fasta}"
stats_file = "${prefix}_extended_stats.json"
report_file = "${prefix}_quality_report.html"
log_file = "${prefix}_qc_log.txt"

sequences = {}
sequence_lengths = []
total_length = 0
chromosome_data = {}

print(f"🔍 Analyzing {input_file}")

# Parse FASTA and collect detailed statistics
with open(input_file, 'r') as f:
    current_seq = []
    current_id = None

    for line_num, line in enumerate(f, 1):
        line = line.strip()
        if line.startswith('>'):
            # Process previous sequence
            if current_id and current_seq:
                seq_data = ''.join(current_seq)
                seq_len = len(seq_data)

                # Detailed analysis
                composition = analyze_sequence_composition(seq_data)

                # Simple repeat detection (sample first 10kb for performance)
                sample_seq = seq_data[:10000] if len(seq_data) > 10000 else seq_data
                repeats = detect_repeats(sample_seq)

                sequences[current_id] = {
                    'length': seq_len,
                    'composition': composition,
                    'repeats': repeats,
                    'is_chromosome': not current_id.startswith('scaffold_')
                }

                sequence_lengths.append(seq_len)
                total_length += seq_len

                if sequences[current_id]['is_chromosome']:
                    chromosome_data[current_id] = sequences[current_id]

            # Start new sequence
            current_id = line[1:].split()[0]
            current_seq = []
        else:
            current_seq.append(line)

    # Process last sequence
    if current_id and current_seq:
        seq_data = ''.join(current_seq)
        seq_len = len(seq_data)
        composition = analyze_sequence_composition(seq_data)
        sample_seq = seq_data[:10000] if len(seq_data) > 10000 else seq_data
        repeats = detect_repeats(sample_seq)

        sequences[current_id] = {
            'length': seq_len,
            'composition': composition,
            'repeats': repeats,
            'is_chromosome': not current_id.startswith('scaffold_')
        }

        sequence_lengths.append(seq_len)
        total_length += seq_len

        if sequences[current_id]['is_chromosome']:
            chromosome_data[current_id] = sequences[current_id]

print(f"📊 Processed {len(sequences)} sequences")

# Calculate assembly statistics
if sequence_lengths:
    n50 = calculate_n50(sequence_lengths)
    n90 = calculate_n50([l for l in sorted(sequence_lengths, reverse=True)[:int(len(sequence_lengths)*0.9)]])

    # L50: number of sequences needed to reach N50
    sorted_lengths = sorted(sequence_lengths, reverse=True)
    cumulative = 0
    l50 = 0
    for length in sorted_lengths:
        cumulative += length
        l50 += 1
        if cumulative >= total_length / 2:
            break
else:
    n50 = n90 = l50 = 0

# Genome-wide composition
if total_length > 0:
    total_composition = {
        'A': sum(s['composition']['A'] for s in sequences.values()),
        'T': sum(s['composition']['T'] for s in sequences.values()),
        'G': sum(s['composition']['G'] for s in sequences.values()),
        'C': sum(s['composition']['C'] for s in sequences.values()),
        'N': sum(s['composition']['N'] for s in sequences.values()),
    }

    genome_gc = ((total_composition['G'] + total_composition['C']) / total_length) * 100
    genome_n = (total_composition['N'] / total_length) * 100
else:
    total_composition = {'A': 0, 'T': 0, 'G': 0, 'C': 0, 'N': 0}
    genome_gc = genome_n = 0

# Quality tier classification
def classify_quality_tier(stats):
    \"\"\"Classify genome quality into tiers A+ to C\"\"\"
    score = 0
    max_score = 100

    # Chromosome completeness (30 points)
    chrom_count = len([s for s in stats['sequences'] if stats['sequences'][s]['is_chromosome']])
    expected_chroms = 27  # 26 autosomes + X
    chrom_completeness = min(chrom_count / expected_chroms, 1.0)
    score += chrom_completeness * 30

    # Assembly contiguity (25 points)
    if stats['assembly_stats']['n50'] > 50_000_000:  # >50Mb
        score += 25
    elif stats['assembly_stats']['n50'] > 10_000_000:  # >10Mb
        score += 20
    elif stats['assembly_stats']['n50'] > 1_000_000:   # >1Mb
        score += 15
    elif stats['assembly_stats']['n50'] > 100_000:    # >100kb
        score += 10
    else:
        score += 5

    # Sequence quality (20 points)
    n_content = stats['genome_composition']['n_content']
    if n_content < 1.0:
        score += 20
    elif n_content < 2.0:
        score += 15
    elif n_content < 5.0:
        score += 10
    elif n_content < 10.0:
        score += 5

    # Genome size appropriateness (15 points)
    genome_size = stats['assembly_stats']['total_length']
    if 2.4e9 <= genome_size <= 3.2e9:  # Sheep genome range
        score += 15
    elif 2.0e9 <= genome_size <= 3.8e9:  # Extended range
        score += 10
    elif 1.5e9 <= genome_size <= 4.5e9:  # Very extended
        score += 5

    # GC content appropriateness (10 points)
    gc_content = stats['genome_composition']['gc_content']
    if 40 <= gc_content <= 45:  # Optimal for sheep
        score += 10
    elif 35 <= gc_content <= 50:  # Acceptable
        score += 8
    elif 30 <= gc_content <= 55:  # Marginal
        score += 5

    # Assign tier
    score_percent = (score / max_score) * 100

    if score_percent >= 90:
        return "A+", score_percent
    elif score_percent >= 80:
        return "A", score_percent
    elif score_percent >= 70:
        return "B+", score_percent
    elif score_percent >= 60:
        return "B", score_percent
    elif score_percent >= 50:
        return "B-", score_percent
    elif score_percent >= 40:
        return "C+", score_percent
    elif score_percent >= 30:
        return "C", score_percent
    else:
        return "C-", score_percent

# Compile comprehensive statistics
extended_stats = {
    "sample_id": "${prefix}",
    "analysis_timestamp": datetime.datetime.now().isoformat(),
    "input_file": "${fasta}",

    "assembly_stats": {
        "total_sequences": len(sequences),
        "total_length": total_length,
        "chromosome_sequences": len(chromosome_data),
        "scaffold_sequences": len(sequences) - len(chromosome_data),
        "longest_sequence": max(sequence_lengths) if sequence_lengths else 0,
        "shortest_sequence": min(sequence_lengths) if sequence_lengths else 0,
        "mean_length": total_length / len(sequences) if sequences else 0,
        "median_length": sorted(sequence_lengths)[len(sequence_lengths)//2] if sequence_lengths else 0,
        "n50": n50,
        "n90": n90,
        "l50": l50
    },

    "genome_composition": {
        **total_composition,
        "gc_content": genome_gc,
        "at_content": 100 - genome_gc - genome_n,
        "n_content": genome_n,
        "total_bases": sum(total_composition.values())
    },

    "chromosome_analysis": {
        "expected_chromosomes": 27,
        "found_chromosomes": len(chromosome_data),
        "completeness_percent": (len(chromosome_data) / 27) * 100,
        "chromosome_details": {
            name: {
                "length": data["length"],
                "gc_content": data["composition"]["gc_content"],
                "n_content": data["composition"]["n_content"]
            } for name, data in chromosome_data.items()
        }
    },

    "sequences": sequences
}

# Quality classification
tier, score = classify_quality_tier(extended_stats)
extended_stats["quality_assessment"] = {
    "tier": tier,
    "score": score,
    "suitable_for_pangenome": tier in ["A+", "A", "B+", "B"],
    "recommended_use": {
        "A+": "Optimal reference candidate",
        "A": "Excellent for pangenome construction",
        "B+": "Good quality, suitable for inclusion",
        "B": "Acceptable quality",
        "B-": "Marginal quality, use with caution",
        "C+": "Poor quality, limited utility",
        "C": "Very poor quality",
        "C-": "Unsuitable for pangenome analysis"
    }.get(tier, "Unknown")
}

# Save detailed JSON statistics
with open(stats_file, 'w') as f:
    json.dump(extended_stats, f, indent=2)

print(f"✅ Extended statistics saved to {stats_file}")

# Generate HTML quality report
html_content = f'''
<!DOCTYPE html>
<html>
<head>
    <title>Genome Quality Report - {extended_stats["sample_id"]}</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 20px; background-color: #f5f5f5; }}
        .container {{ max-width: 1200px; margin: 0 auto; background: white; padding: 20px; border-radius: 10px; }}
        .header {{ background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; padding: 20px; border-radius: 5px; }}
        .tier-{tier.lower().replace('+', 'plus').replace('-', 'minus')} {{
            background: {'linear-gradient(135deg, #4CAF50, #45a049)' if tier.startswith('A') else
                        'linear-gradient(135deg, #2196F3, #1976D2)' if tier.startswith('B') else
                        'linear-gradient(135deg, #FF9800, #F57C00)'};
        }}
        .stats-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(300px, 1fr)); gap: 20px; margin: 20px 0; }}
        .stat-card {{ background: #f8f9fa; padding: 15px; border-radius: 5px; border-left: 4px solid #007bff; }}
        .chromosome-grid {{ display: grid; grid-template-columns: repeat(auto-fill, minmax(200px, 1fr)); gap: 10px; }}
        .chromosome-card {{ background: #e8f5e8; padding: 10px; border-radius: 3px; font-size: 0.9em; }}
        .quality-badge {{ display: inline-block; padding: 5px 15px; border-radius: 20px; color: white; font-weight: bold; }}
        h1, h2 {{ margin-top: 0; }}
        .metric {{ font-size: 1.2em; font-weight: bold; color: #333; }}
        .warning {{ background: #fff3cd; border: 1px solid #ffeaa7; padding: 10px; border-radius: 5px; }}
        .success {{ background: #d4edda; border: 1px solid #c3e6cb; padding: 10px; border-radius: 5px; }}
        table {{ width: 100%; border-collapse: collapse; margin: 10px 0; }}
        th, td {{ padding: 8px 12px; text-align: left; border-bottom: 1px solid #ddd; }}
        th {{ background-color: #f2f2f2; }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header tier-{tier.lower().replace('+', 'plus').replace('-', 'minus')}">
            <h1>🐑 Sheep Genome Quality Report</h1>
            <h2>Sample: {extended_stats["sample_id"]}</h2>
            <p>Quality Tier: <span class="quality-badge">{tier}</span> (Score: {score:.1f}/100)</p>
            <p>{extended_stats["quality_assessment"]["recommended_use"]}</p>
        </div>

        <div class="stats-grid">
            <div class="stat-card">
                <h3>📏 Assembly Statistics</h3>
                <div class="metric">{extended_stats["assembly_stats"]["total_sequences"]:,}</div>
                <div>Total Sequences</div><br>
                <div class="metric">{extended_stats["assembly_stats"]["total_length"]/1e9:.2f} Gb</div>
                <div>Total Length</div><br>
                <div class="metric">{extended_stats["assembly_stats"]["n50"]/1e6:.1f} Mb</div>
                <div>N50</div>
            </div>

            <div class="stat-card">
                <h3>🧬 Composition</h3>
                <div class="metric">{extended_stats["genome_composition"]["gc_content"]:.1f}%</div>
                <div>GC Content</div><br>
                <div class="metric">{extended_stats["genome_composition"]["n_content"]:.2f}%</div>
                <div>N Content</div><br>
                <div class="metric">{extended_stats["genome_composition"]["total_bases"]:,}</div>
                <div>Total Bases</div>
            </div>

            <div class="stat-card">
                <h3>🧩 Chromosomes</h3>
                <div class="metric">{extended_stats["chromosome_analysis"]["found_chromosomes"]}/27</div>
                <div>Chromosomes Found</div><br>
                <div class="metric">{extended_stats["chromosome_analysis"]["completeness_percent"]:.1f}%</div>
                <div>Completeness</div><br>
                <div class="metric">{extended_stats["assembly_stats"]["scaffold_sequences"]:,}</div>
                <div>Scaffolds</div>
            </div>
        </div>

        <h3>📊 Detailed Assembly Metrics</h3>
        <table>
            <tr><th>Metric</th><th>Value</th><th>Assessment</th></tr>
            <tr>
                <td>Genome Size</td>
                <td>{extended_stats["assembly_stats"]["total_length"]/1e9:.2f} Gb</td>
                <td>{'✅ Within sheep range' if 2.4e9 <= extended_stats["assembly_stats"]["total_length"] <= 3.2e9 else '⚠️ Outside typical range'}</td>
            </tr>
            <tr>
                <td>N50</td>
                <td>{extended_stats["assembly_stats"]["n50"]/1e6:.1f} Mb</td>
                <td>{'✅ Excellent' if extended_stats["assembly_stats"]["n50"] > 50e6 else '✅ Good' if extended_stats["assembly_stats"]["n50"] > 10e6 else '⚠️ Fragmented'}</td>
            </tr>
            <tr>
                <td>GC Content</td>
                <td>{extended_stats["genome_composition"]["gc_content"]:.1f}%</td>
                <td>{'✅ Optimal' if 40 <= extended_stats["genome_composition"]["gc_content"] <= 45 else '✅ Acceptable' if 35 <= extended_stats["genome_composition"]["gc_content"] <= 50 else '⚠️ Unusual'}</td>
            </tr>
            <tr>
                <td>Gap Content (N%)</td>
                <td>{extended_stats["genome_composition"]["n_content"]:.2f}%</td>
                <td>{'✅ Excellent' if extended_stats["genome_composition"]["n_content"] < 1 else '✅ Good' if extended_stats["genome_composition"]["n_content"] < 5 else '⚠️ High'}</td>
            </tr>
        </table>

        <h3>🧩 Chromosome Analysis</h3>
        <p>Expected: 26 autosomes (chr1-chr26) + X chromosome + mitochondrial (chrMT)</p>
        <div class="chromosome-grid">
'''

# Add chromosome cards
for i in range(1, 27):
    chrom_name = f"chr{i}"
    if chrom_name in chromosome_data:
        html_content += f'<div class="chromosome-card">✅ {chrom_name}<br>{chromosome_data[chrom_name]["length"]/1e6:.1f} Mb</div>'
    else:
        html_content += f'<div class="chromosome-card" style="background:#ffebee;">❌ {chrom_name}<br>Missing</div>'

# Add X and MT
for chrom_name in ['chrX', 'chrMT']:
    if chrom_name in chromosome_data:
        html_content += f'<div class="chromosome-card">✅ {chrom_name}<br>{chromosome_data[chrom_name]["length"]/1e6:.1f} Mb</div>'
    else:
        html_content += f'<div class="chromosome-card" style="background:#ffebee;">❌ {chrom_name}<br>Missing</div>'

html_content += f'''
        </div>

        <h3>🔬 Quality Assessment Summary</h3>
        <div class="{'success' if extended_stats["quality_assessment"]["suitable_for_pangenome"] else 'warning'}">
            <strong>Pangenome Suitability: {'✅ SUITABLE' if extended_stats["quality_assessment"]["suitable_for_pangenome"] else '⚠️ LIMITED SUITABILITY'}</strong><br>
            Tier {tier} - {extended_stats["quality_assessment"]["recommended_use"]}<br>
            Overall Score: {score:.1f}/100
        </div>

        <p><small>Generated: {extended_stats["analysis_timestamp"]}</small></p>
    </div>
</body>
</html>
'''

with open(report_file, 'w') as f:
    f.write(html_content)

# Generate log file
with open(log_file, 'w') as f:
    f.write(f"Extended Quality Control Report - {extended_stats['sample_id']}\\n")
    f.write("="*60 + "\\n")
    f.write(f"Analysis completed: {extended_stats['analysis_timestamp']}\\n")
    f.write(f"Input file: {input_file}\\n")
    f.write(f"Quality tier: {tier} (Score: {score:.1f}/100)\\n")
    f.write(f"Suitable for pangenome: {extended_stats['quality_assessment']['suitable_for_pangenome']}\\n\\n")

    f.write("ASSEMBLY STATISTICS:\\n")
    f.write(f"  Total sequences: {extended_stats['assembly_stats']['total_sequences']:,}\\n")
    f.write(f"  Total length: {extended_stats['assembly_stats']['total_length']:,} bp ({extended_stats['assembly_stats']['total_length']/1e9:.2f} Gb)\\n")
    f.write(f"  N50: {extended_stats['assembly_stats']['n50']:,} bp ({extended_stats['assembly_stats']['n50']/1e6:.1f} Mb)\\n")
    f.write(f"  L50: {extended_stats['assembly_stats']['l50']:,}\\n")
    f.write(f"  Longest sequence: {extended_stats['assembly_stats']['longest_sequence']:,} bp\\n")
    f.write(f"  Mean length: {extended_stats['assembly_stats']['mean_length']:,.0f} bp\\n\\n")

    f.write("COMPOSITION ANALYSIS:\\n")
    f.write(f"  GC content: {extended_stats['genome_composition']['gc_content']:.2f}%\\n")
    f.write(f"  AT content: {extended_stats['genome_composition']['at_content']:.2f}%\\n")
    f.write(f"  N content: {extended_stats['genome_composition']['n_content']:.2f}%\\n\\n")

    f.write("CHROMOSOME ANALYSIS:\\n")
    f.write(f"  Expected chromosomes: 27 (chr1-26, chrX, chrMT)\\n")
    f.write(f"  Found chromosomes: {extended_stats['chromosome_analysis']['found_chromosomes']}\\n")
    f.write(f"  Completeness: {extended_stats['chromosome_analysis']['completeness_percent']:.1f}%\\n\\n")

    f.write("RECOMMENDATIONS:\\n")
    if extended_stats['quality_assessment']['suitable_for_pangenome']:
        f.write(f"  ✅ This genome is suitable for pangenome construction\\n")
        f.write(f"  ✅ Quality tier {tier} indicates {extended_stats['quality_assessment']['recommended_use'].lower()}\\n")
    else:
        f.write(f"  ⚠️  This genome has limitations for pangenome construction\\n")
        f.write(f"  ⚠️  Consider quality improvement or use with caution\\n")

print(f"✅ Quality control analysis completed")
print(f"   Quality Tier: {tier} (Score: {score:.1f}/100)")
print(f"   Report saved: {report_file}")
print(f"   Log saved: {log_file}")

EOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | cut -d' ' -f2)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_extended_stats.json
    touch ${prefix}_quality_report.html
    touch ${prefix}_qc_log.txt
    touch versions.yml
    """
}