/*
========================================================================================
    Reference Genome Selection Module
========================================================================================
    Description: Select optimal reference genome for PGGB pangenome construction
    Algorithm: Multi-criteria scoring incorporating quality, completeness, and suitability
========================================================================================
*/

process SELECT_REFERENCE {
    label 'process_low'

    container 'python:3.9-slim'

    publishDir "${params.outdir}/02_preprocessing/reference_selection", mode: params.publish_dir_mode

    input:
    val(combined_genome_data)  // List of [meta, genome.fa, stats.json, metadata.json]

    output:
    tuple val(selected_meta), path("selected_reference.fa")  , emit: reference_genome
    path "reference_metadata.json"                          , emit: reference_metadata
    path "reference_selection_report.json"                  , emit: selection_report
    path "versions.yml"                                     , emit: versions

    script:
    """
    #!/usr/bin/env python3

    import json
    import sys
    import shutil
    from pathlib import Path
    from collections import defaultdict

    def load_json_file(json_file):
        \"\"\"Safely load JSON file\"\"\"
        try:
            if Path(json_file).exists():
                with open(json_file, 'r') as f:
                    return json.load(f)
            else:
                print(f"Warning: JSON file not found: {json_file}")
                return {}
        except Exception as e:
            print(f"Error loading JSON file {json_file}: {e}")
            return {}

    def calculate_reference_score(meta, stats, metadata):
        \"\"\"Calculate reference suitability score using multiple criteria\"\"\"

        score = 0
        score_breakdown = {}

        # Criterion 1: Assembly Quality (25 points)
        assembly_quality = 0

        # Assembly level preference
        assembly_level = metadata.get('assemblyInfo', {}).get('assemblyLevel', '').lower()
        if assembly_level == 'chromosome':
            assembly_quality += 15
            score_breakdown['assembly_level'] = 15
        elif assembly_level == 'scaffold':
            assembly_quality += 10
            score_breakdown['assembly_level'] = 10
        elif assembly_level == 'contig':
            assembly_quality += 5
            score_breakdown['assembly_level'] = 5
        else:
            score_breakdown['assembly_level'] = 0

        # Scaffold N50 assessment
        scaffold_n50 = stats.get('assembly_metrics', {}).get('scaffold_n50', 0)
        if scaffold_n50 >= 50e6:
            assembly_quality += 10
            n50_score = 10
        elif scaffold_n50 >= 20e6:
            assembly_quality += 7
            n50_score = 7
        elif scaffold_n50 >= 10e6:
            assembly_quality += 5
            n50_score = 5
        elif scaffold_n50 >= 1e6:
            assembly_quality += 2
            n50_score = 2
        else:
            n50_score = 0

        score_breakdown['scaffold_n50'] = n50_score
        score += assembly_quality

        # Criterion 2: Completeness (25 points)
        completeness_score = 0

        # BUSCO completeness (if available)
        busco_completeness = 0
        if hasattr(meta, 'busco_complete') and meta.busco_complete is not None:
            busco_completeness = float(meta.busco_complete)
        elif 'quality_assessment' in stats:
            # Estimate from quality assessment
            quality_score = stats.get('quality_assessment', {}).get('chromosome_completeness', 0)
            busco_completeness = min(100, quality_score * 1.1)  # Conservative estimate

        if busco_completeness >= 95:
            completeness_score += 20
            busco_score = 20
        elif busco_completeness >= 90:
            completeness_score += 15
            busco_score = 15
        elif busco_completeness >= 80:
            completeness_score += 10
            busco_score = 10
        elif busco_completeness >= 70:
            completeness_score += 5
            busco_score = 5
        else:
            busco_score = 0

        score_breakdown['busco_completeness'] = busco_score

        # Chromosome-level scaffolds
        chr_scaffolds = stats.get('sequence_size_distribution', {}).get('very_large', 0)
        if chr_scaffolds >= 27:  # All sheep chromosomes
            completeness_score += 5
            chr_score = 5
        elif chr_scaffolds >= 20:
            completeness_score += 3
            chr_score = 3
        elif chr_scaffolds >= 15:
            completeness_score += 1
            chr_score = 1
        else:
            chr_score = 0

        score_breakdown['chromosome_coverage'] = chr_score
        score += completeness_score

        # Criterion 3: Sequence Quality (20 points)
        quality_score = 0

        # Gap content (N percentage)
        n_percentage = stats.get('nucleotide_composition', {}).get('n_percentage', 5.0)
        if n_percentage <= 0.5:
            quality_score += 10
            gap_score = 10
        elif n_percentage <= 1.0:
            quality_score += 8
            gap_score = 8
        elif n_percentage <= 2.0:
            quality_score += 6
            gap_score = 6
        elif n_percentage <= 3.0:
            quality_score += 4
            gap_score = 4
        elif n_percentage <= 5.0:
            quality_score += 2
            gap_score = 2
        else:
            gap_score = 0

        score_breakdown['gap_content'] = gap_score

        # Genome size accuracy (expected sheep genome ~2.8 Gb)
        total_length = stats.get('sample_info', {}).get('total_length', 0)
        expected_size = 2.8e9
        size_accuracy = max(0, 1 - abs(total_length - expected_size) / expected_size)

        if size_accuracy >= 0.95:
            quality_score += 10
            size_score = 10
        elif size_accuracy >= 0.90:
            quality_score += 8
            size_score = 8
        elif size_accuracy >= 0.85:
            quality_score += 6
            size_score = 6
        elif size_accuracy >= 0.80:
            quality_score += 4
            size_score = 4
        else:
            size_score = 0

        score_breakdown['size_accuracy'] = size_score
        score += quality_score

        # Criterion 4: Annotation and Metadata (15 points)
        annotation_score = 0

        # Annotation availability
        annotation_info = metadata.get('annotationInfo', {})
        if annotation_info:
            annotation_score += 10
            annot_score = 10
        else:
            annot_score = 0

        score_breakdown['annotation_available'] = annot_score

        # RefSeq status preference
        source_db = metadata.get('sourceDatabase', '').upper()
        if 'REFSEQ' in source_db:
            annotation_score += 5
            refseq_score = 5
        else:
            refseq_score = 0

        score_breakdown['refseq_status'] = refseq_score
        score += annotation_score

        # Criterion 5: Reference Suitability (15 points)
        suitability_score = 0

        # Breed/population representativeness
        breed = metadata.get('organism', {}).get('infraspecificNames', {}).get('breed', '')
        if breed:
            # Prefer widely-used reference breeds
            reference_breeds = ['texel', 'rambouillet', 'dorset', 'suffolk', 'merino']
            if breed.lower() in reference_breeds:
                suitability_score += 8
                breed_score = 8
            else:
                suitability_score += 5
                breed_score = 5
        else:
            breed_score = 0

        score_breakdown['breed_suitability'] = breed_score

        # Assembly method preference (newer methods preferred)
        assembly_method = metadata.get('assemblyInfo', {}).get('assemblyMethod', '').lower()
        if any(method in assembly_method for method in ['hifi', 'pacbio', 'nanopore', 'ont']):
            suitability_score += 4
            method_score = 4
        elif any(method in assembly_method for method in ['illumina', 'sanger']):
            suitability_score += 2
            method_score = 2
        else:
            method_score = 0

        score_breakdown['assembly_method'] = method_score

        # Recent assembly preference
        release_date = metadata.get('assemblyInfo', {}).get('releaseDate', '')
        if release_date:
            try:
                year = int(release_date.split('-')[0])
                if year >= 2020:
                    suitability_score += 3
                    date_score = 3
                elif year >= 2015:
                    suitability_score += 2
                    date_score = 2
                elif year >= 2010:
                    suitability_score += 1
                    date_score = 1
                else:
                    date_score = 0
            except:
                date_score = 0
        else:
            date_score = 0

        score_breakdown['assembly_date'] = date_score
        score += suitability_score

        return score, score_breakdown, {
            'busco_completeness': busco_completeness,
            'scaffold_n50': scaffold_n50,
            'gap_content': n_percentage,
            'total_length': total_length,
            'assembly_level': assembly_level,
            'breed': breed,
            'release_date': release_date
        }

    def select_best_reference(genome_data_list):
        \"\"\"Select the best reference genome from candidates\"\"\"

        candidates = []

        print(f"Evaluating {len(genome_data_list)} genome candidates for reference selection")

        # Score each candidate
        for genome_data in genome_data_list:
            meta, genome_file, stats_file, metadata_file = genome_data

            # Load data files
            stats = load_json_file(stats_file)
            metadata = load_json_file(metadata_file)

            # Calculate reference score
            score, score_breakdown, metrics = calculate_reference_score(meta, stats, metadata)

            candidate = {
                'meta': meta,
                'genome_file': genome_file,
                'sample_id': meta['id'],
                'selection_score': score,
                'score_breakdown': score_breakdown,
                'key_metrics': metrics,
                'stats_data': stats,
                'metadata_data': metadata
            }

            candidates.append(candidate)

            print(f"  {meta['id']}: Score {score}/100")
            print(f"    Assembly level: {metrics.get('assembly_level', 'unknown')}")
            print(f"    N50: {metrics.get('scaffold_n50', 0)/1e6:.1f} Mb")
            print(f"    Gap content: {metrics.get('gap_content', 0):.2f}%")

        # Sort by score (descending)
        candidates.sort(key=lambda x: x['selection_score'], reverse=True)

        if not candidates:
            raise ValueError("No valid candidates found for reference selection")

        selected = candidates[0]
        alternatives = candidates[1:5]  # Top 5 alternatives

        print(f"")
        print(f"🏆 Selected reference: {selected['sample_id']}")
        print(f"   Selection score: {selected['selection_score']}/100")
        print(f"   Key metrics:")
        print(f"     Assembly level: {selected['key_metrics'].get('assembly_level', 'unknown')}")
        print(f"     N50: {selected['key_metrics'].get('scaffold_n50', 0)/1e6:.1f} Mb")
        print(f"     BUSCO completeness: {selected['key_metrics'].get('busco_completeness', 0):.1f}%")
        print(f"     Gap content: {selected['key_metrics'].get('gap_content', 0):.2f}%")

        return selected, alternatives

    # Main processing
    try:
        genome_data_list = ${combined_genome_data}

        print(f"Starting reference genome selection process")
        print(f"Candidates received: {len(genome_data_list)}")

        # Select best reference
        selected_ref, alternatives = select_best_reference(genome_data_list)

        # Copy selected reference genome
        shutil.copy2(selected_ref['genome_file'], 'selected_reference.fa')

        # Create reference metadata
        reference_metadata = {
            'selected_reference': {
                'sample_id': selected_ref['sample_id'],
                'selection_score': selected_ref['selection_score'],
                'score_breakdown': selected_ref['score_breakdown'],
                'key_metrics': selected_ref['key_metrics'],
                'quality_tier': selected_ref.get('quality_tier', 'Unknown'),
                'selection_timestamp': __import__('datetime').datetime.now().isoformat(),
                'original_metadata': selected_ref['metadata_data']
            },
            'selection_criteria': {
                'assembly_quality': 25,
                'completeness': 25,
                'sequence_quality': 20,
                'annotation_metadata': 15,
                'reference_suitability': 15
            },
            'pangenome_suitability': 'Excellent' if selected_ref['selection_score'] >= 85 else 'Good' if selected_ref['selection_score'] >= 70 else 'Adequate'
        }

        with open('reference_metadata.json', 'w') as f:
            json.dump(reference_metadata, f, indent=2)

        # Create selection report
        selection_report = {
            'selection_summary': {
                'total_candidates': len(genome_data_list),
                'selected_reference': selected_ref['sample_id'],
                'selection_score': selected_ref['selection_score'],
                'selection_timestamp': reference_metadata['selected_reference']['selection_timestamp']
            },
            'selected_reference': {
                'sample_id': selected_ref['sample_id'],
                'selection_score': selected_ref['selection_score'],
                'score_breakdown': selected_ref['score_breakdown'],
                'assembly_level': selected_ref['key_metrics'].get('assembly_level', 'unknown'),
                'scaffold_n50': selected_ref['key_metrics'].get('scaffold_n50', 0),
                'busco_completeness': selected_ref['key_metrics'].get('busco_completeness', 0),
                'gap_content': selected_ref['key_metrics'].get('gap_content', 0),
                'genome_size': selected_ref['key_metrics'].get('total_length', 0),
                'breed': selected_ref['key_metrics'].get('breed', ''),
                'release_date': selected_ref['key_metrics'].get('release_date', ''),
                'quality_tier': selected_ref.get('quality_tier', 'Unknown')
            },
            'alternative_candidates': [
                {
                    'sample_id': alt['sample_id'],
                    'selection_score': alt['selection_score'],
                    'quality_tier': alt.get('quality_tier', 'Unknown'),
                    'key_metrics': alt['key_metrics']
                } for alt in alternatives
            ],
            'selection_criteria_weights': reference_metadata['selection_criteria']
        }

        with open('reference_selection_report.json', 'w') as f:
            json.dump(selection_report, f, indent=2)

        print(f"✅ Reference selection completed successfully")
        print(f"   Selected: {selected_ref['sample_id']} (score: {selected_ref['selection_score']}/100)")
        print(f"   Reference file: selected_reference.fa")
        print(f"   Metadata: reference_metadata.json")
        print(f"   Report: reference_selection_report.json")

        # Output selected meta for emit
        selected_meta = selected_ref['meta']

    except Exception as e:
        print(f"❌ Error in reference selection: {str(e)}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """
}