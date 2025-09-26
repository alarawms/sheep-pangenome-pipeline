/*
========================================================================================
    Normalize Headers Module
========================================================================================
    Description: Standardize FASTA sequence headers for consistent pangenome analysis
    Features: Clean headers, remove problematic characters, maintain traceability
========================================================================================
*/

process NORMALIZE_HEADERS {
    tag "$meta.id"
    label 'process_low'

    container 'biocontainers/biopython:1.79'

    publishDir "${params.outdir}/02_preprocessing/standardized", mode: params.publish_dir_mode,
        saveAs: { filename ->
            if (filename.endsWith('.fa')) "${meta.id}_standardized.fa"
            else null
        }

    publishDir "${params.outdir}/02_preprocessing/header_mapping", mode: params.publish_dir_mode,
        saveAs: { filename ->
            if (filename.endsWith('.tsv')) "${meta.id}_header_mapping.tsv"
            else null
        }

    input:
    tuple val(meta), path(renamed_genome)

    output:
    tuple val(meta), path("${meta.id}_normalized.fa")       , emit: normalized_genomes
    tuple val(meta), path("${meta.id}_header_mapping.tsv")  , emit: header_map
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

    def clean_sequence_header(header, sample_id):
        \"\"\"Clean and standardize FASTA headers\"\"\"

        # Extract the sequence identifier (first part)
        seq_id = header.split()[0] if header else f"unknown_{sample_id}"

        # Remove problematic characters for downstream tools
        # Keep only alphanumeric, underscore, hyphen, dot, and pipe
        clean_id = re.sub(r'[^a-zA-Z0-9_\\-\\.\\|]', '_', seq_id)

        # Remove multiple consecutive underscores
        clean_id = re.sub(r'_+', '_', clean_id)

        # Remove leading/trailing underscores
        clean_id = clean_id.strip('_')

        # Ensure ID is not empty
        if not clean_id:
            clean_id = f"seq_{sample_id}"

        # Ensure ID doesn't start with a number (some tools don't like this)
        if clean_id[0].isdigit():
            clean_id = f"seq_{clean_id}"

        # Limit length for compatibility (most tools handle up to 255 chars)
        if len(clean_id) > 200:
            # Keep the beginning and end, indicate truncation
            clean_id = clean_id[:95] + "_TRUNC_" + clean_id[-95:]

        return clean_id

    def extract_description_metadata(description):
        \"\"\"Extract useful metadata from FASTA description\"\"\"

        if not description:
            return {}

        metadata = {}

        # Common patterns in genome assembly descriptions
        patterns = {
            'length': r'length[=:]?(\\d+)',
            'coverage': r'coverage[=:]?([\\d\\.]+)',
            'gc': r'gc[=:]?([\\d\\.]+)',
            'n50': r'n50[=:]?(\\d+)',
            'chromosome': r'chromosome[=:]?([\\w\\d]+)',
            'scaffold': r'scaffold[=:]?([\\w\\d]+)',
            'contig': r'contig[=:]?([\\w\\d]+)',
            'assembly': r'assembly[=:]?([\\w\\d\\.]+)',
            'version': r'version[=:]?([\\w\\d\\.]+)',
            'organism': r'organism[=:]?([\\w\\s]+)',
            'strain': r'strain[=:]?([\\w\\d]+)',
            'breed': r'breed[=:]?([\\w\\s]+)',
        }

        desc_lower = description.lower()
        for key, pattern in patterns.items():
            match = re.search(pattern, desc_lower)
            if match:
                metadata[key] = match.group(1).strip()

        return metadata

    def create_standardized_description(clean_id, original_header, metadata, sample_id):
        \"\"\"Create standardized description with essential information\"\"\"

        desc_parts = [clean_id]

        # Add sample information
        desc_parts.append(f"sample:{sample_id}")

        # Add key metadata if available
        if 'length' in metadata:
            desc_parts.append(f"len:{metadata['length']}")
        if 'chromosome' in metadata:
            desc_parts.append(f"chr:{metadata['chromosome']}")
        if 'scaffold' in metadata:
            desc_parts.append(f"scaffold:{metadata['scaffold']}")

        # Add original identifier for traceability
        original_id = original_header.split()[0] if original_header else "unknown"
        if original_id != clean_id:
            desc_parts.append(f"orig:{original_id}")

        return " | ".join(desc_parts)

    def process_fasta_headers(input_file, output_file, mapping_file, sample_id):
        \"\"\"Process FASTA file to normalize headers\"\"\"

        print(f"Normalizing headers for sample: {sample_id}")
        print(f"Input file: {input_file}")

        normalized_sequences = []
        mapping_records = []

        # Statistics
        total_sequences = 0
        normalized_count = 0
        duplicate_ids = set()
        id_counter = {}

        # Process each sequence
        for record in SeqIO.parse(input_file, "fasta"):
            total_sequences += 1

            original_header = record.description
            original_id = record.id

            # Clean the sequence ID
            clean_id = clean_sequence_header(original_id, sample_id)

            # Handle potential duplicates
            if clean_id in id_counter:
                id_counter[clean_id] += 1
                clean_id_final = f"{clean_id}_{id_counter[clean_id]}"
                duplicate_ids.add(clean_id)
            else:
                id_counter[clean_id] = 0
                clean_id_final = clean_id

            # Extract metadata from description
            metadata = extract_description_metadata(original_header)

            # Create standardized description
            new_description = create_standardized_description(
                clean_id_final, original_header, metadata, sample_id
            )

            # Create normalized record
            normalized_record = SeqRecord(
                record.seq,
                id=clean_id_final,
                description=new_description
            )

            normalized_sequences.append(normalized_record)
            normalized_count += 1

            # Record mapping information
            mapping_record = {
                'original_id': original_id,
                'normalized_id': clean_id_final,
                'original_description': original_header,
                'normalized_description': new_description,
                'sequence_length': len(record.seq),
                'metadata_extracted': len(metadata),
                'was_modified': original_id != clean_id_final
            }

            mapping_records.append(mapping_record)

            # Print progress for significant modifications
            if original_id != clean_id_final:
                print(f"  Modified: {original_id} -> {clean_id_final}")

        # Write normalized sequences
        with open(output_file, 'w') as out_handle:
            SeqIO.write(normalized_sequences, out_handle, "fasta")

        # Write mapping file
        with open(mapping_file, 'w') as map_handle:
            map_handle.write("original_id\\tnormalized_id\\toriginal_description\\tnormalized_description\\tsequence_length\\tmetadata_extracted\\twas_modified\\n")
            for record in mapping_records:
                map_handle.write(f"{record['original_id']}\\t{record['normalized_id']}\\t{record['original_description']}\\t{record['normalized_description']}\\t{record['sequence_length']}\\t{record['metadata_extracted']}\\t{record['was_modified']}\\n")

        # Summary statistics
        modified_count = sum(1 for r in mapping_records if r['was_modified'])
        metadata_count = sum(r['metadata_extracted'] for r in mapping_records)

        print(f"")
        print(f"Header normalization summary for {sample_id}:")
        print(f"  Total sequences processed: {total_sequences}")
        print(f"  Headers modified: {modified_count}")
        print(f"  Headers unchanged: {total_sequences - modified_count}")
        print(f"  Duplicate IDs resolved: {len(duplicate_ids)}")
        print(f"  Metadata elements extracted: {metadata_count}")
        print(f"  Modification rate: {(modified_count/total_sequences*100):.1f}%")

        return {
            'total_sequences': total_sequences,
            'modified_headers': modified_count,
            'unchanged_headers': total_sequences - modified_count,
            'duplicate_ids_resolved': len(duplicate_ids),
            'metadata_extracted': metadata_count,
            'modification_rate': modified_count/total_sequences
        }

    # Main processing
    try:
        input_file = "${renamed_genome}"
        output_file = "${meta.id}_normalized.fa"
        mapping_file = "${meta.id}_header_mapping.tsv"
        sample_id = "${meta.id}"

        result = process_fasta_headers(input_file, output_file, mapping_file, sample_id)

        print(f"")
        print(f"✅ Successfully normalized headers for {sample_id}")
        print(f"   Output: {output_file}")
        print(f"   Mapping: {mapping_file}")

    except Exception as e:
        print(f"❌ Error processing {sample_id}: {str(e)}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        biopython: \$(python -c "import Bio; print(Bio.__version__)")
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """
}