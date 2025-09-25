# Stage 1: Data Acquisition & Preparation

## Overview

Stage 1 implements the critical foundation of the sheep pangenome pipeline: **data acquisition and preparation with comprehensive validation**. This stage ensures all input genomes meet quality standards before proceeding to computationally intensive downstream analyses.

## Scientific Objectives

- **Genome Collection**: Acquire 15-55+ sheep genomes from NCBI and local sources
- **Quality Assessment**: Validate genome completeness, size, and composition
- **Standardization**: Prepare genomes for downstream pangenome construction
- **Validation Gates**: Ensure only high-quality data proceeds to Stage 2

## Input Requirements

### Samplesheet Format

Create a CSV file with the following columns:

```csv
sample,accession,breed,population,geographic_origin
rambouillet_ref,GCF_016772045.1,Rambouillet,European,United_States
texel_sample,GCF_000298735.2,Texel,European,Netherlands
custom_breed,,/path/to/genome.fa,Custom_Breed,Regional,Farm_Location
```

**Required columns:**
- `sample`: Unique sample identifier
- Either `accession` (for NCBI download) OR `fasta` (for local files)

**Optional columns:**
- `breed`: Sheep breed name
- `population`: Population group (European, Asian, African, Wild)
- `geographic_origin`: Geographic location or country

### Mixed Input Support

Stage 1 supports both:
- **NCBI Downloads**: Automatic genome downloading via accession numbers
- **Local Files**: Pre-existing FASTA files on the filesystem
- **Combined**: Mix of downloaded and local genomes in one analysis

## Execution

### Quick Start

```bash
# Basic execution with Docker
nextflow run . --input sheep_samples.csv -profile docker

# SLURM cluster execution
nextflow run . --input sheep_samples.csv -profile slurm,singularity

# Test with small dataset
nextflow run . --input test_samples.csv -profile test,docker
```

### Advanced Options

```bash
nextflow run . \
    --input sheep_samples.csv \
    --outdir results/stage1 \
    --stage 1 \
    --max_download_time 45m \
    --download_retries 3 \
    --validation_strict true \
    -profile slurm,singularity \
    -work-dir /scratch/$USER/nextflow-work
```

## Parameters

### Core Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input` | null | Path to input samplesheet (CSV) |
| `--outdir` | ./results | Output directory |
| `--stage` | 1 | Pipeline stage to run |

### Download Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--max_download_time` | 30.m | Maximum time per genome download |
| `--download_retries` | 3 | Number of download retry attempts |

### Validation Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--validation_strict` | true | Strict validation mode |
| `--genome_size_min` | 2.4e9 | Minimum genome size (bp) |
| `--genome_size_max` | 3.2e9 | Maximum genome size (bp) |
| `--gc_content_min` | 35.0 | Minimum GC content (%) |
| `--gc_content_max` | 50.0 | Maximum GC content (%) |
| `--n_content_max` | 5.0 | Maximum N content (%) |
| `--max_contigs` | 50000 | Maximum number of contigs |

## Validation Criteria

### Sheep-Specific Quality Gates

Stage 1 applies biologically-informed validation:

1. **Size Validation**
   - Sheep genomes: 2.4-3.2 Gb
   - Flags outliers that may indicate assembly issues

2. **Composition Validation**
   - GC content: 35-50% (mammalian range)
   - N content: <5% (gap content)
   - Sequence count: <50,000 contigs

3. **Download Validation**
   - Successful NCBI retrieval
   - File integrity checks
   - Metadata completeness

### Success Metrics

Stage 1 completion requires:
- ✅ **100% download success** for NCBI accessions
- ✅ **≥95% validation pass rate** for all genomes
- ✅ **Complete metadata** for all samples
- ✅ **No duplicate samples** in dataset

## Output Structure

```
results/01_data_preparation/
├── downloaded_genomes/          # FASTA files from NCBI
│   ├── sample1.fa
│   └── sample2.fa
├── metadata/                    # Genome metadata
│   ├── sample1.json
│   └── sample2.json
├── validation/                  # Validation results
│   ├── sample1_validation.json
│   └── sample2_validation.json
├── statistics/                  # Genome statistics
│   ├── sample1_stats.txt
│   └── sample2_stats.txt
└── logs/                       # Download logs
    ├── download_log.txt
    └── pipeline_info/
```

## Quality Control Reports

### Validation Summary

Each genome produces:

1. **Validation JSON**: Pass/fail status with detailed metrics
2. **Statistics Report**: Comprehensive genome statistics
3. **Download Log**: Download success and timing information

### Example Validation Output

```json
{
  "sample_id": "rambouillet_ref",
  "validation_status": "PASS",
  "total_length": 2847123456,
  "total_sequences": 27,
  "gc_content": 42.1,
  "n_content": 2.3,
  "validation_messages": []
}
```

## Error Handling & Troubleshooting

### Common Issues

1. **Download Failures**
   ```bash
   # Check NCBI accession validity
   datasets summary genome accession GCF_016772045.1

   # Increase timeout for slow networks
   --max_download_time 60m
   ```

2. **Validation Failures**
   ```bash
   # Review validation results
   cat results/01_data_preparation/validation/*_validation.json

   # Use relaxed validation for testing
   -profile test
   ```

3. **Resource Issues**
   ```bash
   # Check SLURM job status
   squeue -u $USER

   # Increase memory for large genomes
   --max_memory 16.GB
   ```

## Next Steps: Stage 2 Preparation

Upon successful Stage 1 completion:

1. **Review Validation Summary**
   - Check that all samples passed validation
   - Verify expected genome count and diversity

2. **Prepare for Stage 2**
   - Stage 1 outputs become Stage 2 inputs
   - No manual intervention required for standard cases

3. **Optional Quality Checks**
   ```bash
   # Review genome statistics
   cat results/01_data_preparation/statistics/*.txt

   # Check validation summary
   grep -h "validation_status" results/01_data_preparation/validation/*.json
   ```

## Advanced Usage

### Custom Validation Criteria

For specialized analyses, modify validation parameters:

```bash
# Relaxed validation for ancient DNA or draft assemblies
nextflow run . \
    --input samples.csv \
    --genome_size_min 2.0e9 \
    --genome_size_max 4.0e9 \
    --gc_content_min 30.0 \
    --gc_content_max 55.0 \
    --n_content_max 15.0 \
    --validation_strict false
```

### Batch Processing

For large datasets (50+ genomes):

```bash
# Use SLURM job arrays for efficient parallel downloading
nextflow run . \
    --input large_dataset.csv \
    -profile slurm,singularity \
    --max_cpus 64 \
    --max_memory 128.GB \
    -work-dir /scratch/$USER/sheep-pggb
```

### Resume Failed Runs

Nextflow supports automatic resume:

```bash
# Resume from last successful checkpoint
nextflow run . --input samples.csv -resume

# Force fresh start
nextflow run . --input samples.csv -N  # or --fresh
```

## Stage 1 Completion Checklist

- [ ] All genomes downloaded successfully
- [ ] All validations passed (or reviewed failures)
- [ ] Output directory contains all expected files
- [ ] Validation summary shows acceptable success rate
- [ ] Ready to proceed to Stage 2: Genome Preprocessing

## Support

For Stage 1 issues:
1. Check validation reports in `results/01_data_preparation/validation/`
2. Review SLURM logs: `seff <job_id>`
3. Test with small dataset using `-profile test`
4. Validate NCBI accessions manually before pipeline execution