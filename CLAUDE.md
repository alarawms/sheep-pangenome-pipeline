# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This repository is for developing **nf-core/sheeppangenome**, a Nextflow pipeline for comprehensive sheep pangenome construction and comparative genomics analysis. The pipeline will process 50+ publicly available sheep genome assemblies using PGGB (pangenome graph builder) and perform population genomics analysis.

## Project Status

**Current State**: Planning/Documentation phase
- Only `project.md` exists with comprehensive implementation guide
- No actual pipeline code has been implemented yet
- Repository is initialized but empty (no commits)

## Architecture Plan

The planned pipeline follows nf-core DSL2 structure:

```
nf-core-sheeppangenome/
├── main.nf                     # Pipeline entry point
├── nextflow.config             # Main configuration with sheep-optimized parameters
├── workflows/sheeppangenome.nf # Primary workflow
├── subworkflows/local/         # Custom subworkflows
├── modules/local/              # Custom modules (PGGB, download_genome, etc.)
├── conf/                       # Configuration files (base, slurm, modules)
├── assets/                     # Genome catalog and validation schemas
└── docs/                       # Usage and integration documentation
```

## Key Components to Implement

### Core Modules
- **PGGB**: Pangenome graph construction optimized for 3GB sheep genomes
- **DOWNLOAD_GENOME**: Automatic genome downloading via NCBI datasets
- **ODGI_STATS**: Graph statistics and visualization
- **Population analysis modules**: ADMIXTURE, FST, phylogenetics

### Computational Requirements
- **Memory**: 150-300GB for PGGB depending on genome count
- **CPU**: 32+ cores recommended for PGGB parallelization
- **Storage**: 200GB-800GB+ depending on analysis scale
- **Time**: 24-72+ hours for complete analysis

## SLURM Optimization

The pipeline is designed for HPC clusters with SLURM:
- PGGB processes: 32 CPU, 150GB memory, 48h runtime
- Download processes: 2 CPU, 4GB memory, 2h runtime
- Variant calling: 4 CPU, 24GB memory, 12h runtime
- Work directory: `/scratch/$USER/nextflow-work`

## Development Commands

### Initial Setup
```bash
# Initialize nf-core pipeline structure
nf-core create sheeppangenome

# Install development dependencies
pip install nf-core nextflow
conda install -c bioconda pggb ncbi-datasets-cli
```

### Testing
```bash
# Test with small dataset (3-5 genomes)
nextflow run . --input small_test.csv -profile test,docker

# Validate pipeline structure
nf-core lint .

# Run pipeline tests
nf-core test .
```

### Execution
```bash
# Basic pangenome construction
nextflow run . --input sheep_pangenome.csv --outdir results -profile slurm,singularity

# Full analysis with population genomics
nextflow run . --input extended_sheep_set.csv --enable_variant_calling true --enable_population_analysis true -profile slurm,singularity
```

## Key Parameters

### PGGB Configuration (Sheep-Optimized)
- `segment_length`: 50000 (optimal for 3GB genomes)
- `block_id_min`: 0.95 (95% identity threshold for domestic breeds)
- `map_pct_id`: 0.95 (high-confidence alignments)
- `n_mappings`: 10 (multiple mappings per segment)
- `threads`: 32 (parallel processing)

### Resource Limits
- `max_memory`: 2000.GB
- `max_cpus`: 128
- `max_time`: 240.h

## Sample Input Format

Expected CSV format for mixed local/downloaded genomes:
```csv
sample,accession,fasta,breed,population,geographic_origin
reference,GCF_016772045.1,,Rambouillet,European,United_States
custom_breed,,/data/genomes/breed1.fa,Custom_Breed,Regional,Farm_A
```

## Development Priorities

1. **Core pipeline structure**: Implement main.nf, nextflow.config, and workflow files
2. **PGGB module**: Create optimized PGGB process with sheep-specific parameters
3. **Download module**: Implement NCBI datasets integration for automatic downloading
4. **SLURM configuration**: Set up HPC-optimized resource allocation
5. **Testing framework**: Create test datasets and validation procedures
6. **Documentation**: Convert project.md plans into proper nf-core documentation

## Special Considerations

- **Genome catalog**: Pipeline includes 55+ pre-cataloged sheep genome assemblies with quality metrics
- **Mixed input support**: Handle both local FASTA files and NCBI accession numbers
- **HPC optimization**: Pipeline designed specifically for SLURM cluster environments
- **AI integration**: Planned Claude Code integration for development assistance and result interpretation
- **Memory management**: PGGB requires careful memory allocation for large-scale analyses

## Quality Standards

- Follow nf-core template and best practices
- All modules must include proper error handling and resource requirements
- Container-based execution (Docker/Singularity) required
- Comprehensive testing with multiple genome scales (15, 35, 55+ genomes)
- Proper validation schemas for input files