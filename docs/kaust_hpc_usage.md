# KAUST Ibex HPC Usage Guide

## Overview

This pipeline has been optimized for **KAUST Ibex cluster**. The SLURM configuration follows Ibex-specific resource limits, partition names, and job submission policies.

## KAUST Ibex Cluster

### Ibex Specifications
- **Default partition**: `batch`
- **Memory per node**: ~370GB maximum
- **Max wall time**: 14 days (336 hours)
- **Job limits**: 2000 running jobs, 1300 CPUs, 24 GPUs
- **Default memory**: 2GB per CPU requested
- **Work directory**: `/ibex/user/$USER/nextflow-work`
- **Node preference**: Intel nodes via `--constraint=intel`

## Configuration

### Environment Setup

Set your KAUST account before running:
```bash
export SLURM_ACCOUNT="your_project_account"
# Example: export SLURM_ACCOUNT="k1234"
```

### Pipeline Execution

The pipeline is configured specifically for Ibex:

```bash
nextflow run . --input samples.csv -profile slurm,singularity
```

## Usage Examples

### Stage 1 on Ibex
```bash
# Set your account
export SLURM_ACCOUNT="k1234"

# Run Stage 1
nextflow run . \
    --input sheep_samples.csv \
    --outdir /ibex/user/$USER/sheep-pggb-results \
    --stage 1 \
    -profile slurm,singularity \
    -work-dir /ibex/user/$USER/nextflow-work
```

### Large Scale Analysis (Stage 3 PGGB)
```bash
# For memory-intensive PGGB on Ibex
nextflow run . \
    --input large_sheep_dataset.csv \
    --stage 3 \
    -profile slurm,singularity \
    --max_memory 370.GB \
    --max_time 336.h
```

## Resource Optimization

### Stage-specific Resource Allocation

The pipeline automatically optimizes resources per stage:

| Stage | Partition | CPUs | Memory | Time | Notes |
|-------|-----------|------|--------|------|-------|
| 1: Download | batch | 2 | 4GB | 30m | Network downloads |
| 1: Validate | batch | 1 | 2GB | 1h | Genome validation |
| 3: PGGB | batch | 32 | 150-370GB | 48-336h | High-memory |
| 4: Analysis | batch | 4-8 | 16-32GB | 6-12h | Graph analysis |
| 5: Variants | batch | 4-8 | 16-64GB | 8-24h | Variant calling |
| 6: Population | batch | 8-16 | 32-64GB | 24-48h | Pop genomics |

### Memory Scaling

The pipeline respects Ibex memory limits:
- **Ibex**: Maximum 370GB per node
- **Auto-scaling**: Memory increases on job retries within limits
- **Intel nodes**: Preferred for consistent performance

### Time Limits

Jobs are configured within Ibex time limits:
- **Ibex**: Maximum 14 days (336 hours)
- **Automatic scaling**: Time increases on retries within limits

## Ibex-Specific Features

### Account Management
```bash
# Default account fallback
clusterOptions = '--account=${SLURM_ACCOUNT ?: "k01"}'

# Custom account
export SLURM_ACCOUNT="your_specific_account"
```

### Partition Selection
```bash
# Ibex: Uses batch partition (default)
queue = 'batch'
```

### Intel Node Preference
```bash
# Prefer Intel nodes on Ibex for consistent performance
clusterOptions = '--account=${SLURM_ACCOUNT ?: "k01"} --constraint=intel'
```

### Storage Optimization
```bash
# Ibex user space (1.5TB limit)
workDir = "/ibex/user/$USER/nextflow-work"
```

## Job Monitoring

### Check Job Status
```bash
# Check your jobs
squeue -u $USER

# Detailed job information
scontrol show job <job_id>

# Job efficiency report
seff <job_id>
```

### Monitor Resource Usage
```bash
# Check account usage
sshare -U $USER

# Node information
sinfo -N -l

# Partition information
sinfo -p batch
```

## Troubleshooting

### Common Ibex Issues

**1. Account Not Set**
```bash
# Error: Invalid account specified
# Solution: Set your SLURM account
export SLURM_ACCOUNT="k1234"
```

**2. Memory Limit Exceeded**
```bash
# Error: Job killed for exceeding memory
# Solution: Check job memory usage and adjust
seff <job_id>
# Reduce genome batch size or increase memory limits
```

**3. Time Limit Exceeded**
```bash
# Error: Job killed for exceeding time limit
# Solution: Use longer time limits (up to 14 days on Ibex)
--max_time 336.h
```

**4. Queue Limits Reached**
```bash
# Error: Too many jobs in queue
# Solution: Use job arrays more efficiently
# Pipeline automatically respects Ibex job limits: 2000 jobs/user
```

### Storage Issues

**1. Disk Quota Exceeded**
```bash
# Check Ibex quotas
df -h /ibex/user/$USER

# Clean old work directories
find /ibex/user/$USER/nextflow-work -name "work" -mtime +7 -exec rm -rf {} +

# Use cleanup option to remove intermediate files
nextflow run . --input samples.csv -with-cleanup
```

## Best Practices for Ibex

### Resource Planning

1. **Start Small**: Test with 3-5 genomes before scaling
2. **Stage Progression**: Complete each stage before proceeding
3. **Resource Monitoring**: Check `seff` reports for optimization
4. **Cleanup**: Remove old work directories regularly

### Optimization Tips

- **Use Intel nodes**: More consistent performance for genomics workloads
- **Monitor memory usage**: Stay within 370GB per node limit
- **Batch similar jobs**: Use job arrays for parallel processing

### Account Management

```bash
# Add to your ~/.bashrc
export SLURM_ACCOUNT="your_default_account"

# Project-specific accounts
export SLURM_ACCOUNT="project_specific_account"
```

## Support

### KAUST HPC Support
- **Email**: help@hpc.kaust.edu.sa
- **Documentation**: https://docs.hpc.kaust.edu.sa/
- **Office Hours**: Check KAUST HPC website

### Pipeline-Specific Issues
- **Check logs**: `ls work/*/*/.command.log`
- **Validate syntax**: `nextflow inspect main.nf`
- **Test configuration**: Run with `--profile test` first

## Example Submission Scripts

### Ibex Submission Script
```bash
#!/bin/bash
#SBATCH --job-name=sheep-pangenome
#SBATCH --partition=batch
#SBATCH --account=k1234
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --output=sheep-pangenome-%j.out
#SBATCH --error=sheep-pangenome-%j.err

# Load modules
module load nextflow singularity

# Set environment
export SLURM_ACCOUNT="k1234"
export NXF_SINGULARITY_CACHEDIR="/ibex/user/$USER/.singularity"

# Run pipeline
nextflow run . \
    --input sheep_samples.csv \
    --outdir results \
    --stage 1 \
    -profile slurm,singularity \
    -work-dir /ibex/user/$USER/nextflow-work
```

### High-Memory PGGB Submission Script (Ibex)
```bash
#!/bin/bash
#SBATCH --job-name=sheep-pggb
#SBATCH --partition=batch
#SBATCH --account=k1234
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=168:00:00
#SBATCH --constraint=intel
#SBATCH --output=sheep-pggb-%j.out
#SBATCH --error=sheep-pggb-%j.err

# Load modules
module load nextflow singularity

# Set environment
export SLURM_ACCOUNT="k1234"
export NXF_SINGULARITY_CACHEDIR="/ibex/user/$USER/.singularity"

# Run pipeline
nextflow run . \
    --input large_sheep_dataset.csv \
    --outdir /ibex/user/$USER/results \
    --stage 3 \
    --max_memory 370.GB \
    -profile slurm,singularity \
    -work-dir /ibex/user/$USER/nextflow-work
```

---

This configuration ensures optimal performance on KAUST Ibex cluster while respecting all cluster policies and resource limits.