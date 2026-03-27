# HPC Cluster Parallel Job Submission Guide

## Overview
The COJO conditional analysis now supports parallel job submission to HPC clusters. Each SNP/locus runs as a separate cluster job with 32GB memory, and all jobs are submitted simultaneously.

## Configuration

### Location
- **Profile directory**: `.smk/profiles/hpc/config.yaml`
- **Resource allocation**: Defined in `rules/cojo_condition.smk` 
  - Memory: 32 GB (RESOURCE_MEM_MB = 32000)
  - Wall time: 30 minutes (default, configurable)

### Running with Cluster Submission

#### Step 1: Extract regions (local execution)
```bash
./submit.sh --snakefile Snakefile --cores 2 extract_regions_all
```

#### Step 2: Run COJO with HPC cluster parallelization
```bash
./submit.sh --snakefile Snakefile --profile .smk/profiles/hpc cojo_condition
```

**What happens:**
- Snakemake builds DAG with all SNP/locus combinations
- Each locus COJO analysis submitted as separate cluster job
- All jobs queue simultaneously (up to `jobs: 50` limit)
- Each job: 32 GB memory, auto-determined CPU cores, 30 min wall time
- Snakemake waits for all jobs to complete

### Customizing for Your HPC System

#### SLURM (default)
The profile uses SLURM `sbatch` by default. No changes needed if your cluster uses SLURM.

#### PBS/Torque Systems
Edit `.smk/profiles/hpc/config.yaml`:
```yaml
cluster: "qsub -N {rule}-{wildcards} -l vmem={resources.mem_mb}mb -lselect=1:ncpus={resources.cores}:walltime={resources.time} -j oe -o {log}.pbs.log"
```

#### LSF Systems
Edit `.smk/profiles/hpc/config.yaml`:
```yaml
cluster: "bsub -J {rule}-{wildcards} -M {resources.mem_mb} -n {resources.cores} -W {resources.time} -o {log}.lsf.log"
```

### Adjusting Parallelism

#### Increase/Decrease Parallel Jobs
Edit `.smk/profiles/hpc/config.yaml`:
```yaml
jobs: 50  # Change to number of simultaneous cluster jobs (default: 50)
```

**Example values:**
- `jobs: 25` → Submit max 25 jobs at once
- `jobs: 50` → Submit max 50 jobs at once  
- `jobs: 100` → Submit all at once (if < 100 SNPs)

**Recommendation**: Set `jobs` slightly higher than typical queue capacity to let the scheduler manage the queue.

### Adjusting Memory/Wall-time

#### Per-locus memory
Edit `.smk/profiles/hpc/config.yaml` → update memory in `cluster:` command, OR
Edit `configs/analysis.yml`:
```yaml
target_analyses:
  Hispanic:
    default_resources:
      mem_mb: 64000  # 64GB instead of 32GB
      time: "01:00:00"  # 1 hour instead of 30 min
```

### Monitoring Jobs

#### Check job status
```bash
# For SLURM:
squeue -u $(whoami)

# For PBS:
qstat -u $(whoami)

# For LSF:
bjobs
```

#### Monitor Snakemake progress
While `./submit.sh --profile .smk/profiles/hpc cojo_condition` is running:
- Snakemake prints job submissions to STDOUT
- Check `.snakemake/log/` for detailed Snakemake logs
- Individual job logs stored in `results/direct-cond/log/cojo_condition/<analysis>/<locus>.log.slurm.log`

## Troubleshooting

### Jobs not submitting
- Check cluster module environment: `module avail sbatch` (or qsub, bsub)
- Verify `sbatch` in PATH: `which sbatch`
- Check HPC documentation for cluster submission limits

### Memory errors
- Increase `mem_mb` in `.smk/profiles/hpc/config.yaml`
- Check `sinfo -O memory` (SLURM) to see node memory availability
- Some clusters require node exclusivity for large memory jobs

### Wall-time exceeded
- Increase `time` in profile or `analysis.yml`
- Check GCTA documentation for expected COJO runtime on your data

### Jobs stuck in queue
- May be normal if queue is busy
- Check node availability: `sinfo` (SLURM)
- Try reducing `jobs` to allow faster scheduling of smaller batches

## Running Without Cluster (Local Mode)

To run on local machine again (e.g., for testing):
```bash
./submit.sh --snakefile Snakefile --cores 8 cojo_condition
```

## Example: 25 SNPs workflow
```bash
# Extract regions around each of 25 SNPs (target_snps.txt)
./submit.sh --snakefile Snakefile --cores 2 extract_regions_all
# → Creates 25 folders under results/direct-cond/01_extract_regions/Hispanic/

# Run COJO on all 25 SNPs in parallel on cluster
./submit.sh --snakefile Snakefile --profile .smk/profiles/hpc cojo_condition
# → Submits 25 jobs (one per SNP) to HPC
# → Each job: 32GB memory, 30 min wall time
# → All run in parallel (queue permitting)

# Generate plots (local)
./submit.sh --snakefile Snakefile --cores 2 all
./submit.sh --snakefile Snakefile --cores 2 all_simple_plots
```

## Further Customization

For advanced Snakemake cluster configurations, see:
- [Snakemake Cluster Submission](https://snakemake.readthedocs.io/en/stable/executing/cluster.html)
- [JSON Job Properties](https://snakemake.readthedocs.io/en/stable/executing/cluster.html#job-properties)

To create custom job scripts with more complex logic, add `.smk/profiles/hpc/jobscript.sh` and reference in config.yaml:
```yaml
cluster: "sbatch --job-name={rule}-{wildcards} -J {resources.mem_mb}M ..."
jobscript: ".smk/profiles/hpc/jobscript.sh"
```
