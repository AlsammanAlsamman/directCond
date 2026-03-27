# Parallel Locus COJO Submission Guide

## Overview

The `scripts/submit_cojo_parallel_loci.sh` script submits each locus/SNP COJO analysis as a **separate independent task** to your HPC system. Each job runs with 32GB memory independently, and all jobs queue simultaneously.

## Usage

### Basic (all analyses)
```bash
bash scripts/submit_cojo_parallel_loci.sh
```

Submits each locus from all configured analyses in `configs/analysis.yml`.

### Single analysis
```bash
bash scripts/submit_cojo_parallel_loci.sh --analysis Hispanic
```

Submits only loci from the Hispanic analysis.

### Dry-run (show what would be submitted)
```bash
bash scripts/submit_cojo_parallel_loci.sh --dry-run
```

Shows the commands that would be submitted without actually submitting them.

### Combine options
```bash
bash scripts/submit_cojo_parallel_loci.sh --analysis Hispanic --dry-run
```

## How It Works

1. **Reads configuration** from `configs/analysis.yml`
2. **Extracts loci** from either:
   - `regions_file` (legacy mode) - reads locus IDs from TSV
   - `snpList_file` (new mode) - each line is a locus/SNP ID
3. **For each locus**, submits an independent job:
   ```bash
   ./submit.sh --snakefile Snakefile --cores 2 results/direct-cond/03_cojo/{analysis}/{locus}/cojo.done
   ```
4. **All jobs queue in parallel** on your HPC system

## Example Workflow

For 25 SNPs in `inputs/target_snps.txt`:

```bash
# Step 1: Extract regions around all SNPs
./submit.sh --snakefile Snakefile --cores 2 extract_regions_all

# Step 2: Submit 25 independent COJO jobs (each gets 32GB memory)
bash scripts/submit_cojo_parallel_loci.sh

# Step 3-4: Generate plots (after all COJO jobs complete)
./submit.sh --snakefile Snakefile --cores 2 all
./submit.sh --snakefile Snakefile --cores 2 all_simple_plots
```

### Monitor progress
```bash
# Check HPC job queue
squeue -u $(whoami)   # SLURM
qstat -u $(whoami)     # PBS
bjobs                  # LSF

# Check Snakemake logs
tail -f results/direct-cond/log/cojo_condition/*/*log
```

## Resource Allocation

- **Memory per job**: 32 GB (configured in `configs/analysis.yml` → `default_resources.mem_mb`)
- **Wall time per job**: 30 min (configured in `configs/analysis.yml` → `default_resources.time`)
- **Parallelism**: Limited only by HPC queue capacity (all jobs submitted immediately)

### Adjust memory/time
Edit `configs/analysis.yml`:
```yaml
default_resources:
  mem_mb: 64000       # 64GB instead of 32GB
  cores: 2
  time: "01:00:00"    # 1 hour instead of 30 min
```

## Automatic vs Manual Waiting

The script **submits all jobs immediately** and returns. To wait for all jobs to complete:

### Option A: HPC queue commands
```bash
# Wait for all your jobs to finish
# SLURM:
until ! squeue -u $(whoami) | grep -q cojo; do sleep 5; done

# PBS:
until ! qstat -u $(whoami) | grep -q cojo; do sleep 5; done
```

### Option B: Snakemake dependency tracking
After all COJO jobs queue, run plotting which depends on cojo.done:
```bash
./submit.sh --snakefile Snakefile --cores 2 all
```
Snakemake will wait for all upstream COJO jobs before starting plots.

### Option C: Complete workflow in one session (recommended)
```bash
# Run full steps.sh which orchestrates everything
bash steps.sh
```

## Troubleshooting

### Few jobs submitted (not all loci)
- Check snpList_file is valid: `wc -l inputs/target_snps.txt`
- Verify regions_file format (if using legacy mode)
- Run with `--dry-run` to see what would be submitted

### Jobs not queuing
- Check HPC scheduler is available: `which sbatch` (or `qsub`, `bsub`)
- Verify HPC module available: `module avail` or check your cluster docs
- Check queue status: `sinfo` (SLURM) or equivalent for your HPC

### Out of memory errors in jobs
- Increase `mem_mb` in `configs/analysis.yml`
- Check available memory on nodes: `sinfo -O memory` (SLURM)

## Advanced: Parallel Batching

To submit in smaller batches (e.g., 5 jobs at a time) instead of all at once:

```bash
# Limit jobs per submission
for snp in $(grep -v '^\s*$' inputs/target_snps.txt | head -5); do
    ./submit.sh --snakefile Snakefile --cores 2 results/direct-cond/03_cojo/Hispanic/$snp/cojo.done
done
# Then run remaining SNPs...
```

Or modify the script to add a `--batch-size` option. Contact maintainer for help.

## See Also

- `steps.sh` - Full workflow orchestration
- `configs/analysis.yml` - Configuration and resource allocation
- `rules/cojo_condition.smk` - Snakemake rules (technical details)
