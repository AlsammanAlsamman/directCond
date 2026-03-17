# directCond

Simple Snakemake workflow for regional GWAS conditional analysis (COJO) and plotting.

## What it does

- Step 01: Extract region-specific GWAS records and reference panel subset
- Step 02: Harmonize GWAS alleles with reference panel
- Step 03: Run iterative GCTA-COJO conditioning
- Step 04: Plot before/after conditioning and joint horizontal figure

## Project structure

- `configs/analysis.yml` — project and analysis inputs
- `configs/software.yml` — module/tool configuration
- `rules/*.smk` — Snakemake rules for each step
- `scripts/*` — helper scripts (Python, Bash, R)
- `results/<project_name>/` — all generated outputs

## Requirements

- Snakemake
- PLINK2
- GCTA
- R (with modules from `configs/software.yml`)

## Configuration

1. Edit `configs/analysis.yml`
   - set `project_name`
   - define `target_analyses`
   - point to GWAS and region input files
2. Edit `configs/software.yml`
   - set modules/binaries for `plink2`, `gcta`, and `r`
   - set `r.params.r_libs_user` if needed

## Run

From the project root:

```bash
snakemake --snakefile Snakefile --cores 2
```

Run specific steps:

```bash
# Step 01
snakemake --snakefile Snakefile --cores 2 -R extract_regions

# Step 02
snakemake --snakefile Snakefile --cores 2 -R harmonize_gwas

# Step 03
snakemake --snakefile Snakefile --cores 2 -R cojo_condition

# Step 04
snakemake --snakefile Snakefile --cores 2 -R plot_conditional_iterations
```

## Key outputs

For each analysis/locus under `results/<project_name>/`:

- `03_cojo/<analysis>/<locus>/cojo.ma`
- `03_cojo/<analysis>/<locus>/cojo.cma.cojo`
- `04_plots/<analysis>/<locus>/before_conditioning.pdf`
- `04_plots/<analysis>/<locus>/after_conditioning.pdf`
- `04_plots/<analysis>/<locus>/before_conditioning.png`
- `04_plots/<analysis>/<locus>/after_conditioning.png`
- `04_plots/<analysis>/<locus>/conditional_iterations.pdf`
- `04_plots/<analysis>/<locus>/conditional_iterations.png`

## Notes

- If running on HPC, ensure module names in `configs/software.yml` match your environment.
- Iterative COJO now builds one fixed list of SNPs with `p < 5e-8`, conditions on them in significance order, and only removes candidates that lose significance after each conditioning step. It does not add new SNPs to the candidate list.
- Plotting is robust to missing optional R packages and uses fallbacks when needed.
