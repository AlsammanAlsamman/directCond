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
   - optionally set `snpList_file` to define locus centers for Step 01 when you do not provide a `regions_file`
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

# Excel report of conditioned SNPs per locus
snakemake --snakefile Snakefile --cores 1 conditioned_snp_excel_report
```

## Key outputs

For each analysis/locus under `results/<project_name>/`:

- `03_cojo/<analysis>/<locus>/cojo.ma`
- `03_cojo/<analysis>/<locus>/cojo.cma.cojo`
- `05_reports/<analysis>/conditioned_snps.xlsx`
- `04_plots/<analysis>/<locus>/before_conditioning.pdf`
- `04_plots/<analysis>/<locus>/after_conditioning.pdf`
- `04_plots/<analysis>/<locus>/before_conditioning.png`
- `04_plots/<analysis>/<locus>/after_conditioning.png`
- `04_plots/<analysis>/<locus>/conditional_iterations.pdf`
- `04_plots/<analysis>/<locus>/conditional_iterations.png`

## Notes

- If running on HPC, ensure module names in `configs/software.yml` match your environment.
- The folder name under `results/<project_name>/.../<locus>/` is a locus label inherited from the region-definition input. It is not treated as a required conditioning SNP.
- Step 03 always runs locus-local iterative COJO using the harmonized GWAS and the matching locus PLINK subset. Conditioning stops when no matched SNP in that locus remains significant.
- If no matched SNP in a locus passes the significance threshold, Step 03 now writes an unconditioned final COJO table for that locus so downstream plotting can still run.
- The Excel report rule reads the original GWAS file from `target_analyses.<name>.gwas_file`, joins it to each locus `cojo.cond.snp` list and `cojo.cma.cojo`, and writes a styled workbook per analysis.
- Plotting is robust to missing optional R packages and uses fallbacks when needed.
