#!/usr/bin/env bash
set -euo pipefail

# Step 01: extract GWAS and reference subsets around each locus
./submit.sh --snakefile Snakefile --cores 2 results/direct-cond/01_extract_regions/KLF2_ASIAN/KLF2/extract.done

# Step 02: harmonize GWAS alleles and add rsid using reference subset
./submit.sh --snakefile Snakefile --cores 2 results/direct-cond/02_harmonized/KLF2_ASIAN/KLF2/harmonize.done

# Step 03: perform GCTA-COJO conditional analysis
./submit.sh --snakefile Snakefile --cores 2 results/direct-cond/03_cojo/KLF2_ASIAN/KLF2/cojo.done

# Step 04: plot before/after conditioning and all iterations
./submit.sh --snakefile Snakefile --cores 2 results/direct-cond/04_plots/KLF2_ASIAN/KLF2/plot.done
