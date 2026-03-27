#!/usr/bin/env bash
set -euo pipefail

# Step 01: extract GWAS and reference subsets around each SNP in snpList_file
./submit.sh --jobs 25 --rerun-incomplete extract_regions_all

# Step 02: run locus COJO jobs on HPC
./submit.sh --jobs 25 cojo_condition

# Step 03: generate full before/after iteration plots
./submit.sh --jobs 10 all

# Step 04: generate simple plots (original, each single-SNP conditioning, then all SNPs)
./submit.sh --jobs 10 all_simple_plots

# Step 05: generate Excel report of conditioned SNPs per locus from the original GWAS file
./submit.sh --force --jobs 1 conditioned_snp_excel_report