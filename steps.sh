#!/usr/bin/env bash
set -euo pipefail

# Step 01: extract GWAS and reference subsets around each SNP in snpList_file
./submit.sh --cores 2 extract_regions_all

# Step 02: submit each locus COJO as separate HPC job
# --jobs 25: submit up to 25 parallel jobs (each locus gets 32GB memory)
# All loci run in parallel on HPC cluster
./submit.sh --jobs 25 cojo_condition

# Step 03: generate full before/after iteration plots
./submit.sh --cores 2 all

# Step 04: generate simple plots (original, each single-SNP conditioning, then all SNPs)
./submit.sh --cores 2 all_simple_plots

# Step 05: generate Excel report of conditioned SNPs per locus from the original GWAS file
./submit.sh --jobs 1 conditioned_snp_excel_report