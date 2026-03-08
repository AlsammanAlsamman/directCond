#!/usr/bin/env bash
set -euo pipefail

snakemake --snakefile rules/plot_conditional_iterations.smk --cores 2 "$@"
