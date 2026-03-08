#!/usr/bin/env bash
set -euo pipefail

snakemake --snakefile Snakefile --cores 2 "$@"
