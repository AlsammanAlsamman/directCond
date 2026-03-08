#!/usr/bin/env bash
set -euo pipefail

snakemake --snakefile rules/cojo_condition.smk --cores 2 "$@"
