#!/usr/bin/env bash
################################################################################
# submit.sh
# Wrapper for Snakemake with HPC support
#
# Usage:
#   ./submit.sh --cores 2 <target>              # Local execution with 2 cores
#   ./submit.sh --jobs 25 <target>              # HPC submission: 25 jobs in parallel
#                                               # (each locus = separate HPC task, 32GB each)
#
# Examples:
#   ./submit.sh --cores 2 extract_regions_all
#   ./submit.sh --cores 2 cojo_condition
#   ./submit.sh --jobs 25 cojo_condition        # Submit each locus to HPC
#   ./submit.sh --jobs 25 --snakefile Snakefile cojo_condition
################################################################################

set -euo pipefail

# Check if --jobs flag is present
HAS_JOBS=false

while [[ $# -gt 0 ]]; do
    case "$1" in
        --jobs)
            if [[ $# -lt 2 ]]; then
                echo "Error: --jobs requires a value." >&2
                exit 1
            fi
            HAS_JOBS=true
            break
            ;;
        --jobs=*)
            HAS_JOBS=true
            break
            ;;
    esac
    shift
done

# If --jobs is present, use HPC cluster profile
if [[ "$HAS_JOBS" == "true" ]]; then
    # Apply cluster profile for parallel HPC submission
    snakemake --profile .smk/profiles/hpc "$@"
else
    # Local execution (default behavior)
    snakemake "$@"
fi
