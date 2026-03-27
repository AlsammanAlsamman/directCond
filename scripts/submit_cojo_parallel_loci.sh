#!/usr/bin/env bash
################################################################################
# submit_cojo_parallel_loci.sh
# 
# Submit each locus COJO analysis as a separate independent HPC job.
# Each job runs with 32GB memory and all jobs queue simultaneously.
#
# Usage:
#   bash scripts/submit_cojo_parallel_loci.sh [--analysis ANALYSIS] [--dry-run]
#
# Examples:
#   bash scripts/submit_cojo_parallel_loci.sh                    # All analyses
#   bash scripts/submit_cojo_parallel_loci.sh --analysis Hispanic # Single analysis
#   bash scripts/submit_cojo_parallel_loci.sh --dry-run          # Show commands, don't submit
################################################################################

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
CONFIG_FILE="${PROJECT_ROOT}/configs/analysis.yml"

# Parse arguments
DRY_RUN=false
TARGET_ANALYSIS=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --analysis)
            TARGET_ANALYSIS="$2"
            shift 2
            ;;
        --dry-run)
            DRY_RUN=true
            shift
            ;;
        *)
            echo "Unknown option: $1" >&2
            exit 1
            ;;
    esac
done

# Helper: Extract loci from regions_file (legacy format)
_loci_from_regions_file() {
    local file="$1"
    if [[ ! -f "$file" ]]; then
        echo "ERROR: regions file not found: $file" >&2
        return 1
    fi
    tail -n +2 "$file" | awk '{print $1}' | sort -u
}

# Helper: Extract loci from snpList_file (new format - each line is a SNP/locus)
_loci_from_snp_list() {
    local file="$1"
    if [[ ! -f "$file" ]]; then
        echo "ERROR: SNP list file not found: $file" >&2
        return 1
    fi
    grep -v '^\s*$' "$file" | sed 's/\r$//' | sort -u
}

# Helper: Extract value from YAML config
_yaml_value() {
    local yaml_file="$1"
    local key_path="$2"
    python3 <<EOF
import yaml
import sys

try:
    with open("$yaml_file", "r") as f:
        data = yaml.safe_load(f)
    
    keys = "$key_path".split(".")
    val = data
    for k in keys:
        if isinstance(val, dict):
            val = val.get(k)
        else:
            val = None
        if val is None:
            sys.exit(1)
    
    if isinstance(val, (list, dict)):
        print(yaml.dump(val, default_flow_style=False))
    else:
        print(str(val))
except Exception as e:
    print(f"Error extracting {key_path}: {e}", file=sys.stderr)
    sys.exit(1)
EOF
}

# Main

echo "=========================================="
echo "COJO Parallel Locus Submission"
echo "=========================================="

cd "$PROJECT_ROOT"

# Get all analyses or use target
if [[ ! -f "$CONFIG_FILE" ]]; then
    echo "ERROR: Config file not found: $CONFIG_FILE" >&2
    exit 1
fi

if [[ -z "$TARGET_ANALYSIS" ]]; then
    echo "Reading all analyses from $CONFIG_FILE ..."
    # Use Python to extract analyses
    ANALYSES=$(python3 <<'PYSCRIPT'
import yaml
with open("configs/analysis.yml") as f:
    data = yaml.safe_load(f)
analyses = list(data.get("target_analyses", {}).keys())
for a in analyses:
    print(a)
PYSCRIPT
)
else
    ANALYSES="$TARGET_ANALYSIS"
fi

TOTAL_JOBS=0

# For each analysis
for ANALYSIS in $ANALYSES; do
    echo ""
    echo ">>> Analysis: $ANALYSIS"
    
    # Check if regions_file or snpList_file exists
    REGIONS_FILE=$(python3 <<EOF
import yaml
with open("configs/analysis.yml") as f:
    data = yaml.safe_load(f)
cfg = data.get("target_analyses", {}).get("$ANALYSIS", {})
print(cfg.get("regions_file", ""))
EOF
) || REGIONS_FILE=""

    SNP_LIST_FILE=$(python3 <<EOF
import yaml
with open("configs/analysis.yml") as f:
    data = yaml.safe_load(f)
cfg = data.get("target_analyses", {}).get("$ANALYSIS", {})
print(cfg.get("snpList_file", ""))
EOF
) || SNP_LIST_FILE=""

    # Get loci list
    if [[ -n "$REGIONS_FILE" && -f "$REGIONS_FILE" ]]; then
        echo "  Source: regions_file = $REGIONS_FILE"
        LOCI=$(_loci_from_regions_file "$REGIONS_FILE")
    elif [[ -n "$SNP_LIST_FILE" && -f "$SNP_LIST_FILE" ]]; then
        echo "  Source: snpList_file = $SNP_LIST_FILE"
        LOCI=$(_loci_from_snp_list "$SNP_LIST_FILE")
    else
        echo "  ERROR: No regions_file or snpList_file found for $ANALYSIS" >&2
        exit 1
    fi

    # Count and submit
    LOCUS_COUNT=$(echo "$LOCI" | wc -l)
    echo "  Loci found: $LOCUS_COUNT"
    echo ""

    for LOCUS in $LOCI; do
        echo "    Submitting: $ANALYSIS / $LOCUS"
        
        # Build Snakemake target (use wildcard expansion to single rule)
        TARGET="results/direct-cond/03_cojo/${ANALYSIS}/${LOCUS}/cojo.done"
        
        # Build command
        CMD="./submit.sh --snakefile Snakefile --cores 2 \"$TARGET\""
        
        if [[ "$DRY_RUN" == "true" ]]; then
            echo "      [DRY RUN] $CMD"
        else
            eval "$CMD"
        fi
        
        ((TOTAL_JOBS++))
    done
done

echo ""
echo "=========================================="
echo "Submitted: $TOTAL_JOBS jobs"
if [[ "$DRY_RUN" == "true" ]]; then
    echo "(Dry-run mode - no jobs actually submitted)"
fi
echo "=========================================="
