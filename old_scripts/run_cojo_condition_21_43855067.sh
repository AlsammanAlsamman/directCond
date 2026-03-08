#!/usr/bin/env bash
set -euo pipefail

# Conditional analysis in COJO for chr21:43855067
#
# Usage:
#   bash run_cojo_condition_21_43855067.sh [SUB_GWAS] [BFILE_PREFIX] [N] [OUT_PREFIX] [GCTA_BIN]
#
# Defaults:
#   SUB_GWAS     = chr21_43855067_pm4e+05bp_window.tsv
#   BFILE_PREFIX = ref/g1000_eas (or set arg2 / BFILE_PREFIX env var)
#   N            = 100000
#   OUT_PREFIX   = cojo_chr21_43855067
#   GCTA_BIN     = /usr/local/analysis/gcta/1.26.0/bin/gcta

SUB_GWAS="${1:-chr21_43855067_pm4e+05bp_window.tsv}"
BFILE_PREFIX="${2:-${BFILE_PREFIX:-ref/g1000_eas}}"
N="${3:-100000}"
OUT_PREFIX="${4:-cojo_chr21_43855067}"
GCTA_BIN="${5:-/usr/local/analysis/gcta/1.26.0/bin/gcta}"

TARGET_CHR="21"
TARGET_POS="43855067"
WINDOW_BP="400000"
START_POS="$((TARGET_POS - WINDOW_BP))"
END_POS="$((TARGET_POS + WINDOW_BP))"

if [[ ! -f "$SUB_GWAS" ]]; then
  echo "ERROR: GWAS file not found: $SUB_GWAS" >&2
  exit 1
fi

find_bfile_prefix() {
  local bed_file
  while IFS= read -r bed_file; do
    local prefix="${bed_file%.bed}"
    if [[ -f "${prefix}.bim" && -f "${prefix}.fam" ]]; then
      echo "$prefix"
      return 0
    fi
  done < <(find . -maxdepth 4 -type f -name "*.bed" | sort)
  return 1
}

if [[ -z "$BFILE_PREFIX" ]]; then
  if DETECTED_PREFIX="$(find_bfile_prefix)"; then
    BFILE_PREFIX="$DETECTED_PREFIX"
    echo "Auto-detected BFILE prefix: $BFILE_PREFIX"
  else
    echo "ERROR: Missing BFILE prefix and no local *.bed/*.bim/*.fam set was found." >&2
    echo "Provide as arg2 or env var BFILE_PREFIX." >&2
    echo "Example:" >&2
    echo "  bash $0 $SUB_GWAS /path/to/ld_ref_prefix 200000" >&2
    echo "Expected files:" >&2
    echo "  /path/to/ld_ref_prefix.bed" >&2
    echo "  /path/to/ld_ref_prefix.bim" >&2
    echo "  /path/to/ld_ref_prefix.fam" >&2
    exit 1
  fi
fi

for ext in bed bim fam; do
  if [[ ! -f "${BFILE_PREFIX}.${ext}" ]]; then
    echo "ERROR: Missing ${BFILE_PREFIX}.${ext}" >&2
    exit 1
  fi
done

if ! command -v "$GCTA_BIN" >/dev/null 2>&1; then
  echo "ERROR: GCTA binary not found in PATH: $GCTA_BIN" >&2
  exit 1
fi

MA_FILE="${OUT_PREFIX}.ma"
COND_FILE="${OUT_PREFIX}.cond.snp"
LOG_PREFIX="${OUT_PREFIX}"

echo "Preparing COJO summary file: $MA_FILE"

echo -e "SNP\tA1\tA2\tfreq\tb\tse\tp\tN" > "$MA_FILE"

awk -v chr="$TARGET_CHR" -v start="$START_POS" -v end="$END_POS" -v N="$N" '
BEGIN { FS="\t"; OFS="\t" }
NR==1 {
  for (i=1; i<=NF; i++) {
    h=tolower($i)
    gsub(/^ +| +$/, "", h)
    if (h=="chrom" || h=="chr") c_chrom=i
    else if (h=="pos" || h=="position" || h=="bp") c_pos=i
    else if (h=="ea" || h=="a1" || h=="effect_allele") c_a1=i
    else if (h=="nea" || h=="a2" || h=="other_allele" || h=="non_effect_allele") c_a2=i
    else if (h=="beta" || h=="b" || h=="effect") c_b=i
    else if (h=="se") c_se=i
    else if (h=="p" || h=="pval" || h=="pvalue") c_p=i
    else if (h=="varid" || h=="snp" || h=="rsid") c_snp=i
    else if (h=="freq" || h=="eaf" || h=="af" || h=="maf") c_freq=i
  }
  if (!c_chrom || !c_pos || !c_a1 || !c_a2 || !c_b || !c_se || !c_p) {
    print "ERROR: Required columns not found in header." > "/dev/stderr"
    exit 2
  }
  next
}
{
  row_chr=$c_chrom
  gsub(/^chr/i, "", row_chr)
  row_pos=$c_pos + 0
  if (row_chr != chr) next
  if (row_pos < start || row_pos > end) next

  snp=(c_snp ? $c_snp : row_chr ":" row_pos ":" $c_a2 ":" $c_a1)
  a1=$c_a1
  a2=$c_a2
  b=$c_b
  se=$c_se
  p=$c_p
  freq=(c_freq ? $c_freq : 0.5)

  if (snp=="" || a1=="" || a2=="" || b=="" || se=="" || p=="") next
  print snp, a1, a2, freq, b, se, p, N
}
' "$SUB_GWAS" >> "$MA_FILE"

if [[ $(wc -l < "$MA_FILE") -le 1 ]]; then
  echo "ERROR: No rows found for chr${TARGET_CHR}:${START_POS}-${END_POS} in $SUB_GWAS" >&2
  exit 1
fi

COND_SNP=$(awk -v target="${TARGET_CHR}:${TARGET_POS}:" '
BEGIN{FS="\t"}
NR==1{next}
$1 ~ "^"target {
  if (best=="" || ($7+0) < bestp) { best=$1; bestp=($7+0) }
}
END{print best}
' "$MA_FILE")

if [[ -z "$COND_SNP" ]]; then
  COND_SNP=$(awk -v target="${TARGET_CHR}:${TARGET_POS}$" 'BEGIN{FS="\t"} NR>1 && $1 ~ target {print $1; exit}' "$MA_FILE")
fi

if [[ -z "$COND_SNP" ]]; then
  echo "ERROR: Could not find SNP at ${TARGET_CHR}:${TARGET_POS} in $MA_FILE" >&2
  echo "Check variant IDs in your input (expected varid like 21:43855067:A:G)." >&2
  exit 1
fi

echo "$COND_SNP" > "$COND_FILE"

echo "Conditioning SNP: $COND_SNP"
echo "Running COJO..."

"$GCTA_BIN" \
  --bfile "$BFILE_PREFIX" \
  --chr "$TARGET_CHR" \
  --cojo-file "$MA_FILE" \
  --cojo-cond "$COND_SNP" \
  --out "$LOG_PREFIX"

echo "Done. Outputs prefix: $LOG_PREFIX"
