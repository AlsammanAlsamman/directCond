#!/usr/bin/env bash
# run_indiv_cond_eval.sh
#
# For each SNP in cojo.cond.snp:
#   1. Run GCTA --cojo-cond on that single SNP using the existing cojo.ma
#   2. Count how many SNPs in the resulting .cma.cojo have pC < 5e-8
#   3. Look up the SNP's own pC from the final multi-SNP cojo.cma.cojo
#   4. Write one summary row: snp, original_p, own_pC, n_significant_after_cond, status
#
# Usage:
#   bash scripts/run_indiv_cond_eval.sh \
#     --cond-snp    <cojo.cond.snp>        \
#     --eval-snp-list <list_of_snps_to_evaluate> \
#     --cojo-ma     <cojo.ma>              \
#     --cojo-final  <cojo.cma.cojo>        \
#     --bfile-prefix <ref_subset_prefix>   \
#     --out-dir     <working_dir>          \
#     --out-summary <summary.tsv>          \
#     --out-done    <done_marker>          \
#     [--gcta-bin   <gcta>]

set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  bash scripts/run_indiv_cond_eval.sh \
    --cond-snp    <cojo.cond.snp>      \
    --eval-snp-list <list_of_snps_to_evaluate> \
    --cojo-ma     <cojo.ma>            \
    --cojo-final  <cojo.cma.cojo>      \
    --bfile-prefix <ref_subset_prefix> \
    --out-dir     <working_dir>        \
    --out-summary <summary.tsv>        \
    --out-done    <done_marker>        \
    [--gcta-bin   <gcta>]
EOF
}

COND_SNP_FILE=""
EVAL_SNP_LIST=""
COJO_MA=""
COJO_FINAL=""
BFILE_PREFIX=""
OUT_DIR=""
OUT_SUMMARY=""
OUT_DONE=""
GCTA_BIN="gcta"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --cond-snp)    COND_SNP_FILE="$2"; shift 2 ;;
    --eval-snp-list) EVAL_SNP_LIST="$2"; shift 2 ;;
    --cojo-ma)     COJO_MA="$2";       shift 2 ;;
    --cojo-final)  COJO_FINAL="$2";    shift 2 ;;
    --bfile-prefix) BFILE_PREFIX="$2"; shift 2 ;;
    --out-dir)     OUT_DIR="$2";       shift 2 ;;
    --out-summary) OUT_SUMMARY="$2";   shift 2 ;;
    --out-done)    OUT_DONE="$2";      shift 2 ;;
    --gcta-bin)    GCTA_BIN="$2";      shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage; exit 1 ;;
  esac
done

echo "[DEBUG] Parsed arguments:" >&2
echo "  COND_SNP_FILE = ${COND_SNP_FILE}" >&2
echo "  EVAL_SNP_LIST = ${EVAL_SNP_LIST}" >&2
echo "  COJO_MA       = ${COJO_MA}" >&2
echo "  COJO_FINAL    = ${COJO_FINAL}" >&2
echo "  BFILE_PREFIX  = ${BFILE_PREFIX}" >&2
echo "  OUT_DIR       = ${OUT_DIR}" >&2
echo "  OUT_SUMMARY   = ${OUT_SUMMARY}" >&2
echo "  OUT_DONE      = ${OUT_DONE}" >&2
echo "  GCTA_BIN      = ${GCTA_BIN}" >&2

for var_name in COND_SNP_FILE COJO_MA COJO_FINAL BFILE_PREFIX OUT_DIR OUT_SUMMARY OUT_DONE; do
  if [[ -z "${!var_name}" ]]; then
    echo "Missing required argument: --${var_name//_/-}" >&2
    usage
    exit 1
  fi
done

if [[ -z "$EVAL_SNP_LIST" ]]; then
  EVAL_SNP_LIST="$COND_SNP_FILE"
fi

echo "[DEBUG] Checking required input files exist..." >&2
for f in "$COND_SNP_FILE" "$EVAL_SNP_LIST" "$COJO_MA" "$COJO_FINAL"; do
  if [[ ! -f "$f" ]]; then
    echo "[DEBUG] FAIL - File not found: $f" >&2
    exit 1
  else
    echo "[DEBUG]   OK: $f" >&2
  fi
done

echo "[DEBUG] Checking reference panel files..." >&2
for ext in bed bim fam; do
  if [[ ! -f "${BFILE_PREFIX}.${ext}" ]]; then
    echo "[DEBUG] FAIL - Missing reference panel file: ${BFILE_PREFIX}.${ext}" >&2
    exit 1
  else
    echo "[DEBUG]   OK: ${BFILE_PREFIX}.${ext}" >&2
  fi
done

echo "[DEBUG] Checking GCTA binary..." >&2
if ! command -v "$GCTA_BIN" >/dev/null 2>&1; then
  echo "[DEBUG] FAIL - GCTA executable not found: $GCTA_BIN" >&2
  exit 1
fi
echo "[DEBUG]   OK: GCTA found at $(command -v "$GCTA_BIN")" >&2

echo "[DEBUG] Creating output directories..." >&2
mkdir -p "$OUT_DIR" "$(dirname "$OUT_DONE")" "$(dirname "$OUT_SUMMARY")"
echo "[DEBUG]   OUT_DIR created: $OUT_DIR" >&2

P_THRESHOLD="5e-8"
echo "[DEBUG] P_THRESHOLD = $P_THRESHOLD" >&2

# ── build lookup: snp -> original_p from cojo.ma ────────────────────────────
build_p_lookup() {
  local ma_file="$1"
  awk '
BEGIN { FS="\t" }
NR==1 { next }
{ if ($1 != "" && $7 != "") print $1"\t"$7 }
' "$ma_file"
}

# ── build lookup: snp -> pC from the multi-SNP final cojo.cma.cojo ──────────
build_pc_lookup() {
  local cma_file="$1"
  awk '
BEGIN { FS="\t" }
NR==1 {
  for (i=1; i<=NF; i++) {
    h=tolower($i); gsub(/^ +| +$/, "", h)
    if (h=="snp") c_snp=i
    else if (h=="pc") c_pc=i
  }
  next
}
{
  if (!c_snp || !c_pc) next
  snp=$c_snp; pc=$c_pc
  if (snp != "" && pc != "") print snp"\t"pc
}
' "$cma_file"
}

echo "[DEBUG] Building P lookup from COJO_MA..." >&2
P_LOOKUP=$(build_p_lookup "$COJO_MA")
echo "[DEBUG]   P_LOOKUP entries: $(echo "$P_LOOKUP" | wc -l)" >&2

echo "[DEBUG] Building PC lookup from COJO_FINAL..." >&2
PC_LOOKUP=$(build_pc_lookup "$COJO_FINAL")
echo "[DEBUG]   PC_LOOKUP entries: $(echo "$PC_LOOKUP" | wc -l)" >&2

# ── write summary header ─────────────────────────────────────────────────────
echo "[DEBUG] Writing summary header to: $OUT_SUMMARY" >&2
printf "snp\toriginal_p\town_pC_in_final\tn_significant_after_cond\tstatus\n" > "$OUT_SUMMARY"

# ── read GCTA chromosome from BIM (col 1 of first data line) ─────────────────
BIM_CHR=$(awk 'NR==1 {print $1; exit}' "${BFILE_PREFIX}.bim")
echo "[DEBUG] BIM_CHR = $BIM_CHR" >&2
echo "[DEBUG] Number of conditioning SNPs: $(wc -l < "$COND_SNP_FILE")" >&2
echo "[DEBUG] Contents of cond.snp file:" >&2
cat "$COND_SNP_FILE" >&2
echo "" >&2

# ── loop over each conditioned SNP ───────────────────────────────────────────
while IFS= read -r SNP || [[ -n "$SNP" ]]; do
  SNP=$(echo "$SNP" | tr -d '[:space:]')
  [[ -z "$SNP" ]] && continue

  # Use here-strings instead of echo|awk to avoid pipefail+SIGPIPE exits.
  ORIG_P=$(awk -v s="$SNP" '$1==s {print $2; exit}' <<< "$P_LOOKUP")
  OWN_PC=$(awk -v s="$SNP" '$1==s {print $2; exit}' <<< "$PC_LOOKUP")
  [[ -z "$ORIG_P" ]] && ORIG_P="NA"
  [[ -z "$OWN_PC" ]] && OWN_PC="NA"

  echo "[DEBUG] ── Processing SNP: ${SNP} ──" >&2
  echo "[DEBUG]   ORIG_P = ${ORIG_P}" >&2
  echo "[DEBUG]   OWN_PC = ${OWN_PC}" >&2

  ITER_DIR="${OUT_DIR}/indiv_${SNP}"
  mkdir -p "$ITER_DIR"
  ITER_PREFIX="${ITER_DIR}/cojo"
  SINGLE_COND_FILE="${ITER_DIR}/single.cond.snp"
  GCTA_LOG="${ITER_DIR}/gcta.log"

  printf "%s\n" "$SNP" > "$SINGLE_COND_FILE"
  echo "[DEBUG]   Single cond file written: ${SINGLE_COND_FILE}" >&2

  echo "Running individual conditioning on SNP: ${SNP}" >&2

  STATUS="ok"
  N_SIG="NA"

  # run GCTA with this single SNP as the condition
  echo "[DEBUG]   Running: $GCTA_BIN --bfile $BFILE_PREFIX --chr $BIM_CHR --cojo-file $COJO_MA --cojo-cond $SINGLE_COND_FILE --out $ITER_PREFIX" >&2
  if ! "$GCTA_BIN" \
      --bfile "$BFILE_PREFIX" \
      --chr   "$BIM_CHR" \
      --cojo-file "$COJO_MA" \
      --cojo-cond "$SINGLE_COND_FILE" \
      --out "$ITER_PREFIX" > "$GCTA_LOG" 2>&1; then

    if grep -qi "collinearity problem" "$GCTA_LOG"; then
      STATUS="collinear"
      echo "[DEBUG]   SNP ${SNP}: collinearity — skipping count." >&2
    else
      STATUS="gcta_failed"
      echo "[DEBUG]   SNP ${SNP}: GCTA failed — skipping count." >&2
      echo "[DEBUG]   GCTA log contents:" >&2
      cat "$GCTA_LOG" >&2
    fi
  else
    echo "[DEBUG]   GCTA completed successfully for SNP: ${SNP}" >&2
  fi

  if [[ "$STATUS" == "ok" ]]; then
    echo "[DEBUG]   Checking for output file: ${ITER_PREFIX}.cma.cojo" >&2
    if [[ -f "${ITER_PREFIX}.cma.cojo" ]]; then
      echo "[DEBUG]   Output file found, counting significant SNPs only in eval list..." >&2
      N_SIG=$(awk -v thr="$P_THRESHOLD" '
BEGIN { FS="\t"; n=0 }
FNR==NR {
  s=$1
  gsub(/^[[:space:]]+|[[:space:]]+$/, "", s)
  if (s!="") keep[s]=1
  next
}
FNR==1 {
  for (i=1; i<=NF; i++) {
    h=tolower($i); gsub(/^ +| +$/, "", h)
    if (h=="pc") c_pc=i
    else if (h=="snp") c_snp=i
  }
  next
}
{
  if (!c_pc || !c_snp) next
  snp=$c_snp
  pc=$c_pc
  if (!(snp in keep)) next
  if (pc=="" || toupper(pc)=="NA") next
  if ((pc+0) < thr) n++
}
END { print n }
' "$EVAL_SNP_LIST" "${ITER_PREFIX}.cma.cojo")
      echo "[DEBUG]   N_SIG = ${N_SIG} for SNP: ${SNP}" >&2
    else
      STATUS="missing_output"
      echo "[DEBUG]   FAIL - SNP ${SNP}: expected COJO output missing: ${ITER_PREFIX}.cma.cojo" >&2
      echo "[DEBUG]   Files in ITER_DIR:" >&2
      ls -la "$ITER_DIR" >&2 || true
    fi
  fi

  echo "[DEBUG]   Writing summary row: SNP=${SNP} STATUS=${STATUS}" >&2
  printf "%s\t%s\t%s\t%s\t%s\n" "$SNP" "$ORIG_P" "$OWN_PC" "$N_SIG" "$STATUS" >> "$OUT_SUMMARY"

  # clean up intermediate GCTA files to save disk space
  echo "[DEBUG]   Cleaning up intermediate files for SNP: ${SNP}" >&2
  rm -f "${ITER_PREFIX}.cma.cojo" "${ITER_PREFIX}.jma.cojo" \
        "${ITER_PREFIX}.ldr.cojo" "${ITER_PREFIX}.given.cojo"
  echo "[DEBUG]   Done with SNP: ${SNP}" >&2

done < "$COND_SNP_FILE"

echo "" >&2
echo "[DEBUG] All SNPs processed. Writing final outputs..." >&2
echo "Individual conditioning summary written to: ${OUT_SUMMARY}" >&2
echo "SNPs evaluated: $(tail -n +2 "$OUT_SUMMARY" | wc -l)" >&2
echo "SNPs with n_significant_after_cond < threshold available in summary." >&2

echo "[DEBUG] Writing done marker: ${OUT_DONE}" >&2
printf "ok\tsummary=%s\tn_snps_evaluated=%s\n" \
  "$OUT_SUMMARY" \
  "$(tail -n +2 "$OUT_SUMMARY" | wc -l)" > "$OUT_DONE"
echo "[DEBUG] Script completed successfully." >&2
