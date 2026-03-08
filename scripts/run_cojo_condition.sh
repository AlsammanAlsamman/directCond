#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  bash scripts/run_cojo_condition.sh \
    --gwas-file <gwas_window.tsv> \
    --bfile-prefix <ref_subset_prefix> \
    --target-chr <chr> \
    --target-pos <pos> \
    --out-prefix <output_prefix_without_extension> \
    --out-done <done_file> \
    [--gcta-bin <gcta_binary>]

Required arguments:
  --gwas-file    GWAS window TSV file (must contain chr/pos/allele/beta/se/p columns)
  --bfile-prefix PLINK prefix with .bed/.bim/.fam
  --target-chr   Target chromosome for condition SNP
  --target-pos   Target position for condition SNP
  --out-prefix   Output prefix for COJO outputs (e.g., .../cojo)
  --out-done     Done marker file path

Optional arguments:
  --gcta-bin     GCTA executable (default: gcta)
EOF
}

GWAS_FILE=""
BFILE_PREFIX=""
TARGET_CHR=""
TARGET_POS=""
OUT_PREFIX=""
OUT_DONE=""
GCTA_BIN="gcta"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --gwas-file) GWAS_FILE="$2"; shift 2 ;;
    --bfile-prefix) BFILE_PREFIX="$2"; shift 2 ;;
    --target-chr) TARGET_CHR="$2"; shift 2 ;;
    --target-pos) TARGET_POS="$2"; shift 2 ;;
    --out-prefix) OUT_PREFIX="$2"; shift 2 ;;
    --out-done) OUT_DONE="$2"; shift 2 ;;
    --gcta-bin) GCTA_BIN="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *)
      echo "Unknown argument: $1" >&2
      usage
      exit 1
      ;;
  esac
done

if [[ -z "$GWAS_FILE" || -z "$BFILE_PREFIX" || -z "$TARGET_CHR" || -z "$TARGET_POS" || -z "$OUT_PREFIX" || -z "$OUT_DONE" ]]; then
  echo "Missing required arguments." >&2
  usage
  exit 1
fi

if [[ ! -f "$GWAS_FILE" ]]; then
  echo "GWAS file not found: $GWAS_FILE" >&2
  exit 1
fi

for ext in bed bim fam; do
  if [[ ! -f "${BFILE_PREFIX}.${ext}" ]]; then
    echo "Missing reference panel file: ${BFILE_PREFIX}.${ext}" >&2
    exit 1
  fi
done

if ! command -v "$GCTA_BIN" >/dev/null 2>&1; then
  echo "GCTA executable not found: $GCTA_BIN" >&2
  exit 1
fi

OUT_DIR="$(dirname "$OUT_PREFIX")"
mkdir -p "$OUT_DIR" "$(dirname "$OUT_DONE")"

MA_FILE="${OUT_PREFIX}.ma"
COND_FILE="${OUT_PREFIX}.cond.snp"

printf "SNP\tA1\tA2\tfreq\tb\tse\tp\tN\n" > "$MA_FILE"

awk '
BEGIN { FS="\t"; OFS="\t" }
NR==1 {
  for (i=1; i<=NF; i++) {
    h=tolower($i)
    gsub(/^ +| +$/, "", h)
    if (h=="chrom" || h=="chr") c_chr=i
    else if (h=="pos" || h=="position" || h=="bp") c_pos=i
    else if (h=="ea" || h=="a1" || h=="effect_allele" || h=="effectallele") c_a1=i
    else if (h=="nea" || h=="a2" || h=="other_allele" || h=="otherallele" || h=="non_effect_allele") c_a2=i
    else if (h=="beta" || h=="b" || h=="effect") c_b=i
    else if (h=="se") c_se=i
    else if (h=="p" || h=="pvalue" || h=="pval") c_p=i
    else if (h=="rsid") c_rsid=i
    else if (h=="varid" || h=="snp") c_snp=i
    else if (h=="freq" || h=="eaf" || h=="af" || h=="maf") c_freq=i
    else if (h=="n" || h=="samples") c_n=i
    else if (h=="harmonization_status") c_hstatus=i
  }
  if (!c_chr || !c_pos || !c_a1 || !c_a2 || !c_b || !c_se || !c_p) {
    print "Missing required columns in GWAS file header." > "/dev/stderr"
    exit 2
  }
  next
}
{
  chr=$c_chr
  gsub(/^chr/i, "", chr)
  pos=$c_pos + 0
  a1=$c_a1
  a2=$c_a2
  b=$c_b
  se=$c_se
  p=$c_p
  if (c_hstatus && $c_hstatus=="no_ref_match") next

  snp=""
  if (c_rsid && $c_rsid!="") {
    snp=$c_rsid
  } else if (c_snp && $c_snp!="") {
    snp=$c_snp
  } else {
    snp=chr ":" pos ":" a2 ":" a1
  }
  freq=(c_freq ? $c_freq : 0.5)
  n=(c_n ? $c_n : 100000)

  if (snp=="" || a1=="" || a2=="" || b=="" || se=="" || p=="") next
  print snp, a1, a2, freq, b, se, p, n
}
' "$GWAS_FILE" >> "$MA_FILE"

if [[ $(wc -l < "$MA_FILE") -le 1 ]]; then
  echo "No valid records written to MA file: $MA_FILE" >&2
  exit 1
fi

FIRST_SNP=$(awk '
BEGIN { FS="\t" }
NR==1 { next }
{
  snp=$1
  p=$7 + 0
  if (snp=="") next
  if (best=="" || p < bestp) {
    best=snp
    bestp=p
  }
}
END { print best }
' "$MA_FILE")

if [[ -z "$FIRST_SNP" ]]; then
  echo "Could not determine initial conditioning SNP from $MA_FILE" >&2
  exit 1
fi

printf "%s\n" "$FIRST_SNP" > "$COND_FILE"
echo "Initial conditioning SNP: $FIRST_SNP" >&2

P_THRESHOLD="5e-8"
MAX_ITERS=50
ITER=1
FINAL_PREFIX=""

while [[ $ITER -le $MAX_ITERS ]]; do
  ITER_PREFIX="${OUT_PREFIX}.iter${ITER}"

  "$GCTA_BIN" \
    --bfile "$BFILE_PREFIX" \
    --chr "$TARGET_CHR" \
    --cojo-file "$MA_FILE" \
    --cojo-cond "$COND_FILE" \
    --out "$ITER_PREFIX"

  if [[ ! -f "${ITER_PREFIX}.cma.cojo" ]]; then
    echo "Expected COJO output missing: ${ITER_PREFIX}.cma.cojo" >&2
    exit 1
  fi

  FINAL_PREFIX="$ITER_PREFIX"

  NEXT_LINE=$(awk -v thr="$P_THRESHOLD" '
BEGIN { FS="\t"; OFS="\t" }
NR==1 {
  for (i=1; i<=NF; i++) {
    h=tolower($i)
    gsub(/^ +| +$/, "", h)
    if (h=="snp") c_snp=i
    else if (h=="pc") c_pc=i
  }
  next
}
{
  if (!c_snp || !c_pc) next
  snp=$c_snp
  pc=$c_pc + 0
  if (snp=="") next
  if (pc < thr) {
    if (best=="" || pc < best_pc) {
      best=snp
      best_pc=pc
    }
  }
}
END {
  if (best!="") print best, best_pc
}
' "${ITER_PREFIX}.cma.cojo")

  if [[ -z "$NEXT_LINE" ]]; then
    echo "Stopping iterative conditioning at iteration ${ITER}: no SNP with pC < ${P_THRESHOLD}." >&2
    break
  fi

  NEXT_SNP=$(echo "$NEXT_LINE" | awk '{print $1}')
  NEXT_PC=$(echo "$NEXT_LINE" | awk '{print $2}')

  if grep -Fxq "$NEXT_SNP" "$COND_FILE"; then
    echo "Stopping iterative conditioning at iteration ${ITER}: next SNP ${NEXT_SNP} already conditioned." >&2
    break
  fi

  echo "$NEXT_SNP" >> "$COND_FILE"
  echo "Iteration ${ITER}: adding SNP ${NEXT_SNP} with pC=${NEXT_PC}" >&2
  ITER=$((ITER + 1))
done

if [[ -z "$FINAL_PREFIX" ]]; then
  echo "COJO did not produce any iteration output." >&2
  exit 1
fi

cp -f "${FINAL_PREFIX}.cma.cojo" "${OUT_PREFIX}.cma.cojo"
[[ -f "${FINAL_PREFIX}.jma.cojo" ]] && cp -f "${FINAL_PREFIX}.jma.cojo" "${OUT_PREFIX}.jma.cojo"
[[ -f "${FINAL_PREFIX}.ldr.cojo" ]] && cp -f "${FINAL_PREFIX}.ldr.cojo" "${OUT_PREFIX}.ldr.cojo"
if [[ "$MA_FILE" != "${OUT_PREFIX}.ma" ]]; then
  cp -f "$MA_FILE" "${OUT_PREFIX}.ma"
fi

# Keep conditioned SNPs in final cma output and mark conditional fields as NA.
FINAL_CMA="${OUT_PREFIX}.cma.cojo"
TMP_CMA="${OUT_PREFIX}.cma.with_cond.tmp"

awk -F'\t' 'BEGIN{OFS="\t"}
FILENAME==ARGV[1] {
  if (FNR==1) { print; next }
  existing[$2]=1
  print
  next
}
FILENAME==ARGV[2] {
  if (FNR==1) next
  ma_refA[$1]=toupper($2)
  ma_freq[$1]=$4
  ma_b[$1]=$5
  ma_se[$1]=$6
  ma_p[$1]=$7
  ma_n[$1]=$8
  next
}
FILENAME==ARGV[3] {
  bim_chr[$2]=$1
  bim_bp[$2]=$4
  next
}
FILENAME==ARGV[4] {
  snp=$1
  if (snp=="" || (snp in existing) || (snp in added)) next
  chr=(snp in bim_chr)?bim_chr[snp]:"NA"
  bp=(snp in bim_bp)?bim_bp[snp]:"NA"
  refA=(snp in ma_refA)?ma_refA[snp]:"NA"
  freq=(snp in ma_freq)?ma_freq[snp]:"NA"
  b=(snp in ma_b)?ma_b[snp]:"NA"
  se=(snp in ma_se)?ma_se[snp]:"NA"
  p=(snp in ma_p)?ma_p[snp]:"NA"
  n=(snp in ma_n)?ma_n[snp]:"NA"
  print chr,snp,bp,refA,freq,b,se,p,n,"NA","NA","NA","NA"
  added[snp]=1
  next
}
' "$FINAL_CMA" "$MA_FILE" "${BFILE_PREFIX}.bim" "$COND_FILE" > "$TMP_CMA"

mv -f "$TMP_CMA" "$FINAL_CMA"

LAST_COND=$(tail -n 1 "$COND_FILE")
printf "ok\tchr=%s\tpos=%s\tlast_cond_snp=%s\tn_cond=%s\n" "$TARGET_CHR" "$TARGET_POS" "$LAST_COND" "$(wc -l < "$COND_FILE")" > "$OUT_DONE"
