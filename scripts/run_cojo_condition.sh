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
    [--cond-snps-file <file_with_one_snp_per_line>] \
    --out-prefix <output_prefix_without_extension> \
    [--out-indiv-summary <summary_tsv>] \
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
  --cond-snps-file File with one SNP ID per line. If provided, run one COJO
                   pass for each SNP individually first, then run one combined
                   conditioning pass on all listed SNPs.
  --out-indiv-summary TSV summary for the per-SNP individual conditioning runs.
                   Default: <out-prefix>.individual.summary.tsv
  --gcta-bin     GCTA executable (default: gcta)
EOF
}

GWAS_FILE=""
BFILE_PREFIX=""
TARGET_CHR=""
TARGET_POS=""
COND_SNPS_FILE=""
OUT_PREFIX=""
OUT_INDIV_SUMMARY=""
OUT_DONE=""
GCTA_BIN="gcta"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --gwas-file) GWAS_FILE="$2"; shift 2 ;;
    --bfile-prefix) BFILE_PREFIX="$2"; shift 2 ;;
    --target-chr) TARGET_CHR="$2"; shift 2 ;;
    --target-pos) TARGET_POS="$2"; shift 2 ;;
    --cond-snps-file) COND_SNPS_FILE="$2"; shift 2 ;;
    --out-prefix) OUT_PREFIX="$2"; shift 2 ;;
    --out-indiv-summary) OUT_INDIV_SUMMARY="$2"; shift 2 ;;
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

if [[ -n "$COND_SNPS_FILE" && ! -f "$COND_SNPS_FILE" ]]; then
  echo "Condition SNP list file not found: $COND_SNPS_FILE" >&2
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
if [[ -z "$OUT_INDIV_SUMMARY" ]]; then
  OUT_INDIV_SUMMARY="${OUT_PREFIX}.individual.summary.tsv"
fi

sanitize_for_path() {
  printf '%s' "$1" | sed 's/[^A-Za-z0-9._-]/_/g'
}

count_significant_snps_in_cma() {
  local cma_file="$1"
  local eval_list_file="$2"

  awk -v thr="5e-8" '
BEGIN { FS="\t"; n=0 }
FNR==NR {
  snp=$1
  gsub(/^[[:space:]]+|[[:space:]]+$/, "", snp)
  if (snp!="") keep[snp]=1
  next
}
FNR==1 {
  for (i=1; i<=NF; i++) {
    h=tolower($i)
    gsub(/^ +| +$/, "", h)
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
' "$eval_list_file" "$cma_file"
}

run_individual_list_conditioning() {
  local eval_list_file="$1"
  local indiv_dir="${OUT_DIR}/individual"
  local bim_chr
  local n_evaluated=0

  mkdir -p "$indiv_dir" "$(dirname "$OUT_INDIV_SUMMARY")"
  printf "snp\tn_significant_after_cond\tstatus\tcma_file\n" > "$OUT_INDIV_SUMMARY"

  bim_chr=$(awk 'NR==1 {print $1; exit}' "${BFILE_PREFIX}.bim")
  if [[ -z "$bim_chr" ]]; then
    echo "Unable to determine chromosome from BIM file: ${BFILE_PREFIX}.bim" >&2
    exit 1
  fi

  while IFS= read -r snp || [[ -n "$snp" ]]; do
    snp=$(printf '%s' "$snp" | tr -d '\r[:space:]')
    [[ -z "$snp" ]] && continue

    local safe_snp
    local snp_dir
    local single_cond_file
    local indiv_prefix
    local gcta_log
    local status="ok"
    local n_sig="NA"

    safe_snp=$(sanitize_for_path "$snp")
    snp_dir="${indiv_dir}/${safe_snp}"
    single_cond_file="${snp_dir}/single.cond.snp"
    indiv_prefix="${snp_dir}/cojo"
    gcta_log="${snp_dir}/gcta.log"

    mkdir -p "$snp_dir"
    printf "%s\n" "$snp" > "$single_cond_file"

    if ! "$GCTA_BIN" \
      --bfile "$BFILE_PREFIX" \
      --chr "$bim_chr" \
      --cojo-file "$MA_FILE" \
      --cojo-cond "$single_cond_file" \
      --out "$indiv_prefix" > "$gcta_log" 2>&1; then
      if grep -qi "collinearity problem" "$gcta_log"; then
        status="collinear"
      else
        status="gcta_failed"
      fi
    elif [[ -f "${indiv_prefix}.cma.cojo" ]]; then
      n_sig=$(count_significant_snps_in_cma "${indiv_prefix}.cma.cojo" "$eval_list_file")
    else
      status="missing_output"
    fi

    printf "%s\t%s\t%s\t%s\n" "$snp" "$n_sig" "$status" "${indiv_prefix}.cma.cojo" >> "$OUT_INDIV_SUMMARY"
    n_evaluated=$((n_evaluated + 1))
  done < "$eval_list_file"

  if [[ $n_evaluated -eq 0 ]]; then
    echo "No SNPs were available for individual conditioning from: $eval_list_file" >&2
    exit 1
  fi
}

append_conditioned_snps_to_final_cma() {
  local final_cma="$1"
  local tmp_cma="${final_cma}.tmp"

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
' "$final_cma" "$MA_FILE" "${BFILE_PREFIX}.bim" "$COND_FILE" > "$tmp_cma"

  mv -f "$tmp_cma" "$final_cma"
}

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

if [[ -n "$COND_SNPS_FILE" ]]; then
  awk '{
    gsub(/\r/, "", $0)
    if (NF > 0) {
      snp=$1
      gsub(/\r/, "", snp)
      print snp
    }
  }' "$COND_SNPS_FILE" > "$COND_FILE"

  if [[ ! -s "$COND_FILE" ]]; then
    echo "Condition SNP list is empty after filtering blank lines: $COND_SNPS_FILE" >&2
    exit 1
  fi

  run_individual_list_conditioning "$COND_FILE"

  if ! "$GCTA_BIN" \
    --bfile "$BFILE_PREFIX" \
    --chr "$TARGET_CHR" \
    --cojo-file "$MA_FILE" \
    --cojo-cond "$COND_FILE" \
    --out "$OUT_PREFIX"; then
    echo "GCTA failed while conditioning on SNP list: $COND_SNPS_FILE" >&2
    exit 1
  fi

  if [[ ! -f "${OUT_PREFIX}.cma.cojo" ]]; then
    echo "Expected COJO output missing: ${OUT_PREFIX}.cma.cojo" >&2
    exit 1
  fi

  append_conditioned_snps_to_final_cma "${OUT_PREFIX}.cma.cojo"

  printf "ok\tmode=list\tchr=%s\tpos=%s\tn_cond=%s\tcond_source=%s\tindiv_summary=%s\n" \
    "$TARGET_CHR" "$TARGET_POS" "$(wc -l < "$COND_FILE")" "$COND_SNPS_FILE" "$OUT_INDIV_SUMMARY" > "$OUT_DONE"
  exit 0
fi

printf "snp\tn_significant_after_cond\tstatus\tcma_file\n" > "$OUT_INDIV_SUMMARY"

P_THRESHOLD="5e-8"
MAX_ITERS=50
ITER=1
FINAL_PREFIX=""
TOTAL_SKIPPED_FOR_COLLINEARITY=0
INITIAL_CANDIDATES_FILE="${OUT_PREFIX}.initial_candidates.tsv"
REMAINING_CANDIDATES_FILE="${OUT_PREFIX}.remaining_candidates.tsv"

build_initial_candidate_list() {
  local ma_file="$1"
  local out_file="$2"

  awk -v thr="$P_THRESHOLD" '
BEGIN { FS="\t"; OFS="\t" }
NR==1 { next }
{
  snp=$1
  raw_p=$7
  if (snp=="" || raw_p=="" || toupper(raw_p)=="NA") next
  p=raw_p + 0
  if (p < thr) print snp, p
}
' "$ma_file" | sort -k2,2g > "$out_file"
}

prune_remaining_candidates() {
  local remaining_file="$1"
  local cma_file="$2"
  local tmp_file="${remaining_file}.tmp"

  awk -v thr="$P_THRESHOLD" '
BEGIN { FS="\t"; OFS="\t" }
FNR==NR {
  keep[$1]=1
  next
}
FNR==1 {
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
  raw_pc=$c_pc
  if (!(snp in keep) || raw_pc=="" || toupper(raw_pc)=="NA") next
  pc=raw_pc + 0
  if (pc < thr) print snp, pc
}
' "$remaining_file" "$cma_file" | sort -k2,2g > "$tmp_file"

  mv -f "$tmp_file" "$remaining_file"
}

drop_snp_from_file() {
  local file_path="$1"
  local snp_to_drop="$2"
  local tmp_file="${file_path}.tmp"

  awk -v snp="$snp_to_drop" '$1 != snp { print }' "$file_path" > "$tmp_file"
  mv -f "$tmp_file" "$file_path"
}

peek_next_candidate() {
  local file_path="$1"
  awk 'NF > 0 { print $1; exit }' "$file_path"
}

run_gcta_attempt() {
  local cond_file="$1"
  local iter_prefix="$2"
  local attempt_log="${iter_prefix}.gcta.log"

  rm -f "${iter_prefix}.cma.cojo" "${iter_prefix}.jma.cojo" "${iter_prefix}.ldr.cojo" "${iter_prefix}.given.cojo"

  if "$GCTA_BIN" \
    --bfile "$BFILE_PREFIX" \
    --chr "$TARGET_CHR" \
    --cojo-file "$MA_FILE" \
    --cojo-cond "$cond_file" \
    --out "$iter_prefix" > "$attempt_log" 2>&1; then
    cat "$attempt_log"
    return 0
  fi

  cat "$attempt_log"
  return 1
}

build_initial_candidate_list "$MA_FILE" "$INITIAL_CANDIDATES_FILE"

if [[ ! -s "$INITIAL_CANDIDATES_FILE" ]]; then
  echo "No SNPs with p < ${P_THRESHOLD} found in $MA_FILE" >&2
  exit 1
fi

cp -f "$INITIAL_CANDIDATES_FILE" "$REMAINING_CANDIDATES_FILE"

FIRST_SNP=$(peek_next_candidate "$REMAINING_CANDIDATES_FILE")
if [[ -z "$FIRST_SNP" ]]; then
  echo "Could not determine initial conditioning SNP from $INITIAL_CANDIDATES_FILE" >&2
  exit 1
fi

printf "%s\n" "$FIRST_SNP" > "$COND_FILE"
drop_snp_from_file "$REMAINING_CANDIDATES_FILE" "$FIRST_SNP"
LAST_ADDED_SNP=""
echo "Initial conditioning SNP: $FIRST_SNP" >&2

while [[ $ITER -le $MAX_ITERS ]]; do
  ITER_PREFIX="${OUT_PREFIX}.iter${ITER}"

  ATTEMPT_OK=0
  if run_gcta_attempt "$COND_FILE" "$ITER_PREFIX"; then
    ATTEMPT_OK=1
  fi

  if [[ ! -f "${ITER_PREFIX}.cma.cojo" ]]; then
    if grep -qi "collinearity problem" "${ITER_PREFIX}.gcta.log"; then
      if [[ -n "$LAST_ADDED_SNP" ]]; then
        echo "Iteration ${ITER}: skipping SNP ${LAST_ADDED_SNP} due to collinearity." >&2
        TOTAL_SKIPPED_FOR_COLLINEARITY=$((TOTAL_SKIPPED_FOR_COLLINEARITY + 1))
        grep -Fxv "$LAST_ADDED_SNP" "$COND_FILE" > "${COND_FILE}.tmp"
        mv -f "${COND_FILE}.tmp" "$COND_FILE"
        LAST_ADDED_SNP=""

        NEXT_SNP=$(peek_next_candidate "$REMAINING_CANDIDATES_FILE")
        if [[ -z "$NEXT_SNP" ]]; then
          echo "Stopping iterative conditioning at iteration ${ITER}: no remaining non-conditioned SNPs from the initial significant list." >&2
          break
        fi

        drop_snp_from_file "$REMAINING_CANDIDATES_FILE" "$NEXT_SNP"
        printf "%s\n" "$NEXT_SNP" >> "$COND_FILE"
        LAST_ADDED_SNP="$NEXT_SNP"
        echo "Iteration ${ITER}: trying next SNP ${NEXT_SNP} from the fixed candidate list." >&2
        continue
      fi

      echo "Stopping iterative conditioning at iteration ${ITER}: current conditioned SNP set is collinear." >&2
      break
    fi

    if [[ $ATTEMPT_OK -eq 0 ]]; then
      echo "GCTA failed for iteration ${ITER}." >&2
    fi
    echo "Expected COJO output missing: ${ITER_PREFIX}.cma.cojo" >&2
    exit 1
  fi

  FINAL_PREFIX="$ITER_PREFIX"
  LAST_ADDED_SNP=""

  if [[ -s "$REMAINING_CANDIDATES_FILE" ]]; then
    prune_remaining_candidates "$REMAINING_CANDIDATES_FILE" "${ITER_PREFIX}.cma.cojo"
  fi

  NEXT_SNP=$(peek_next_candidate "$REMAINING_CANDIDATES_FILE")
  if [[ -z "$NEXT_SNP" ]]; then
    echo "Stopping iterative conditioning at iteration ${ITER}: no SNP from the initial significant list remains below ${P_THRESHOLD}." >&2
    break
  fi

  drop_snp_from_file "$REMAINING_CANDIDATES_FILE" "$NEXT_SNP"
  printf "%s\n" "$NEXT_SNP" >> "$COND_FILE"
  LAST_ADDED_SNP="$NEXT_SNP"
  echo "Iteration ${ITER}: adding next SNP ${NEXT_SNP} from the fixed candidate list." >&2
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

append_conditioned_snps_to_final_cma "${OUT_PREFIX}.cma.cojo"

LAST_COND=$(tail -n 1 "$COND_FILE")
printf "ok\tmode=iterative\tchr=%s\tpos=%s\tlast_cond_snp=%s\tn_cond=%s\tskipped_collinear=%s\tinitial_candidates=%s\tremaining_candidates=%s\n" "$TARGET_CHR" "$TARGET_POS" "$LAST_COND" "$(wc -l < "$COND_FILE")" "$TOTAL_SKIPPED_FOR_COLLINEARITY" "$(wc -l < "$INITIAL_CANDIDATES_FILE")" "$(wc -l < "$REMAINING_CANDIDATES_FILE")" > "$OUT_DONE"
