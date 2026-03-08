#!/usr/bin/env bash
set -euo pipefail

TARGET_CHR=21
TARGET_POS=43855067
REF_PREFIX="ref/g1000_eas_chr21"
SNP_LIST="snp100_21_43855067.list"
OUT_PREFIX="ld_21_43855067_100snps"

awk -v chr="$TARGET_CHR" -v pos="$TARGET_POS" '
	$1==chr {
		d=$4-pos
		if (d<0) d=-d
		print d "\t" $2
	}
' "${REF_PREFIX}.bim" | sort -n | head -10000 | cut -f2 > "$SNP_LIST"

plink --bfile "$REF_PREFIX" \
	--extract "$SNP_LIST" \
	--r square gz \
	--out "$OUT_PREFIX"

echo "Done: ${OUT_PREFIX}.ld.gz"
