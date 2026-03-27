import os
import re
import sys
import csv
import hashlib

sys.path.append("utils")
from bioconfigme import get_analysis_value, get_results_dir, get_software_module, get_software_value


def _load_regions(path):
    rows = {}
    with open(path, "r", encoding="utf-8") as handle:
        lines = [line.strip() for line in handle if line.strip()]
        if not lines:
            raise ValueError(f"Empty regions file: {path}")

        header = re.split(r"\s+", lines[0])
        expected_old = ["locus_id", "chr", "pos", "window_bp"]
        expected_new = ["locus_id", "chr", "start", "end"]
        if header not in (expected_old, expected_new):
            raise ValueError("regions.tsv header must be exactly: locus_id chr start end (or legacy: locus_id chr pos window_bp)")

        for line in lines[1:]:
            parts = re.split(r"\s+", line)
            if len(parts) < 4:
                continue
            locus_id = parts[0].strip()
            chrom = str(parts[1]).strip()
            if header == expected_new:
                start_bp = int(parts[2])
                end_bp = int(parts[3])
                pos = (start_bp + end_bp) // 2
            else:
                pos = int(parts[2])

            rows[locus_id] = {
                "chr": chrom,
                "pos": pos,
            }
    return rows


def _find_col(fieldnames, candidates):
    lookup = {str(name).strip().lower(): name for name in fieldnames}
    for cand in candidates:
        if cand in lookup:
            return lookup[cand]
    return None


def _load_snp_list(path):
    rows = []
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            snp = line.strip().replace("\r", "")
            if snp:
                rows.append(snp)
    if not rows:
        raise ValueError(f"Empty snpList_file: {path}")
    return rows


def _regions_from_snp_list(gwas_file, snp_list_file):
    targets = _load_snp_list(snp_list_file)
    target_set = set(targets)
    found = {}

    with open(gwas_file, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"GWAS file has no header: {gwas_file}")

        snp_col = _find_col(reader.fieldnames, ["snp", "rsid", "varid", "marker", "id"])
        chr_col = _find_col(reader.fieldnames, ["chr", "chrom", "chromosome"])
        pos_col = _find_col(reader.fieldnames, ["pos", "position", "bp"])
        if snp_col is None or chr_col is None or pos_col is None:
            raise ValueError("GWAS header must contain SNP/CHR/POS columns for snpList-based extraction")

        for row in reader:
            snp = str(row.get(snp_col, "")).strip()
            if snp not in target_set or snp in found:
                continue
            chrom = str(row.get(chr_col, "")).strip()
            try:
                pos = int(float(row.get(pos_col, "")))
            except Exception:
                continue
            if chrom:
                found[snp] = {"chr": chrom, "pos": pos}
            if len(found) == len(target_set):
                break

    missing = [s for s in targets if s not in found]
    if missing:
        preview = ", ".join(missing[:10])
        suffix = " ..." if len(missing) > 10 else ""
        raise ValueError(f"SNPs from snpList_file not found in GWAS: {preview}{suffix}")

    return {snp: found[snp] for snp in targets}


def _regions_from_locus_chr_startend(path):
    regions = {}
    with open(path, "r", encoding="utf-8") as handle:
        lines = [line.strip() for line in handle if line.strip()]
    if not lines:
        raise ValueError(f"Empty region file: {path}")

    for idx, line in enumerate(lines):
        parts = re.split(r"\s+", line)
        if len(parts) < 4:
            continue

        if idx == 0:
            h0 = parts[0].strip().lower()
            h1 = parts[1].strip().lower()
            h2 = parts[2].strip().lower()
            h3 = parts[3].strip().lower()
            if h0 in {"locus", "locus_id"} and h1 == "chr" and h2 == "start" and h3 == "end":
                continue

        locus_id = parts[0].strip()
        chrom = str(parts[1]).strip()
        start = int(float(parts[2]))
        end = int(float(parts[3]))
        if start > end:
            start, end = end, start
        pos = (start + end) // 2

        key = locus_id
        if key in regions:
            key = f"{locus_id}_{chrom}_{start}_{end}"

        regions[key] = {"chr": chrom, "pos": pos}

    if not regions:
        raise ValueError(f"No valid regions parsed from: {path}")
    return regions


def _cache_key(*parts):
    text = "|".join(str(p) for p in parts)
    return hashlib.sha1(text.encode("utf-8")).hexdigest()[:12]


TARGET_ANALYSES = get_analysis_value(["target_analyses"])
RESULTS_DIR = get_results_dir()
PROJECT_NAME = str(get_analysis_value(["project_name"]))
RESULTS_BASE = os.path.join(RESULTS_DIR, PROJECT_NAME)

ANALYSIS_REGIONS = {}
for analysis_name, analysis_cfg in TARGET_ANALYSES.items():
    if analysis_cfg.get("regions_file"):
        regions = _load_regions(str(analysis_cfg["regions_file"]))
    elif analysis_cfg.get("snpList_file"):
        snp_list_file = str(analysis_cfg["snpList_file"])
        region_file_format = str(analysis_cfg.get("region_file_format", "")).strip().lower()
        if region_file_format == "locus_chr_startend":
            regions = _regions_from_locus_chr_startend(snp_list_file)
        else:
            gwas_file = str(analysis_cfg["gwas_file"])
            region_width = int(analysis_cfg.get("region_width", 500000))
            cache_id = _cache_key(analysis_name, gwas_file, snp_list_file, region_width)
            cache_file = os.path.join(RESULTS_BASE, "00_cache", "auto_regions", f"{analysis_name}.{cache_id}.tsv")
            if os.path.exists(cache_file):
                regions = _load_regions(cache_file)
            else:
                regions = _regions_from_snp_list(gwas_file, snp_list_file)
    else:
        raise ValueError(f"Analysis {analysis_name} must define regions_file or snpList_file")

    ANALYSIS_REGIONS[analysis_name] = regions

R_MODULE = get_software_module("r")
R_LIBS_USER = str(get_software_value("r", ["params", "r_libs_user"], default=""))

ALL_CONDITIONED_SNP_REPORT_XLSX = [
    os.path.join(RESULTS_BASE, "05_reports", analysis_name, "conditioned_snps.xlsx")
    for analysis_name in sorted(ANALYSIS_REGIONS.keys())
]


rule conditioned_snp_excel_report:
    input:
        ALL_CONDITIONED_SNP_REPORT_XLSX


rule conditioned_snp_excel_report_per_analysis:
    input:
        gwas_file=lambda wc: str(TARGET_ANALYSES[wc.analysis]["gwas_file"]),
        cojo_done=lambda wc: [
            os.path.join(RESULTS_BASE, "03_cojo", wc.analysis, locus_id, "cojo.done")
            for locus_id in sorted(ANALYSIS_REGIONS[wc.analysis].keys())
        ],
        cojo_final=lambda wc: [
            os.path.join(RESULTS_BASE, "03_cojo", wc.analysis, locus_id, "cojo.cma.cojo")
            for locus_id in sorted(ANALYSIS_REGIONS[wc.analysis].keys())
        ],
        cond_snp=lambda wc: [
            os.path.join(RESULTS_BASE, "03_cojo", wc.analysis, locus_id, "cojo.cond.snp")
            for locus_id in sorted(ANALYSIS_REGIONS[wc.analysis].keys())
        ],
    output:
        xlsx=os.path.join(RESULTS_BASE, "05_reports", "{analysis}", "conditioned_snps.xlsx"),
    params:
        analysis_name=lambda wc: wc.analysis,
        cojo_dir=os.path.join(RESULTS_BASE, "03_cojo", "{analysis}"),
        r_module=R_MODULE,
        r_libs_user=R_LIBS_USER,
    log:
        os.path.join(RESULTS_BASE, "log", "conditioned_snp_report", "{analysis}.log"),
    resources:
        mem_mb=64000,
        cores=1,
        time="01:00:00",
    wildcard_constraints:
        analysis=r"[A-Za-z0-9_\-\.]+",
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.xlsx})" "$(dirname {log})"
        if command -v module >/dev/null 2>&1; then module load {params.r_module}; fi
        if [ -n "{params.r_libs_user}" ]; then export R_LIBS_USER='{params.r_libs_user}'; fi
        Rscript scripts/export_conditioned_snps_excel.R \
          --analysis {params.analysis_name} \
          --gwas-file {input.gwas_file} \
          --cojo-dir {params.cojo_dir} \
          --out-xlsx {output.xlsx} \
          > {log} 2>&1
        """