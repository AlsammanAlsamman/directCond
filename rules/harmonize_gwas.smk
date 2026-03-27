import os
import re
import sys
import csv
import hashlib

sys.path.append("utils")
from bioconfigme import get_analysis_value, get_results_dir


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
            rows[locus_id] = {
                "chr": str(parts[1]).strip(),
                "pos": int(parts[2]),
            }
    return rows


def _resource_value(name, default_value):
    try:
        return int(get_analysis_value(["default_resources", name])) if name in {"mem_mb", "cores"} else str(get_analysis_value(["default_resources", name]))
    except Exception:
        return default_value


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


def _regions_from_snp_list(gwas_file, snp_list_file, region_width):
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
            if not chrom:
                continue
            half = int(region_width) // 2
            found[snp] = {"chr": chrom, "pos": pos}
            if len(found) == len(target_set):
                break

    missing = [s for s in targets if s not in found]
    if missing:
        preview = ", ".join(missing[:10])
        suffix = " ..." if len(missing) > 10 else ""
        raise ValueError(f"SNPs from snpList_file not found in GWAS: {preview}{suffix}")

    return {snp: found[snp] for snp in targets}


def _cache_key(*parts):
    text = "|".join(str(p) for p in parts)
    return hashlib.sha1(text.encode("utf-8")).hexdigest()[:12]


TARGET_ANALYSES = get_analysis_value(["target_analyses"])
RESULTS_DIR = get_results_dir()
PROJECT_NAME = str(get_analysis_value(["project_name"]))
RESULTS_BASE = os.path.join(RESULTS_DIR, PROJECT_NAME)
RESOURCE_MEM_MB = _resource_value("mem_mb", 32000)
RESOURCE_CORES = _resource_value("cores", 2)
RESOURCE_TIME = _resource_value("time", "00:30:00")

ANALYSIS_REGIONS = {}
for analysis_name, analysis_cfg in TARGET_ANALYSES.items():
    if analysis_cfg.get("regions_file"):
        ANALYSIS_REGIONS[analysis_name] = _load_regions(str(analysis_cfg["regions_file"]))
    elif analysis_cfg.get("snpList_file"):
        gwas_file = str(analysis_cfg["gwas_file"])
        snp_list_file = str(analysis_cfg["snpList_file"])
        region_width = int(analysis_cfg.get("region_width", 500000))

        cache_id = _cache_key(analysis_name, gwas_file, snp_list_file, region_width)
        cache_file = os.path.join(RESULTS_BASE, "00_cache", "auto_regions", f"{analysis_name}.{cache_id}.tsv")

        if os.path.exists(cache_file):
            ANALYSIS_REGIONS[analysis_name] = _load_regions(cache_file)
        else:
            ANALYSIS_REGIONS[analysis_name] = _regions_from_snp_list(gwas_file, snp_list_file, region_width)
    else:
        raise ValueError(f"Analysis {analysis_name} must define regions_file or snpList_file")

ALL_HARMONIZE_DONE_TARGETS = [
    os.path.join(RESULTS_BASE, "02_harmonized", analysis_name, locus_id, "harmonize.done")
    for analysis_name, regions in ANALYSIS_REGIONS.items()
    for locus_id in sorted(regions.keys())
]


rule harmonize_gwas_with_ref:
    input:
        extract_done=os.path.join(RESULTS_BASE, "01_extract_regions", "{analysis}", "{locus_id}", "extract.done"),
        gwas_window=os.path.join(RESULTS_BASE, "01_extract_regions", "{analysis}", "{locus_id}", "gwas_window.tsv"),
        ref_bim=os.path.join(RESULTS_BASE, "01_extract_regions", "{analysis}", "{locus_id}", "ref_subset.bim"),
    output:
        done=os.path.join(RESULTS_BASE, "02_harmonized", "{analysis}", "{locus_id}", "harmonize.done"),
        gwas_harmonized=os.path.join(RESULTS_BASE, "02_harmonized", "{analysis}", "{locus_id}", "gwas_harmonized.tsv"),
    log:
        os.path.join(RESULTS_BASE, "log", "harmonize_gwas", "{analysis}", "{locus_id}.log"),
    resources:
        mem_mb=RESOURCE_MEM_MB,
        cores=RESOURCE_CORES,
        time=RESOURCE_TIME,
    wildcard_constraints:
        analysis=r"[A-Za-z0-9_\-\.]+",
        locus_id=r"[A-Za-z0-9_\-\.]+",
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.done})" "$(dirname {log})"
        python scripts/harmonize_gwas_with_ref.py \
          --gwas-window {input.gwas_window} \
          --ref-bim {input.ref_bim} \
          --out-gwas {output.gwas_harmonized} \
          --out-done {output.done} \
          > {log} 2>&1
        """
