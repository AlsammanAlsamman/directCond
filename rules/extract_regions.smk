import os
import re
import sys
import csv
import hashlib

sys.path.append("utils")
from bioconfigme import get_analysis_value, get_results_dir, get_software_module


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
            else:
                pos = int(parts[2])
                window_bp = int(parts[3])
                half_window = window_bp // 2
                start_bp = max(1, pos - half_window)
                end_bp = pos + half_window

            rows[locus_id] = {
                "chr": chrom,
                "start": start_bp,
                "end": end_bp,
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
    snps = []
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            snp = line.strip().replace("\r", "")
            if not snp:
                continue
            snps.append(snp)
    if not snps:
        raise ValueError(f"Empty snpList_file: {path}")
    return snps


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
            raise ValueError(
                "GWAS header must contain SNP/CHR/POS columns for snpList-based extraction"
            )

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

            half_window = int(region_width) // 2
            found[snp] = {
                "chr": chrom,
                "start": max(1, pos - half_window),
                "end": pos + half_window,
            }

            if len(found) == len(target_set):
                break

    missing = [s for s in targets if s not in found]
    if missing:
        preview = ", ".join(missing[:10])
        suffix = " ..." if len(missing) > 10 else ""
        raise ValueError(f"SNPs from snpList_file not found in GWAS: {preview}{suffix}")

    # Keep folder/output ordering consistent with snp list order.
    return {snp: found[snp] for snp in targets}


def _cache_key(*parts):
    text = "|".join(str(p) for p in parts)
    return hashlib.sha1(text.encode("utf-8")).hexdigest()[:12]


def _load_cached_regions(path):
    if not os.path.exists(path):
        return None
    return _load_regions(path)


def _write_cached_regions(path, regions):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as handle:
        handle.write("locus_id\tchr\tstart\tend\n")
        for locus_id, row in regions.items():
            handle.write(f"{locus_id}\t{row['chr']}\t{row['start']}\t{row['end']}\n")


TARGET_ANALYSES = get_analysis_value(["target_analyses"])
RESULTS_DIR = get_results_dir()
PROJECT_NAME = str(get_analysis_value(["project_name"]))
RESULTS_BASE = os.path.join(RESULTS_DIR, PROJECT_NAME)
RESOURCE_MEM_MB = _resource_value("mem_mb", 32000)
RESOURCE_CORES = _resource_value("cores", 2)
RESOURCE_TIME = _resource_value("time", "00:30:00")

ANALYSIS_REGIONS = {}
REGION_INDEX = {}
for analysis_name, analysis_cfg in TARGET_ANALYSES.items():
    if "regions_file" in analysis_cfg and analysis_cfg["regions_file"]:
        regions = _load_regions(str(analysis_cfg["regions_file"]))
    elif "snpList_file" in analysis_cfg and analysis_cfg["snpList_file"]:
        if "gwas_file" not in analysis_cfg or not analysis_cfg["gwas_file"]:
            raise ValueError(f"Analysis {analysis_name}: gwas_file is required when using snpList_file")
        region_width = int(analysis_cfg.get("region_width", 500000))
        gwas_file = str(analysis_cfg["gwas_file"])
        snp_list_file = str(analysis_cfg["snpList_file"])

        cache_id = _cache_key(analysis_name, gwas_file, snp_list_file, region_width)
        cache_file = os.path.join(RESULTS_BASE, "00_cache", "auto_regions", f"{analysis_name}.{cache_id}.tsv")

        use_cache = False
        if os.path.exists(cache_file):
            try:
                cache_mtime = os.path.getmtime(cache_file)
                gwas_mtime = os.path.getmtime(gwas_file)
                snplist_mtime = os.path.getmtime(snp_list_file)
                use_cache = cache_mtime >= max(gwas_mtime, snplist_mtime)
            except Exception:
                use_cache = False

        if use_cache:
            regions = _load_cached_regions(cache_file)
        else:
            regions = _regions_from_snp_list(
                gwas_file=gwas_file,
                snp_list_file=snp_list_file,
                region_width=region_width,
            )
            _write_cached_regions(cache_file, regions)

        if not regions:
            raise ValueError(f"Analysis {analysis_name}: no regions generated from snpList_file")
    else:
        raise ValueError(
            f"Analysis {analysis_name} must define either regions_file or snpList_file (+ gwas_file)"
        )

    ANALYSIS_REGIONS[analysis_name] = regions
    for locus_id, locus_data in regions.items():
        REGION_INDEX[(analysis_name, locus_id)] = {
            "chr": locus_data["chr"],
            "start": locus_data["start"],
            "end": locus_data["end"],
            "gwas_file": str(analysis_cfg["gwas_file"]),
            "reference_population": str(analysis_cfg["reference_population"]),
        }

ALL_DONE_TARGETS = [
    os.path.join(RESULTS_BASE, "01_extract_regions", analysis_name, locus_id, "extract.done")
    for analysis_name, regions in ANALYSIS_REGIONS.items()
    for locus_id in sorted(regions.keys())
]

PLINK_MODULE = get_software_module("plink2")


rule extract_regions_all:
    input:
        ALL_DONE_TARGETS


rule extract_regions:
    output:
        done=os.path.join(RESULTS_BASE, "01_extract_regions", "{analysis}", "{locus_id}", "extract.done"),
        gwas_window=os.path.join(RESULTS_BASE, "01_extract_regions", "{analysis}", "{locus_id}", "gwas_window.tsv"),
        bed=os.path.join(RESULTS_BASE, "01_extract_regions", "{analysis}", "{locus_id}", "ref_subset.bed"),
        bim=os.path.join(RESULTS_BASE, "01_extract_regions", "{analysis}", "{locus_id}", "ref_subset.bim"),
        fam=os.path.join(RESULTS_BASE, "01_extract_regions", "{analysis}", "{locus_id}", "ref_subset.fam"),
    params:
        gwas_file=lambda wc: REGION_INDEX[(wc.analysis, wc.locus_id)]["gwas_file"],
        ref_prefix=lambda wc: get_analysis_value(["ref_panel", REGION_INDEX[(wc.analysis, wc.locus_id)]["reference_population"]]),
        chrom=lambda wc: REGION_INDEX[(wc.analysis, wc.locus_id)]["chr"],
        start=lambda wc: REGION_INDEX[(wc.analysis, wc.locus_id)]["start"],
        end=lambda wc: REGION_INDEX[(wc.analysis, wc.locus_id)]["end"],
        plink_module=PLINK_MODULE,
    log:
        os.path.join(RESULTS_BASE, "log", "extract_regions", "{analysis}", "{locus_id}.log"),
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
        if command -v module >/dev/null 2>&1; then module load {params.plink_module}; fi
        python scripts/extract_region.py \
          --gwas-file {params.gwas_file} \
          --ref-prefix {params.ref_prefix} \
          --chr {params.chrom} \
                    --start {params.start} \
                    --end {params.end} \
          --out-gwas-window {output.gwas_window} \
          --out-ref-prefix "$(dirname {output.bed})/ref_subset" \
          --out-done {output.done} \
          > {log} 2>&1
        """
