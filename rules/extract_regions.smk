import os
import re
import sys

sys.path.append("utils")
from bioconfigme import get_analysis_value, get_results_dir, get_software_module


def _load_regions(path):
    rows = {}
    with open(path, "r", encoding="utf-8") as handle:
        lines = [line.strip() for line in handle if line.strip()]
        if not lines:
            raise ValueError(f"Empty regions file: {path}")

        header = re.split(r"\s+", lines[0])
        expected = ["locus_id", "chr", "pos", "window_bp"]
        if header != expected:
            raise ValueError("regions.tsv header must be exactly: locus_id chr pos window_bp")

        for line in lines[1:]:
            parts = re.split(r"\s+", line)
            if len(parts) < 4:
                continue
            locus_id = parts[0].strip()
            rows[locus_id] = {
                "chr": str(parts[1]).strip(),
                "pos": int(parts[2]),
                "window_bp": int(parts[3]),
            }
    return rows


def _resource_value(name, default_value):
    try:
        return int(get_analysis_value(["default_resources", name])) if name in {"mem_mb", "cores"} else str(get_analysis_value(["default_resources", name]))
    except Exception:
        return default_value


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
    regions_file = analysis_cfg["regions_file"]
    regions = _load_regions(regions_file)
    ANALYSIS_REGIONS[analysis_name] = regions
    for locus_id, locus_data in regions.items():
        REGION_INDEX[(analysis_name, locus_id)] = {
            "chr": locus_data["chr"],
            "pos": locus_data["pos"],
            "window_bp": locus_data["window_bp"],
            "gwas_file": str(analysis_cfg["gwas_file"]),
            "reference_population": str(analysis_cfg["reference_population"]),
        }

ALL_DONE_TARGETS = [
    os.path.join(RESULTS_BASE, "01_extract_regions", analysis_name, locus_id, "extract.done")
    for analysis_name, regions in ANALYSIS_REGIONS.items()
    for locus_id in sorted(regions.keys())
]

PLINK_MODULE = get_software_module("plink2")


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
        pos=lambda wc: REGION_INDEX[(wc.analysis, wc.locus_id)]["pos"],
        window_bp=lambda wc: REGION_INDEX[(wc.analysis, wc.locus_id)]["window_bp"],
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
          --pos {params.pos} \
          --window-bp {params.window_bp} \
          --out-gwas-window {output.gwas_window} \
          --out-ref-prefix "$(dirname {output.bed})/ref_subset" \
          --out-done {output.done} \
          > {log} 2>&1
        """
