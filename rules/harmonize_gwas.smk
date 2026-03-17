import os
import re
import sys

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


TARGET_ANALYSES = get_analysis_value(["target_analyses"])
RESULTS_DIR = get_results_dir()
PROJECT_NAME = str(get_analysis_value(["project_name"]))
RESULTS_BASE = os.path.join(RESULTS_DIR, PROJECT_NAME)
RESOURCE_MEM_MB = _resource_value("mem_mb", 32000)
RESOURCE_CORES = _resource_value("cores", 2)
RESOURCE_TIME = _resource_value("time", "00:30:00")

ANALYSIS_REGIONS = {}
for analysis_name, analysis_cfg in TARGET_ANALYSES.items():
    regions_file = analysis_cfg["regions_file"]
    ANALYSIS_REGIONS[analysis_name] = _load_regions(regions_file)

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
