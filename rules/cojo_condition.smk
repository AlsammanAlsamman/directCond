import os
import re
import sys

sys.path.append("utils")
from bioconfigme import get_analysis_value, get_results_dir, get_software_module, get_software_value


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
COJO_INDEX = {}
for analysis_name, analysis_cfg in TARGET_ANALYSES.items():
    regions_file = analysis_cfg["regions_file"]
    regions = _load_regions(regions_file)
    ANALYSIS_REGIONS[analysis_name] = regions
    for locus_id, locus_data in regions.items():
        COJO_INDEX[(analysis_name, locus_id)] = {
            "chr": locus_data["chr"],
            "pos": locus_data["pos"],
            "reference_population": str(analysis_cfg["reference_population"]),
        }

ALL_COJO_DONE_TARGETS = [
    os.path.join(RESULTS_BASE, "03_cojo", analysis_name, locus_id, "cojo.done")
    for analysis_name, regions in ANALYSIS_REGIONS.items()
    for locus_id in sorted(regions.keys())
]

GCTA_MODULE = get_software_module("gcta")
GCTA_BIN = str(get_software_value("gcta", ["binary"], default="gcta"))


rule cojo_condition:
    input:
        harmonize_done=os.path.join(RESULTS_BASE, "02_harmonized", "{analysis}", "{locus_id}", "harmonize.done"),
        gwas_harmonized=os.path.join(RESULTS_BASE, "02_harmonized", "{analysis}", "{locus_id}", "gwas_harmonized.tsv"),
        bed=os.path.join(RESULTS_BASE, "01_extract_regions", "{analysis}", "{locus_id}", "ref_subset.bed"),
        bim=os.path.join(RESULTS_BASE, "01_extract_regions", "{analysis}", "{locus_id}", "ref_subset.bim"),
        fam=os.path.join(RESULTS_BASE, "01_extract_regions", "{analysis}", "{locus_id}", "ref_subset.fam"),
    output:
        done=os.path.join(RESULTS_BASE, "03_cojo", "{analysis}", "{locus_id}", "cojo.done"),
        cojo=os.path.join(RESULTS_BASE, "03_cojo", "{analysis}", "{locus_id}", "cojo.cma.cojo"),
        ma=os.path.join(RESULTS_BASE, "03_cojo", "{analysis}", "{locus_id}", "cojo.ma"),
        cond_snp=os.path.join(RESULTS_BASE, "03_cojo", "{analysis}", "{locus_id}", "cojo.cond.snp"),
    params:
        target_chr=lambda wc: COJO_INDEX[(wc.analysis, wc.locus_id)]["chr"],
        target_pos=lambda wc: COJO_INDEX[(wc.analysis, wc.locus_id)]["pos"],
        gcta_module=GCTA_MODULE,
        gcta_bin=GCTA_BIN,
    log:
        os.path.join(RESULTS_BASE, "log", "cojo_condition", "{analysis}", "{locus_id}.log"),
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
        if command -v module >/dev/null 2>&1; then module load {params.gcta_module}; fi
        bash scripts/run_cojo_condition.sh \
                    --gwas-file {input.gwas_harmonized} \
          --bfile-prefix "$(dirname {input.bed})/ref_subset" \
          --target-chr {params.target_chr} \
          --target-pos {params.target_pos} \
          --out-prefix "$(dirname {output.cojo})/cojo" \
          --out-done {output.done} \
          --gcta-bin {params.gcta_bin} \
          > {log} 2>&1
        """
