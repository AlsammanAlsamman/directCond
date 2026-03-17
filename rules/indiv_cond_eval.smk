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

INDIV_ANALYSIS_REGIONS = {}
for analysis_name, analysis_cfg in TARGET_ANALYSES.items():
    regions_file = analysis_cfg["regions_file"]
    INDIV_ANALYSIS_REGIONS[analysis_name] = _load_regions(regions_file)

ALL_INDIV_COND_DONE_TARGETS = [
    os.path.join(RESULTS_BASE, "05_indiv_cond", analysis_name, locus_id, "indiv_cond.done")
    for analysis_name, regions in INDIV_ANALYSIS_REGIONS.items()
    for locus_id in sorted(regions.keys())
]

GCTA_MODULE_INDIV = get_software_module("gcta")
GCTA_BIN_INDIV = str(get_software_value("gcta", ["binary"], default="gcta"))


rule indiv_cond_eval:
    input:
        cojo_done=os.path.join(RESULTS_BASE, "03_cojo", "{analysis}", "{locus_id}", "cojo.done"),
        cond_snp=os.path.join(RESULTS_BASE, "03_cojo", "{analysis}", "{locus_id}", "cojo.cond.snp"),
        eval_snp_list=lambda wc: str(
            TARGET_ANALYSES[wc.analysis].get(
                "snpList_file",
                os.path.join(RESULTS_BASE, "03_cojo", wc.analysis, wc.locus_id, "cojo.cond.snp"),
            )
        ),
        cojo_ma=os.path.join(RESULTS_BASE, "03_cojo", "{analysis}", "{locus_id}", "cojo.ma"),
        cojo_final=os.path.join(RESULTS_BASE, "03_cojo", "{analysis}", "{locus_id}", "cojo.cma.cojo"),
        bed=os.path.join(RESULTS_BASE, "01_extract_regions", "{analysis}", "{locus_id}", "ref_subset.bed"),
        bim=os.path.join(RESULTS_BASE, "01_extract_regions", "{analysis}", "{locus_id}", "ref_subset.bim"),
        fam=os.path.join(RESULTS_BASE, "01_extract_regions", "{analysis}", "{locus_id}", "ref_subset.fam"),
    output:
        done=os.path.join(RESULTS_BASE, "05_indiv_cond", "{analysis}", "{locus_id}", "indiv_cond.done"),
        summary=os.path.join(RESULTS_BASE, "05_indiv_cond", "{analysis}", "{locus_id}", "indiv_cond.summary.tsv"),
    params:
        gcta_module=GCTA_MODULE_INDIV,
        gcta_bin=GCTA_BIN_INDIV,
        out_dir=os.path.join(RESULTS_BASE, "05_indiv_cond", "{analysis}", "{locus_id}"),
    log:
        os.path.join(RESULTS_BASE, "log", "indiv_cond_eval", "{analysis}", "{locus_id}.log"),
    resources:
        mem_mb=131072,
        cores=RESOURCE_CORES,
        time=RESOURCE_TIME,
    wildcard_constraints:
        analysis=r"[A-Za-z0-9_\-\.]+",
        locus_id=r"[A-Za-z0-9_\-\.]+",
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{params.out_dir}" "$(dirname {log})"
        if command -v module >/dev/null 2>&1; then module load {params.gcta_module}; fi
        bash scripts/run_indiv_cond_eval.sh \
          --cond-snp   {input.cond_snp} \
                    --eval-snp-list {input.eval_snp_list} \
          --cojo-ma    {input.cojo_ma} \
          --cojo-final {input.cojo_final} \
          --bfile-prefix "$(dirname {input.bed})/ref_subset" \
          --out-dir    {params.out_dir} \
          --out-summary {output.summary} \
          --out-done   {output.done} \
          --gcta-bin   {params.gcta_bin} \
          > {log} 2>&1
        """
