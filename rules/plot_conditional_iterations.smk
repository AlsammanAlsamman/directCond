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
for analysis_name, analysis_cfg in TARGET_ANALYSES.items():
    regions_file = analysis_cfg["regions_file"]
    ANALYSIS_REGIONS[analysis_name] = _load_regions(regions_file)

PLOT_INDEX = {}
for analysis_name, regions in ANALYSIS_REGIONS.items():
    for locus_id, locus_data in regions.items():
        PLOT_INDEX[(analysis_name, locus_id)] = {
            "chr": locus_data["chr"],
            "pos": locus_data["pos"],
        }

ALL_PLOT_DONE_TARGETS = [
    os.path.join(RESULTS_BASE, "04_plots", analysis_name, locus_id, "plot.done")
    for analysis_name, regions in ANALYSIS_REGIONS.items()
    for locus_id in sorted(regions.keys())
]

R_MODULE = get_software_module("r")
R_LIBS_USER = str(get_software_value("r", ["params", "r_libs_user"], default=""))


rule plot_conditional_iterations:
    input:
        cojo_done=os.path.join(RESULTS_BASE, "03_cojo", "{analysis}", "{locus_id}", "cojo.done"),
        cojo_final=os.path.join(RESULTS_BASE, "03_cojo", "{analysis}", "{locus_id}", "cojo.cma.cojo"),
        cojo_ma=os.path.join(RESULTS_BASE, "03_cojo", "{analysis}", "{locus_id}", "cojo.ma"),
        ref_bed=os.path.join(RESULTS_BASE, "01_extract_regions", "{analysis}", "{locus_id}", "ref_subset.bed"),
    output:
        done=os.path.join(RESULTS_BASE, "04_plots", "{analysis}", "{locus_id}", "plot.done"),
        pdf=os.path.join(RESULTS_BASE, "04_plots", "{analysis}", "{locus_id}", "conditional_iterations.pdf"),
        png=os.path.join(RESULTS_BASE, "04_plots", "{analysis}", "{locus_id}", "conditional_iterations.png"),
        before_pdf=os.path.join(RESULTS_BASE, "04_plots", "{analysis}", "{locus_id}", "before_conditioning.pdf"),
        after_pdf=os.path.join(RESULTS_BASE, "04_plots", "{analysis}", "{locus_id}", "after_conditioning.pdf"),
        before_png=os.path.join(RESULTS_BASE, "04_plots", "{analysis}", "{locus_id}", "before_conditioning.png"),
        after_png=os.path.join(RESULTS_BASE, "04_plots", "{analysis}", "{locus_id}", "after_conditioning.png"),
    params:
        r_module=R_MODULE,
        r_libs_user=R_LIBS_USER,
        cojo_dir=os.path.join(RESULTS_BASE, "03_cojo", "{analysis}", "{locus_id}"),
        target_chr=lambda wc: PLOT_INDEX[(wc.analysis, wc.locus_id)]["chr"],
        target_pos=lambda wc: PLOT_INDEX[(wc.analysis, wc.locus_id)]["pos"],
    log:
        os.path.join(RESULTS_BASE, "log", "plot_conditioning", "{analysis}", "{locus_id}.log"),
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
        if command -v module >/dev/null 2>&1; then module load {params.r_module}; fi
        if [ -n "{params.r_libs_user}" ]; then export R_LIBS_USER='{params.r_libs_user}'; fi
        Rscript scripts/plot_conditional_iterations_clean.R \
          --cojo-ma {input.cojo_ma} \
          --cojo-final {input.cojo_final} \
          --cojo-dir {params.cojo_dir} \
                    --ref-prefix "$(dirname {input.ref_bed})/ref_subset" \
                    --target-chr {params.target_chr} \
                    --target-pos {params.target_pos} \
          --analysis {wildcards.analysis} \
          --locus-id {wildcards.locus_id} \
          --out-pdf {output.pdf} \
                    --out-png {output.png} \
                    --out-before-pdf {output.before_pdf} \
                    --out-after-pdf {output.after_pdf} \
                                        --out-before-png {output.before_png} \
                                        --out-after-png {output.after_png} \
          --out-done {output.done} \
          > {log} 2>&1
        """
