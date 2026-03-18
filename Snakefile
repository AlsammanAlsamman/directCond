configfile: "configs/analysis.yml"

include: "rules/extract_regions.smk"
include: "rules/harmonize_gwas.smk"
include: "rules/cojo_condition.smk"
include: "rules/plot_conditional_iterations.smk"
include: "rules/indiv_cond_eval.smk"

rule all:
    input:
        ALL_PLOT_DONE_TARGETS


rule all_simple_plots:
    input:
        ALL_PLOT_SIMPLE_DONE_TARGETS
