#!/usr/bin/env Rscript

# Extract variants within a window around a target SNP position.
# Defaults:
#   input file: ASIANq2.tsv
#   target SNP: 21:43855067
#   window: 400000 bp (±400 kb)

args <- commandArgs(trailingOnly = TRUE)

input_file <- if (length(args) >= 1) args[[1]] else "ASIANq2.tsv"
target_chr <- if (length(args) >= 2) args[[2]] else "21"
target_pos <- if (length(args) >= 3) as.numeric(args[[3]]) else 43855067
window_bp  <- if (length(args) >= 4) as.numeric(args[[4]]) else 400000

start_pos <- target_pos - window_bp
end_pos   <- target_pos + window_bp

message("Input file: ", input_file)
message("Target: ", target_chr, ":", target_pos)
message("Window: ±", window_bp, " bp (", start_pos, " to ", end_pos, ")")

if (requireNamespace("data.table", quietly = TRUE)) {
  dt <- data.table::fread(input_file, sep = "\t", header = TRUE, showProgress = TRUE)
  cn <- names(dt)

  if (!("chrom" %in% cn && "pos" %in% cn)) {
    stop("Expected columns 'chrom' and 'pos' were not found.")
  }

  dt[, chrom_clean := gsub("^chr", "", as.character(chrom), ignore.case = TRUE)]
  dt[, pos_num := as.numeric(pos)]

  result <- dt[chrom_clean == as.character(target_chr) & pos_num >= start_pos & pos_num <= end_pos]
  before <- dt[chrom_clean == as.character(target_chr) & pos_num >= start_pos & pos_num < target_pos]
  after  <- dt[chrom_clean == as.character(target_chr) & pos_num > target_pos & pos_num <= end_pos]

  out_prefix <- paste0("chr", target_chr, "_", target_pos, "_pm", window_bp, "bp")
  data.table::fwrite(result[, !c("chrom_clean", "pos_num")], paste0(out_prefix, "_window.tsv"), sep = "\t")
  data.table::fwrite(before[, !c("chrom_clean", "pos_num")], paste0(out_prefix, "_before.tsv"), sep = "\t")
  data.table::fwrite(after[, !c("chrom_clean", "pos_num")], paste0(out_prefix, "_after.tsv"), sep = "\t")

  message("Rows in full window: ", nrow(result))
  message("Rows before target: ", nrow(before))
  message("Rows after target: ", nrow(after))
  message("Output files:")
  message("  ", paste0(out_prefix, "_window.tsv"))
  message("  ", paste0(out_prefix, "_before.tsv"))
  message("  ", paste0(out_prefix, "_after.tsv"))
} else {
  stop("Package 'data.table' is required for this large file workflow. Install with: install.packages('data.table')")
}

