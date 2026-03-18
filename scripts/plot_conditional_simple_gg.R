#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ggplot2))

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, required = TRUE, default = NA_character_) {
  idx <- which(args == flag)
  if (length(idx) == 0) {
    if (required) stop(paste0("Missing required argument: ", flag))
    return(default)
  }
  if (idx[1] == length(args)) stop(paste0("Missing value for argument: ", flag))
  args[idx[1] + 1]
}

cojo_ma <- get_arg("--cojo-ma")
cojo_dir <- get_arg("--cojo-dir")
cojo_final <- get_arg("--cojo-final")
analysis_name <- get_arg("--analysis", required = FALSE, default = "analysis")
locus_id <- get_arg("--locus-id", required = FALSE, default = "locus")
out_dir <- get_arg("--out-dir")
out_done <- get_arg("--out-done")

if (!file.exists(cojo_ma)) stop(paste("cojo.ma not found:", cojo_ma))
if (!dir.exists(cojo_dir)) stop(paste("cojo_dir not found:", cojo_dir))
if (!file.exists(cojo_final)) stop(paste("cojo final not found:", cojo_final))

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_done), recursive = TRUE, showWarnings = FALSE)

read_tsv <- function(path) {
  read.table(path, header = TRUE, sep = "\t", quote = "", comment.char = "", stringsAsFactors = FALSE, check.names = FALSE)
}

before <- read_tsv(cojo_ma)
required_before <- c("SNP", "p")
missing_before <- setdiff(required_before, names(before))
if (length(missing_before) > 0) {
  stop("cojo.ma missing required columns: ", paste(missing_before, collapse = ", "))
}

before <- before[, c("SNP", "p")]
before$p <- suppressWarnings(as.numeric(before$p))

iter_files <- list.files(cojo_dir, pattern = "^cojo\\.iter[0-9]+\\.cma\\.cojo$", full.names = TRUE)
iter_ids <- integer(0)
if (length(iter_files) > 0) {
  iter_ids <- as.integer(sub("^cojo\\.iter([0-9]+)\\.cma\\.cojo$", "\\1", basename(iter_files)))
  iter_files <- iter_files[order(iter_ids)]
}

all_after_files <- c(iter_files, cojo_final)
all_labels <- c(
  if (length(iter_ids) > 0) paste0("iter", iter_ids) else character(0),
  "final"
)

n_plots <- 0
for (i in seq_along(all_after_files)) {
  after_file <- all_after_files[i]
  label <- all_labels[i]

  after <- read_tsv(after_file)
  required_after <- c("SNP", "bp", "pC")
  missing_after <- setdiff(required_after, names(after))
  if (length(missing_after) > 0) next

  after <- after[, c("SNP", "bp", "pC")]
  names(after) <- c("SNP", "pos", "p_after")
  after$pos <- suppressWarnings(as.numeric(after$pos))
  after$p_after <- suppressWarnings(as.numeric(after$p_after))

  merged <- merge(after, before, by = "SNP", all.x = TRUE)
  names(merged)[names(merged) == "p"] <- "p_before"

  before_df <- merged[, c("SNP", "pos", "p_before")]
  names(before_df) <- c("SNP", "pos", "p")
  before_df$state <- "before"

  after_df <- merged[, c("SNP", "pos", "p_after")]
  names(after_df) <- c("SNP", "pos", "p")
  after_df$state <- paste0("after_", label)

  plot_df <- rbind(before_df, after_df)
  plot_df <- plot_df[is.finite(plot_df$pos) & is.finite(plot_df$p) & plot_df$p > 0, ]
  if (nrow(plot_df) == 0) next

  p <- ggplot(plot_df, aes(x = pos, y = -log10(p))) +
    geom_point(size = 0.35, alpha = 0.7) +
    geom_hline(yintercept = -log10(5e-8), color = "red", linewidth = 0.5) +
    facet_wrap(~state, nrow = 1, scales = "free_y") +
    labs(
      title = paste0(analysis_name, " | ", locus_id, " | ", label),
      x = "Position (bp)",
      y = "-log10(p)"
    ) +
    theme_bw(base_size = 10)

  out_png <- file.path(out_dir, paste0("before_after_", label, ".png"))
  out_pdf <- file.path(out_dir, paste0("before_after_", label, ".pdf"))

  ggsave(out_png, p, width = 11, height = 4.5, dpi = 180)
  ggsave(out_pdf, p, width = 11, height = 4.5)
  n_plots <- n_plots + 1
}

writeLines(
  paste0("ok\tplots_generated=", n_plots, "\tout_dir=", out_dir),
  con = out_done
)
cat("Generated", n_plots, "simple before/after iteration plots in", out_dir, "\n")
