#!/usr/bin/env Rscript

apply_env_libpaths <- function() {
  libs <- Sys.getenv("R_LIBS_USER", unset = "")
  if (!nzchar(libs)) return(invisible(NULL))
  parts <- strsplit(libs, .Platform$path.sep, fixed = TRUE)[[1]]
  parts <- parts[nzchar(parts)]
  parts <- path.expand(parts)
  .libPaths(unique(c(parts, .libPaths())))
  invisible(NULL)
}

apply_env_libpaths()

HAVE_ENSDB <- requireNamespace("EnsDb.Hsapiens.v75", quietly = TRUE) &&
  requireNamespace("AnnotationFilter", quietly = TRUE) &&
  requireNamespace("ensembldb", quietly = TRUE)

suppressPackageStartupMessages(library(ggplot2))
if (HAVE_ENSDB) {
  suppressPackageStartupMessages({
    library(AnnotationFilter)
    library(EnsDb.Hsapiens.v75)
  })
}

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
cond_snp <- get_arg("--cond-snp")
indiv_summary <- get_arg("--indiv-summary", required = FALSE, default = NA_character_)
analysis_name <- get_arg("--analysis", required = FALSE, default = "analysis")
locus_id <- get_arg("--locus-id", required = FALSE, default = "locus")
out_dir <- get_arg("--out-dir")
out_done <- get_arg("--out-done")

if (!file.exists(cojo_ma)) stop(paste("cojo.ma not found:", cojo_ma))
if (!dir.exists(cojo_dir)) stop(paste("cojo_dir not found:", cojo_dir))
if (!file.exists(cojo_final)) stop(paste("cojo final not found:", cojo_final))
if (!file.exists(cond_snp)) stop(paste("cond.snp not found:", cond_snp))
if (!is.na(indiv_summary) && nzchar(indiv_summary) && !file.exists(indiv_summary)) {
  stop(paste("individual summary not found:", indiv_summary))
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_done), recursive = TRUE, showWarnings = FALSE)

read_tsv <- function(path) {
  read.table(
    path,
    header = TRUE,
    sep = "\t",
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE,
    check.names = FALSE,
    fill = TRUE,
    blank.lines.skip = TRUE
  )
}

pick_first_col <- function(nms, candidates) {
  hit <- intersect(candidates, nms)
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

assign_gene_tracks <- function(gene_df) {
  if (nrow(gene_df) == 0) return(gene_df)

  gene_df <- gene_df[order(gene_df$start, gene_df$end), , drop = FALSE]
  track_ends <- numeric(0)
  gene_df$track <- integer(nrow(gene_df))

  for (i in seq_len(nrow(gene_df))) {
    placed <- FALSE
    if (length(track_ends) > 0) {
      for (track_idx in seq_along(track_ends)) {
        if (gene_df$start[i] > track_ends[track_idx]) {
          gene_df$track[i] <- track_idx
          track_ends[track_idx] <- gene_df$end[i]
          placed <- TRUE
          break
        }
      }
    }
    if (!placed) {
      track_ends <- c(track_ends, gene_df$end[i])
      gene_df$track[i] <- length(track_ends)
    }
  }

  gene_df
}

fetch_gene_track <- function(chrom, xmin, xmax) {
  if (!HAVE_ENSDB || !is.finite(chrom) || !is.finite(xmin) || !is.finite(xmax)) {
    return(data.frame())
  }

  gene_df <- tryCatch({
    raw_genes <- ensembldb::genes(
      EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75,
      filter = AnnotationFilter::SeqNameFilter(as.character(chrom)),
      return.type = "DataFrame"
    )
    as.data.frame(raw_genes)
  }, error = function(e) {
    data.frame()
  })

  if (nrow(gene_df) == 0) return(gene_df)

  start_col <- pick_first_col(names(gene_df), c("start", "gene_seq_start", "tx_start"))
  end_col <- pick_first_col(names(gene_df), c("end", "gene_seq_end", "tx_end"))
  name_col <- pick_first_col(names(gene_df), c("gene_name", "symbol", "gene", "gene_id"))
  strand_col <- pick_first_col(names(gene_df), c("strand", "gene_strand"))

  if (is.na(start_col) || is.na(end_col) || is.na(name_col)) return(data.frame())

  out <- data.frame(
    gene = as.character(gene_df[[name_col]]),
    start = suppressWarnings(as.numeric(gene_df[[start_col]])),
    end = suppressWarnings(as.numeric(gene_df[[end_col]])),
    strand = if (!is.na(strand_col)) as.character(gene_df[[strand_col]]) else "+",
    stringsAsFactors = FALSE
  )
  out <- out[is.finite(out$start) & is.finite(out$end) & nzchar(out$gene), , drop = FALSE]
  out <- out[out$end >= xmin & out$start <= xmax, , drop = FALSE]
  if (nrow(out) == 0) return(out)

  out$start <- pmax(out$start, xmin)
  out$end <- pmin(out$end, xmax)
  out$mid <- (out$start + out$end) / 2
  assign_gene_tracks(out)
}

before <- read_tsv(cojo_ma)
required_before <- c("SNP", "p")
missing_before <- setdiff(required_before, names(before))
if (length(missing_before) > 0) {
  stop("cojo.ma missing required columns: ", paste(missing_before, collapse = ", "))
}

position_map <- read_tsv(cojo_final)
required_position <- c("SNP", "Chr", "bp")
missing_position <- setdiff(required_position, names(position_map))
if (length(missing_position) > 0) {
  stop("cojo final missing required columns for positions: ", paste(missing_position, collapse = ", "))
}

position_map <- unique(position_map[, c("SNP", "Chr", "bp")])
position_map$Chr <- suppressWarnings(as.numeric(position_map$Chr))
position_map$bp <- suppressWarnings(as.numeric(position_map$bp))
locus_chrom <- unique(position_map$Chr[is.finite(position_map$Chr)])
if (length(locus_chrom) == 0) {
  locus_chrom <- NA_real_
} else {
  locus_chrom <- locus_chrom[1]
}

before <- merge(before[, c("SNP", "p")], position_map[, c("SNP", "bp")], by = "SNP", all.x = TRUE)
before$p <- suppressWarnings(as.numeric(before$p))
before$pos <- suppressWarnings(as.numeric(before$bp))
if (any(!is.finite(before$pos))) {
  fallback_idx <- which(!is.finite(before$pos))
  before$pos[fallback_idx] <- seq_len(length(fallback_idx))
}
before$bp <- NULL

cond_snps <- readLines(cond_snp, warn = FALSE)
cond_snps <- trimws(cond_snps)
cond_snps <- cond_snps[nzchar(cond_snps)]

build_plot_df <- function(data, p_col, state_label) {
  plot_df <- data[, c("SNP", "pos", p_col)]
  names(plot_df) <- c("SNP", "pos", "p")
  plot_df$state <- state_label
  plot_df[is.finite(plot_df$pos) & is.finite(plot_df$p) & plot_df$p > 0, ]
}

load_after_file <- function(path) {
  after <- read_tsv(path)
  required_after <- c("SNP", "bp", "pC")
  missing_after <- setdiff(required_after, names(after))
  if (length(missing_after) > 0) {
    stop(paste("COJO result missing required columns:", paste(missing_after, collapse = ", "), "in", path))
  }
  after <- after[, c("SNP", "bp", "pC")]
  names(after) <- c("SNP", "pos", "p_after")
  after$pos <- suppressWarnings(as.numeric(after$pos))
  after$p_after <- suppressWarnings(as.numeric(after$p_after))
  after
}

plot_entries <- list()

before_only <- before
before_only$is_condition_snp <- before_only$SNP %in% cond_snps
plot_entries[[length(plot_entries) + 1]] <- list(
  label = "00_original",
  title = paste0(analysis_name, " | ", locus_id, " | original GWAS"),
  data = build_plot_df(before_only, "p", "original")
)

if (!is.na(indiv_summary) && nzchar(indiv_summary) && file.exists(indiv_summary) && length(cond_snps) > 0) {
  indiv_tbl <- read_tsv(indiv_summary)
  indiv_tbl$snp <- as.character(indiv_tbl$snp)

  for (i in seq_along(cond_snps)) {
    snp_id <- cond_snps[i]
    row_idx <- which(indiv_tbl$snp == snp_id)
    if (length(row_idx) == 0) next

    cma_path <- indiv_tbl$cma_file[row_idx[1]]
    status <- indiv_tbl$status[row_idx[1]]
    if (!is.na(status) && status != "ok") next
    if (is.na(cma_path) || !nzchar(cma_path) || !file.exists(cma_path)) next

    merged <- load_after_file(cma_path)
    plot_entries[[length(plot_entries) + 1]] <- list(
      label = sprintf("%02d_single_%s", i, gsub("[^A-Za-z0-9._-]", "_", snp_id)),
      title = paste0(analysis_name, " | ", locus_id, " | conditioned on ", snp_id),
      data = build_plot_df(merged, "p_after", paste0("cond_", snp_id))
    )
  }
} else {
  iter_files <- list.files(cojo_dir, pattern = "^cojo\\.iter[0-9]+\\.cma\\.cojo$", full.names = TRUE)
  iter_ids <- integer(0)
  if (length(iter_files) > 0) {
    iter_ids <- as.integer(sub("^cojo\\.iter([0-9]+)\\.cma\\.cojo$", "\\1", basename(iter_files)))
    iter_files <- iter_files[order(iter_ids)]
  }

  for (i in seq_along(iter_files)) {
    merged <- load_after_file(iter_files[i])
    plot_entries[[length(plot_entries) + 1]] <- list(
      label = sprintf("%02d_iter_%s", i, iter_ids[i]),
      title = paste0(analysis_name, " | ", locus_id, " | iterative step ", iter_ids[i]),
      data = build_plot_df(merged, "p_after", paste0("iter", iter_ids[i]))
    )
  }
}

merged_final <- load_after_file(cojo_final)
plot_entries[[length(plot_entries) + 1]] <- list(
  label = sprintf("%02d_all", length(plot_entries)),
  title = paste0(analysis_name, " | ", locus_id, " | all conditioned SNPs"),
  data = build_plot_df(merged_final, "p_after", "all")
)
n_plots <- 0
for (entry in plot_entries) {
  plot_df <- entry$data
  if (nrow(plot_df) == 0) next

  y_vals <- -log10(plot_df$p)
  y_vals <- y_vals[is.finite(y_vals)]
  if (length(y_vals) == 0) next
  y_top <- max(y_vals) * 1.05
  if (!is.finite(y_top) || y_top <= 0) y_top <- 1

  gene_df <- data.frame()
  if (all(is.finite(plot_df$pos))) {
    gene_df <- fetch_gene_track(locus_chrom, min(plot_df$pos, na.rm = TRUE), max(plot_df$pos, na.rm = TRUE))
    if (nrow(gene_df) > 0) {
      track_step <- max(y_top * 0.08, 0.35)
      gene_df$y <- -track_step * gene_df$track
      gene_df$label_y <- gene_df$y - track_step * 0.28
      y_bottom <- min(gene_df$label_y) - track_step * 0.45
    } else {
      y_bottom <- 0
    }
  } else {
    y_bottom <- 0
  }

  p <- ggplot(plot_df, aes(x = pos, y = -log10(p))) +
    geom_point(size = 0.35, alpha = 0.7) +
    geom_hline(yintercept = -log10(5e-8), color = "red", linewidth = 0.5) +
    labs(
      title = entry$title,
      x = "Position (bp)",
      y = "-log10(p)"
    ) +
    coord_cartesian(ylim = c(y_bottom, y_top), clip = "off") +
    theme_bw(base_size = 10) +
    theme(plot.margin = margin(5.5, 5.5, 30, 5.5))

  if (nrow(gene_df) > 0) {
    p <- p +
      geom_segment(
        data = gene_df,
        aes(x = start, xend = end, y = y, yend = y),
        inherit.aes = FALSE,
        color = "#2f5d80",
        linewidth = 0.7
      ) +
      geom_segment(
        data = gene_df,
        aes(x = start, xend = start, y = y - abs(y_top) * 0.015, yend = y + abs(y_top) * 0.015),
        inherit.aes = FALSE,
        color = "#2f5d80",
        linewidth = 0.5
      ) +
      geom_segment(
        data = gene_df,
        aes(x = end, xend = end, y = y - abs(y_top) * 0.015, yend = y + abs(y_top) * 0.015),
        inherit.aes = FALSE,
        color = "#2f5d80",
        linewidth = 0.5
      ) +
      geom_text(
        data = gene_df,
        aes(x = mid, y = label_y, label = gene),
        inherit.aes = FALSE,
        color = "#2f5d80",
        size = 2.5,
        check_overlap = TRUE
      )
  }

  if (identical(entry$label, "00_original")) {
    highlight_df <- before_only[before_only$is_condition_snp & is.finite(before_only$p) & before_only$p > 0, c("SNP", "pos", "p")]
    if (nrow(highlight_df) > 0) {
      p <- p +
        geom_point(
          data = highlight_df,
          aes(x = pos, y = -log10(p)),
          inherit.aes = FALSE,
          color = "#d55e00",
          size = 1.8,
          alpha = 0.95
        ) +
        geom_text(
          data = highlight_df,
          aes(x = pos, y = -log10(p), label = SNP),
          inherit.aes = FALSE,
          color = "#d55e00",
          size = 2.7,
          vjust = -0.6,
          check_overlap = TRUE
        )
    }
  }

  out_png <- file.path(out_dir, paste0(entry$label, ".png"))
  out_pdf <- file.path(out_dir, paste0(entry$label, ".pdf"))

  ggsave(out_png, p, width = 11, height = 4.5, dpi = 180)
  ggsave(out_pdf, p, width = 11, height = 4.5)
  n_plots <- n_plots + 1
}

writeLines(
  paste0("ok\tplots_generated=", n_plots, "\tout_dir=", out_dir),
  con = out_done
)
cat("Generated", n_plots, "simple before/after iteration plots in", out_dir, "\n")
