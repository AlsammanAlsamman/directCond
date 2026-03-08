#!/usr/bin/env Rscript

ensure_pkg <- function(pkg, bioc = FALSE) {
  if (requireNamespace(pkg, quietly = TRUE)) return(invisible(TRUE))

  if (bioc) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager", repos = "https://cloud.r-project.org")
    }
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  } else {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(paste0("Failed to install required package: ", pkg))
  }
  invisible(TRUE)
}

ensure_pkg("AnnotationFilter", bioc = TRUE)
ensure_pkg("EnsDb.Hsapiens.v75", bioc = TRUE)
ensure_pkg("data.table")
ensure_pkg("bigsnpr")
ensure_pkg("locuszoomr")

suppressPackageStartupMessages({
  library(data.table)
  library(locuszoomr)
  library(bigsnpr)
  library(AnnotationFilter)
  library(EnsDb.Hsapiens.v75)
})
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
cojo_final <- get_arg("--cojo-final")
cojo_dir <- get_arg("--cojo-dir")
ref_prefix <- get_arg("--ref-prefix")
target_chr <- as.integer(get_arg("--target-chr"))
target_pos <- as.numeric(get_arg("--target-pos"))
analysis_name <- get_arg("--analysis", required = FALSE, default = "analysis")
locus_id <- get_arg("--locus-id", required = FALSE, default = "locus")
out_pdf <- get_arg("--out-pdf")
out_done <- get_arg("--out-done")
if (!file.exists(cojo_final)) stop(paste("cojo final file not found:", cojo_final))
if (!file.exists(cojo_ma)) stop(paste("cojo.ma not found:", cojo_ma))
if (!dir.exists(cojo_dir)) stop(paste("cojo directory not found:", cojo_dir))
if (!file.exists(paste0(ref_prefix, ".bed"))) stop(paste("Reference BED not found:", paste0(ref_prefix, ".bed")))

dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_done), recursive = TRUE, showWarnings = FALSE)

require_columns <- function(dt, cols) {
  miss <- setdiff(cols, names(dt))
  if (length(miss) > 0) stop("Missing required columns: ", paste(miss, collapse = ", "))
}

pick_first_col <- function(nms, candidates) {
  hit <- intersect(candidates, nms)
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

load_reference_panel <- function(ref_prefix_path) {
  bed_file <- paste0(ref_prefix_path, ".bed")
  rds_file <- paste0(ref_prefix_path, ".rds")
  if (!file.exists(rds_file)) snp_readBed(bed_file)
  snp_attach(rds_file)
}

make_ref_map <- function(ref_obj) {
  map <- as.data.table(ref_obj$map)
  nms <- names(map)

  chr_col <- pick_first_col(nms, c("chromosome", "chr", "chrom"))
  pos_col <- pick_first_col(nms, c("physical.pos", "bp", "pos", "position"))
  id_col <- pick_first_col(nms, c("marker.ID", "rsid", "snp", "id", "ID"))

  if (is.na(chr_col) || is.na(pos_col) || is.na(id_col)) {
    stop("Could not detect chr/pos/id columns in bigsnpr map.")
  }

  data.table(
  ref_index = seq_len(nrow(map)),
  chr = as.integer(map[[chr_col]]),
  bp = as.numeric(map[[pos_col]]),
  ref_id = as.character(map[[id_col]])
  )
}

annotate_ref_index <- function(plot_df, ref_map) {
  out <- copy(plot_df)
  out[, ref_index := as.integer(NA)]

  id_match <- match(out$rsid, ref_map$ref_id)
  id_hit <- which(!is.na(id_match))
  if (length(id_hit) > 0) out$ref_index[id_hit] <- ref_map$ref_index[id_match[id_hit]]

  na_idx <- which(is.na(out$ref_index))
  if (length(na_idx) > 0) {
    pos_key <- paste0(out$chrom[na_idx], ":", out$pos[na_idx])
    ref_key <- paste0(ref_map$chr, ":", ref_map$bp)
    pos_match <- match(pos_key, ref_key)
    pos_hit <- which(!is.na(pos_match))
    if (length(pos_hit) > 0) out$ref_index[na_idx[pos_hit]] <- ref_map$ref_index[pos_match[pos_hit]]
  }

  out
}

compute_ld_r2 <- function(plot_df, ref_obj, target_chr_num, target_pos_num) {
  out <- copy(plot_df)
  out[, r2 := as.numeric(0)]

  valid_idx <- which(!is.na(out$ref_index))
  if (length(valid_idx) == 0) return(out)

  target_candidates <- which(out$chrom == target_chr_num & out$pos == target_pos_num & !is.na(out$ref_index))
  if (length(target_candidates) > 0) {
    target_ref_index <- out$ref_index[target_candidates[1]]
  } else {
    mapped <- out[!is.na(ref_index)]
    nearest <- which.min(abs(mapped$pos - target_pos_num))
    target_ref_index <- mapped$ref_index[nearest]
  }

  G <- ref_obj$genotypes
  target_g <- G[, target_ref_index]

  unique_ref <- unique(out$ref_index[valid_idx])
  r2_lookup <- rep(0, length(unique_ref))

  for (i in seq_along(unique_ref)) {
    snp_g <- G[, unique_ref[i]]
    r <- suppressWarnings(cor(target_g[], snp_g[], use = "pairwise.complete.obs"))
    if (is.finite(r)) r2_lookup[i] <- r^2
  }

  out$r2[valid_idx] <- r2_lookup[match(out$ref_index[valid_idx], unique_ref)]
  out[ref_index == target_ref_index, r2 := 1]
  out
}

make_locus_plot <- function(plot_df, plot_title, y_col, ylim_vals) {
  local_df <- copy(plot_df)
  local_df[, p := as.numeric(get(y_col))]
  local_df <- local_df[is.finite(p) & p > 0]
  if (nrow(local_df) == 0) stop("No valid p-values for plotting")

  loc <- locus(
    data = local_df,
    seqname = unique(local_df$chrom),
    xrange = c(min(local_df$pos), max(local_df$pos)),
    LD = "r2",
    ens_db = "EnsDb.Hsapiens.v75"
  )

  index_rsid <- local_df$rsid[which.min(local_df$p)]
  locus_plot(loc, labels = index_rsid, border = TRUE, main = plot_title, ylim = ylim_vals)
  abline(h = -log10(5e-8), lty = 2, col = "red", lwd = 2)
}

cojo <- data.table::fread(cojo_final)
require_columns(cojo, c("Chr", "SNP", "bp", "p", "pC"))
cojo[, `:=`(Chr = as.integer(Chr), bp = as.numeric(bp), p = as.numeric(p), pC = as.numeric(pC))]

chr_data <- cojo[Chr == target_chr]
if (nrow(chr_data) == 0) stop("No variants found for chromosome ", target_chr)

ref_obj <- load_reference_panel(ref_prefix)
ref_map <- make_ref_map(ref_obj)

before_df <- chr_data[, .(rsid = SNP, chrom = Chr, pos = bp, p = p)]
after_df <- chr_data[, .(rsid = SNP, chrom = Chr, pos = bp, pC = pC)]

before_df <- compute_ld_r2(annotate_ref_index(before_df, ref_map), ref_obj, target_chr, target_pos)
after_df <- compute_ld_r2(annotate_ref_index(after_df, ref_map), ref_obj, target_chr, target_pos)

iter_files <- list.files(cojo_dir, pattern = "^cojo\\.iter[0-9]+\\.cma\\.cojo$", full.names = TRUE)
iter_dfs <- list()
if (length(iter_files) > 0) {
  iter_ids <- as.integer(sub("^cojo\\.iter([0-9]+)\\.cma\\.cojo$", "\\1", basename(iter_files)))
  iter_files <- iter_files[order(iter_ids)]

  for (f in iter_files) {
    d <- data.table::fread(f)
    if (!all(c("bp", "pC") %in% names(d))) next
    d <- d[, .(bp = as.numeric(bp), pC = as.numeric(pC))]
    d <- d[is.finite(bp) & is.finite(pC) & pC > 0]
    if (nrow(d) == 0) next
    d[, iter := sub("\\.cma\\.cojo$", "", basename(f))]
    iter_dfs[[length(iter_dfs) + 1]] <- d
  }
}

all_p <- c(before_df$p, after_df$pC)
all_p <- all_p[is.finite(all_p) & all_p > 0]
shared_ylim <- c(0, max(-log10(all_p), na.rm = TRUE) * 1.05)

pdf(out_pdf, width = 12, height = 15)
par(mfrow = c(3, 1), mar = c(4, 5, 4, 2))

make_locus_plot(before_df, paste0(analysis_name, " | ", locus_id, " | Before conditioning"), "p", shared_ylim)
make_locus_plot(after_df, paste0(analysis_name, " | ", locus_id, " | After final conditioning"), "pC", shared_ylim)

plot(NA, xlim = range(before_df$pos, na.rm = TRUE), ylim = shared_ylim,
     xlab = "Position (bp)", ylab = "-log10(pC)",
     main = paste0(analysis_name, " | ", locus_id, " | All iterations"))
abline(h = -log10(5e-8), lty = 2, col = "red", lwd = 2)

if (length(iter_dfs) > 0) {
  iter_all <- data.table::rbindlist(iter_dfs)
  labels <- unique(iter_all$iter)
  cols <- setNames(rainbow(length(labels)), labels)
  for (lab in labels) {
    d <- iter_all[iter == lab]
    points(d$bp, -log10(d$pC), pch = 16, cex = 0.5, col = cols[[lab]])
  }
  legend("topright", legend = labels, pch = 16, col = unname(cols[labels]), bty = "n", cex = 0.75)
} else {
  points(after_df$pos, -log10(after_df$pC), pch = 16, cex = 0.5, col = "firebrick")
  legend("topright", legend = "final", pch = 16, col = "firebrick", bty = "n", cex = 0.75)
}

dev.off()
writeLines(paste0("ok\tplot=", out_pdf), con = out_done)
cat("Saved plot:", out_pdf, "\n")
#!/usr/bin/env Rscript

ensure_pkg <- function(pkg, bioc = FALSE) {
  if (requireNamespace(pkg, quietly = TRUE)) return(invisible(TRUE))

  if (bioc) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager", repos = "https://cloud.r-project.org")
    }
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  } else {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }

  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(paste0("Failed to install required package: ", pkg))
  }
  invisible(TRUE)
}

ensure_pkg("AnnotationFilter", bioc = TRUE)
ensure_pkg("EnsDb.Hsapiens.v75", bioc = TRUE)
ensure_pkg("data.table")
ensure_pkg("bigsnpr")
ensure_pkg("locuszoomr")

suppressPackageStartupMessages({
  library(data.table)
  library(locuszoomr)
  library(bigsnpr)
  library(AnnotationFilter)
  library(EnsDb.Hsapiens.v75)
})

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
cojo_final <- get_arg("--cojo-final")
cojo_dir <- get_arg("--cojo-dir")
ref_prefix <- get_arg("--ref-prefix")
target_chr <- as.integer(get_arg("--target-chr"))
target_pos <- as.numeric(get_arg("--target-pos"))
analysis_name <- get_arg("--analysis", required = FALSE, default = "analysis")
locus_id <- get_arg("--locus-id", required = FALSE, default = "locus")
out_pdf <- get_arg("--out-pdf")
out_done <- get_arg("--out-done")

if (!file.exists(cojo_final)) stop(paste("cojo final file not found:", cojo_final))
if (!file.exists(cojo_ma)) stop(paste("cojo.ma not found:", cojo_ma))
if (!dir.exists(cojo_dir)) stop(paste("cojo directory not found:", cojo_dir))
if (!file.exists(paste0(ref_prefix, ".bed"))) stop(paste("Reference BED not found:", paste0(ref_prefix, ".bed")))

dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_done), recursive = TRUE, showWarnings = FALSE)

require_columns <- function(dt, cols) {
  miss <- setdiff(cols, names(dt))
  if (length(miss) > 0) stop("Missing required columns: ", paste(miss, collapse = ", "))
}

pick_first_col <- function(nms, candidates) {
  hit <- intersect(candidates, nms)
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

load_reference_panel <- function(ref_prefix_path) {
  bed_file <- paste0(ref_prefix_path, ".bed")
  rds_file <- paste0(ref_prefix_path, ".rds")
  if (!file.exists(rds_file)) {
    cat("Converting PLINK to bigsnpr backing files...\n")
    snp_readBed(bed_file)
  }
  snp_attach(rds_file)
}

make_ref_map <- function(ref_obj) {
  map <- as.data.table(ref_obj$map)
  nms <- names(map)

  chr_col <- pick_first_col(nms, c("chromosome", "chr", "chrom"))
  pos_col <- pick_first_col(nms, c("physical.pos", "bp", "pos", "position"))
  id_col <- pick_first_col(nms, c("marker.ID", "rsid", "snp", "id", "ID"))

  if (is.na(chr_col) || is.na(pos_col) || is.na(id_col)) stop("Could not detect chr/pos/id columns in bigsnpr map.")

  data.table(
    ref_index = seq_len(nrow(map)),
    chr = as.integer(map[[chr_col]]),
    bp = as.numeric(map[[pos_col]]),
    ref_id = as.character(map[[id_col]])
  )
}

annotate_ref_index <- function(plot_df, ref_map) {
  out <- copy(plot_df)
  out[, ref_index := as.integer(NA)]

  id_match <- match(out$rsid, ref_map$ref_id)
  id_hit <- which(!is.na(id_match))
  if (length(id_hit) > 0) out$ref_index[id_hit] <- ref_map$ref_index[id_match[id_hit]]

  na_idx <- which(is.na(out$ref_index))
  if (length(na_idx) > 0) {
    pos_key <- paste0(out$chrom[na_idx], ":", out$pos[na_idx])
    ref_key <- paste0(ref_map$chr, ":", ref_map$bp)
    pos_match <- match(pos_key, ref_key)
    pos_hit <- which(!is.na(pos_match))
    if (length(pos_hit) > 0) out$ref_index[na_idx[pos_hit]] <- ref_map$ref_index[pos_match[pos_hit]]
  }

  out
}

compute_ld_r2 <- function(plot_df, ref_obj, target_chr_num, target_pos_num) {
  out <- copy(plot_df)
  out[, r2 := as.numeric(0)]

  valid_idx <- which(!is.na(out$ref_index))
  if (length(valid_idx) == 0) return(out)

  target_candidates <- which(out$chrom == target_chr_num & out$pos == target_pos_num & !is.na(out$ref_index))
  if (length(target_candidates) > 0) {
    target_ref_index <- out$ref_index[target_candidates[1]]
  } else {
    mapped <- out[!is.na(ref_index)]
    nearest <- which.min(abs(mapped$pos - target_pos_num))
    target_ref_index <- mapped$ref_index[nearest]
    message("Target SNP not directly mapped; using nearest mapped SNP at bp=", mapped$pos[nearest])
  }

  G <- ref_obj$genotypes
  target_g <- G[, target_ref_index]

  unique_ref <- unique(out$ref_index[valid_idx])
  r2_lookup <- rep(0, length(unique_ref))

  for (i in seq_along(unique_ref)) {
    ref_i <- unique_ref[i]
    snp_g <- G[, ref_i]
    r <- suppressWarnings(cor(target_g[], snp_g[], use = "pairwise.complete.obs"))
    if (is.finite(r)) r2_lookup[i] <- r^2
  }

  out$r2[valid_idx] <- r2_lookup[match(out$ref_index[valid_idx], unique_ref)]
  out[ref_index == target_ref_index, r2 := 1]
  out
}

make_locus_plot <- function(plot_df, plot_title, y_col = "p", ylim = NULL) {
  local_df <- copy(plot_df)
  local_df[, p := as.numeric(get(y_col))]
  local_df <- local_df[is.finite(p) & p > 0]
  if (nrow(local_df) == 0) stop("No valid p-values for plotting")

  seqname <- unique(local_df$chrom)
  if (length(seqname) != 1) stop("Plot data must contain one chromosome only")

  x_min <- min(local_df$pos)
  x_max <- max(local_df$pos)

  loc <- locus(
    data = local_df,
    seqname = seqname,
    xrange = c(x_min, x_max),
    LD = "r2",
    ens_db = "EnsDb.Hsapiens.v75"
  )

  index_idx <- which.min(local_df$p)
  index_rsid <- local_df$rsid[index_idx]

  locus_plot(loc, labels = index_rsid, border = TRUE, main = plot_title, ylim = ylim)
  abline(h = -log10(5e-8), lty = 2, col = "red", lwd = 2)
}

cat("Reading final COJO:", cojo_final, "\n")
cojo <- data.table::fread(cojo_final)
require_columns(cojo, c("Chr", "SNP", "bp", "p", "pC"))

cojo[, Chr := as.integer(Chr)]
cojo[, bp := as.numeric(bp)]
cojo[, p := as.numeric(p)]
cojo[, pC := as.numeric(pC)]

chr_data <- cojo[Chr == target_chr]
if (nrow(chr_data) == 0) stop("No variants found for chromosome ", target_chr)

cat("Loading reference panel:", ref_prefix, "\n")
ref_obj <- load_reference_panel(ref_prefix)
ref_map <- make_ref_map(ref_obj)

before_df <- chr_data[, .(rsid = SNP, chrom = Chr, pos = bp, p = p)]
before_df <- annotate_ref_index(before_df, ref_map)
before_df <- compute_ld_r2(before_df, ref_obj, target_chr, target_pos)

after_df <- chr_data[, .(rsid = SNP, chrom = Chr, pos = bp, pC = pC)]
after_df <- annotate_ref_index(after_df, ref_map)
after_df <- compute_ld_r2(after_df, ref_obj, target_chr, target_pos)

iter_files <- list.files(cojo_dir, pattern = "^cojo\\.iter[0-9]+\\.cma\\.cojo$", full.names = TRUE)
iter_dfs <- list()
if (length(iter_files) > 0) {
  iter_ids <- as.integer(sub("^cojo\\.iter([0-9]+)\\.cma\\.cojo$", "\\1", basename(iter_files)))
  iter_files <- iter_files[order(iter_ids)]

  for (f in iter_files) {
    d <- data.table::fread(f)
    if (!all(c("bp", "pC") %in% names(d))) next
    d <- d[, .(bp = as.numeric(bp), pC = as.numeric(pC))]
    d <- d[is.finite(bp) & is.finite(pC) & pC > 0]
    if (nrow(d) == 0) next
    d[, iter := sub("\\.cma\\.cojo$", "", basename(f))]
    iter_dfs[[length(iter_dfs) + 1]] <- d
  }
}

all_p <- c(before_df$p, after_df$pC)
all_p <- all_p[is.finite(all_p) & all_p > 0]
shared_ymax <- max(-log10(all_p), na.rm = TRUE)
shared_ylim <- c(0, shared_ymax * 1.05)

pdf(out_pdf, width = 12, height = 15)
par(mfrow = c(3, 1), mar = c(4, 5, 4, 2))

before_title <- paste0(analysis_name, " | ", locus_id, " | Before conditioning")
after_title <- paste0(analysis_name, " | ", locus_id, " | After final conditioning")
make_locus_plot(before_df, before_title, y_col = "p", ylim = shared_ylim)
make_locus_plot(after_df, after_title, y_col = "pC", ylim = shared_ylim)

plot(NA, xlim = range(before_df$pos, na.rm = TRUE), ylim = shared_ylim,
     xlab = "Position (bp)", ylab = "-log10(pC)",
     main = paste0(analysis_name, " | ", locus_id, " | All iterations"))
abline(h = -log10(5e-8), lty = 2, col = "red", lwd = 2)

if (length(iter_dfs) > 0) {
  iter_all <- data.table::rbindlist(iter_dfs)
  labels <- unique(iter_all$iter)
  cols <- setNames(rainbow(length(labels)), labels)
  for (lab in labels) {
    d <- iter_all[iter == lab]
    points(d$bp, -log10(d$pC), pch = 16, cex = 0.5, col = cols[[lab]])
  }
  legend("topright", legend = labels, pch = 16, col = unname(cols[labels]), bty = "n", cex = 0.75)
} else {
  points(after_df$pos, -log10(after_df$pC), pch = 16, cex = 0.5, col = "firebrick")
  legend("topright", legend = "final", pch = 16, col = "firebrick", bty = "n", cex = 0.75)
}

dev.off()

writeLines(paste0("ok\tplot=", out_pdf), con = out_done)
cat("Saved plot:", out_pdf, "\n")
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(locuszoomr)
  library(bigsnpr)
  library(EnsDb.Hsapiens.v75)
})

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
cojo_final <- get_arg("--cojo-final")
cojo_dir <- get_arg("--cojo-dir")
ref_prefix <- get_arg("--ref-prefix")
target_chr <- as.integer(get_arg("--target-chr"))
target_pos <- as.numeric(get_arg("--target-pos"))
analysis_name <- get_arg("--analysis", required = FALSE, default = "analysis")
locus_id <- get_arg("--locus-id", required = FALSE, default = "locus")
out_pdf <- get_arg("--out-pdf")
out_done <- get_arg("--out-done")

if (!file.exists(cojo_final)) stop(paste("cojo final file not found:", cojo_final))
if (!dir.exists(cojo_dir)) stop(paste("cojo directory not found:", cojo_dir))
if (!file.exists(paste0(ref_prefix, ".bed"))) stop(paste("Reference BED not found:", paste0(ref_prefix, ".bed")))
if (!file.exists(cojo_ma)) stop(paste("cojo.ma not found:", cojo_ma))

dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_done), recursive = TRUE, showWarnings = FALSE)

require_columns <- function(dt, cols) {
  miss <- setdiff(cols, names(dt))
  if (length(miss) > 0) stop("Missing required columns: ", paste(miss, collapse = ", "))
}

pick_first_col <- function(nms, candidates) {
  hit <- intersect(candidates, nms)
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

load_reference_panel <- function(ref_prefix_path) {
  bed_file <- paste0(ref_prefix_path, ".bed")
  rds_file <- paste0(ref_prefix_path, ".rds")
  if (!file.exists(rds_file)) {
    cat("Converting PLINK to bigsnpr backing files...\n")
    snp_readBed(bed_file)
  }
  snp_attach(rds_file)
}

make_ref_map <- function(ref_obj) {
  map <- as.data.table(ref_obj$map)
  nms <- names(map)

  chr_col <- pick_first_col(nms, c("chromosome", "chr", "chrom"))
  pos_col <- pick_first_col(nms, c("physical.pos", "bp", "pos", "position"))
  id_col <- pick_first_col(nms, c("marker.ID", "rsid", "snp", "id", "ID"))

  if (is.na(chr_col) || is.na(pos_col) || is.na(id_col)) {
    stop("Could not detect chr/pos/id columns in bigsnpr map.")
  }

  data.table(
    ref_index = seq_len(nrow(map)),
    chr = as.integer(map[[chr_col]]),
    bp = as.numeric(map[[pos_col]]),
    ref_id = as.character(map[[id_col]])
  )
}

annotate_ref_index <- function(plot_df, ref_map) {
  out <- copy(plot_df)
  out[, ref_index := as.integer(NA)]

  id_match <- match(out$rsid, ref_map$ref_id)
  id_hit <- which(!is.na(id_match))
  if (length(id_hit) > 0) out$ref_index[id_hit] <- ref_map$ref_index[id_match[id_hit]]

  na_idx <- which(is.na(out$ref_index))
  if (length(na_idx) > 0) {
    pos_key <- paste0(out$chrom[na_idx], ":", out$pos[na_idx])
    ref_key <- paste0(ref_map$chr, ":", ref_map$bp)
    pos_match <- match(pos_key, ref_key)
    pos_hit <- which(!is.na(pos_match))
    if (length(pos_hit) > 0) out$ref_index[na_idx[pos_hit]] <- ref_map$ref_index[pos_match[pos_hit]]
  }

  out
}

compute_ld_r2 <- function(plot_df, ref_obj, target_chr_num, target_pos_num) {
  out <- copy(plot_df)
  out[, r2 := as.numeric(0)]

  valid_idx <- which(!is.na(out$ref_index))
  if (length(valid_idx) == 0) return(out)

  target_candidates <- which(out$chrom == target_chr_num & out$pos == target_pos_num & !is.na(out$ref_index))
  if (length(target_candidates) > 0) {
    target_ref_index <- out$ref_index[target_candidates[1]]
  } else {
    mapped <- out[!is.na(ref_index)]
    nearest <- which.min(abs(mapped$pos - target_pos_num))
    target_ref_index <- mapped$ref_index[nearest]
    message("Target SNP not directly mapped; using nearest mapped SNP at bp=", mapped$pos[nearest])
  }

  G <- ref_obj$genotypes
  target_g <- G[, target_ref_index]

  unique_ref <- unique(out$ref_index[valid_idx])
  r2_lookup <- rep(0, length(unique_ref))

  for (i in seq_along(unique_ref)) {
    ref_i <- unique_ref[i]
    snp_g <- G[, ref_i]
    r <- suppressWarnings(cor(target_g[], snp_g[], use = "pairwise.complete.obs"))
    if (is.finite(r)) r2_lookup[i] <- r^2
  }

  out$r2[valid_idx] <- r2_lookup[match(out$ref_index[valid_idx], unique_ref)]
  out[ref_index == target_ref_index, r2 := 1]
  out
}

make_locus_plot <- function(plot_df, plot_title, y_col = "p", ylim = NULL) {
  local_df <- copy(plot_df)
  local_df[, p := as.numeric(get(y_col))]
  local_df <- local_df[is.finite(p) & p > 0]
  if (nrow(local_df) == 0) stop("No valid p-values for plotting")

  seqname <- unique(local_df$chrom)
  if (length(seqname) != 1) stop("Plot data must contain one chromosome only")

  x_min <- min(local_df$pos)
  x_max <- max(local_df$pos)

  loc <- locus(
    data = local_df,
    seqname = seqname,
    xrange = c(x_min, x_max),
    LD = "r2",
    ens_db = "EnsDb.Hsapiens.v75"
  )

  index_idx <- which.min(local_df$p)
  index_rsid <- local_df$rsid[index_idx]

  locus_plot(
    loc,
    labels = index_rsid,
    border = TRUE,
    main = plot_title,
    ylim = ylim
  )

  abline(h = -log10(5e-8), lty = 2, col = "red", lwd = 2)
}

cat("Reading final COJO:", cojo_final, "\n")
cojo <- data.table::fread(cojo_final)
require_columns(cojo, c("Chr", "SNP", "bp", "p", "pC"))

cojo[, Chr := as.integer(Chr)]
cojo[, bp := as.numeric(bp)]
cojo[, p := as.numeric(p)]
cojo[, pC := as.numeric(pC)]

chr_data <- cojo[Chr == target_chr]
if (nrow(chr_data) == 0) stop("No variants found for chromosome ", target_chr)

cat("Loading reference panel:", ref_prefix, "\n")
ref_obj <- load_reference_panel(ref_prefix)
ref_map <- make_ref_map(ref_obj)

before_df <- chr_data[, .(rsid = SNP, chrom = Chr, pos = bp, p = p)]
before_df <- annotate_ref_index(before_df, ref_map)
before_df <- compute_ld_r2(before_df, ref_obj, target_chr, target_pos)

after_df <- chr_data[, .(rsid = SNP, chrom = Chr, pos = bp, pC = pC)]
after_df <- annotate_ref_index(after_df, ref_map)
after_df <- compute_ld_r2(after_df, ref_obj, target_chr, target_pos)

iter_files <- list.files(cojo_dir, pattern = "^cojo\\.iter[0-9]+\\.cma\\.cojo$", full.names = TRUE)
iter_dfs <- list()
if (length(iter_files) > 0) {
  iter_ids <- as.integer(sub("^cojo\\.iter([0-9]+)\\.cma\\.cojo$", "\\1", basename(iter_files)))
  iter_files <- iter_files[order(iter_ids)]

  for (f in iter_files) {
    d <- data.table::fread(f)
    if (!all(c("bp", "pC") %in% names(d))) next
    d <- d[, .(bp = as.numeric(bp), pC = as.numeric(pC))]
    d <- d[is.finite(bp) & is.finite(pC) & pC > 0]
    if (nrow(d) == 0) next
    d[, iter := sub("\\.cma\\.cojo$", "", basename(f))]
    iter_dfs[[length(iter_dfs) + 1]] <- d
  }
}

all_p <- c(before_df$p, after_df$pC)
all_p <- all_p[is.finite(all_p) & all_p > 0]
shared_ymax <- max(-log10(all_p), na.rm = TRUE)
shared_ylim <- c(0, shared_ymax * 1.05)

pdf(out_pdf, width = 12, height = 15)
par(mfrow = c(3, 1), mar = c(4, 5, 4, 2))

before_title <- paste0(analysis_name, " | ", locus_id, " | Before conditioning")
after_title <- paste0(analysis_name, " | ", locus_id, " | After final conditioning")
make_locus_plot(before_df, before_title, y_col = "p", ylim = shared_ylim)
make_locus_plot(after_df, after_title, y_col = "pC", ylim = shared_ylim)

plot(
  NA,
  xlim = range(before_df$pos, na.rm = TRUE),
  ylim = shared_ylim,
  xlab = "Position (bp)",
  ylab = "-log10(pC)",
  main = paste0(analysis_name, " | ", locus_id, " | All iterations")
)
abline(h = -log10(5e-8), lty = 2, col = "red", lwd = 2)

if (length(iter_dfs) > 0) {
  iter_all <- data.table::rbindlist(iter_dfs)
  labels <- unique(iter_all$iter)
  cols <- setNames(rainbow(length(labels)), labels)
  for (lab in labels) {
    d <- iter_all[iter == lab]
    points(d$bp, -log10(d$pC), pch = 16, cex = 0.5, col = cols[[lab]])
  }
  legend("topright", legend = labels, pch = 16, col = unname(cols[labels]), bty = "n", cex = 0.75)
} else {
  points(after_df$pos, -log10(after_df$pC), pch = 16, cex = 0.5, col = "firebrick")
  legend("topright", legend = "final", pch = 16, col = "firebrick", bty = "n", cex = 0.75)
}

dev.off()

writeLines(paste0("ok\tplot=", out_pdf), con = out_done)
cat("Saved plot:", out_pdf, "\n")
