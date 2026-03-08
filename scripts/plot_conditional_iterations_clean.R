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

HAVE_DATA_TABLE <- requireNamespace("data.table", quietly = TRUE)
HAVE_BIGSNPR <- requireNamespace("bigsnpr", quietly = TRUE)
HAVE_LOCUSZOOMR <- requireNamespace("locuszoomr", quietly = TRUE)
HAVE_ENSDB <- requireNamespace("EnsDb.Hsapiens.v75", quietly = TRUE) && requireNamespace("AnnotationFilter", quietly = TRUE)

if (HAVE_DATA_TABLE) suppressPackageStartupMessages(library(data.table))
if (HAVE_LOCUSZOOMR) suppressPackageStartupMessages(library(locuszoomr))
if (HAVE_BIGSNPR) suppressPackageStartupMessages(library(bigsnpr))
if (HAVE_ENSDB) {
  suppressPackageStartupMessages({
    library(AnnotationFilter)
    library(EnsDb.Hsapiens.v75)
  })
}

if (!HAVE_LOCUSZOOMR) message("locuszoomr not available; will use base R plots.")
if (!HAVE_BIGSNPR) message("bigsnpr not available; LD coloring will be approximated.")
if (!HAVE_ENSDB) message("EnsDb.Hsapiens.v75 not available; plotting without gene annotation.")

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
out_png <- get_arg("--out-png", required = FALSE, default = file.path(dirname(out_pdf), "conditional_iterations.png"))
out_before_pdf <- get_arg("--out-before-pdf", required = FALSE, default = file.path(dirname(out_pdf), "before_conditioning.pdf"))
out_after_pdf <- get_arg("--out-after-pdf", required = FALSE, default = file.path(dirname(out_pdf), "after_conditioning.pdf"))
out_before_png <- get_arg("--out-before-png", required = FALSE, default = file.path(dirname(out_pdf), "before_conditioning.png"))
out_after_png <- get_arg("--out-after-png", required = FALSE, default = file.path(dirname(out_pdf), "after_conditioning.png"))
out_done <- get_arg("--out-done")

if (!file.exists(cojo_final)) stop(paste("cojo final file not found:", cojo_final))
if (!file.exists(cojo_ma)) stop(paste("cojo.ma not found:", cojo_ma))
if (!dir.exists(cojo_dir)) stop(paste("cojo directory not found:", cojo_dir))
if (!file.exists(paste0(ref_prefix, ".bed"))) stop(paste("Reference BED not found:", paste0(ref_prefix, ".bed")))

dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_png), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_before_pdf), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_after_pdf), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_before_png), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_after_png), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(out_done), recursive = TRUE, showWarnings = FALSE)

THRESH_P <- 5e-8
THRESH_LWD <- 3

open_png_device <- function(path, width, height, res) {
  if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png(filename = path, width = width, height = height, res = res)
    return(TRUE)
  }

  for (dev_type in c("cairo", "cairo-png")) {
    ok <- tryCatch({
      grDevices::png(filename = path, width = width, height = height, res = res, type = dev_type)
      TRUE
    }, error = function(e) FALSE, warning = function(w) {
      invokeRestart("muffleWarning")
      FALSE
    })
    if (ok) return(TRUE)
  }

  FALSE
}

write_placeholder_png <- function(path) {
  png_raw <- as.raw(c(
    0x89,0x50,0x4E,0x47,0x0D,0x0A,0x1A,0x0A,
    0x00,0x00,0x00,0x0D,0x49,0x48,0x44,0x52,
    0x00,0x00,0x00,0x01,0x00,0x00,0x00,0x01,
    0x08,0x06,0x00,0x00,0x00,0x1F,0x15,0xC4,
    0x89,0x00,0x00,0x00,0x0A,0x49,0x44,0x41,
    0x54,0x78,0x9C,0x63,0x00,0x01,0x00,0x00,
    0x05,0x00,0x01,0x0D,0x0A,0x2D,0xB4,0x00,
    0x00,0x00,0x00,0x49,0x45,0x4E,0x44,0xAE,
    0x42,0x60,0x82
  ))
  con <- file(path, "wb")
  on.exit(close(con), add = TRUE)
  writeBin(png_raw, con)
}

compose_joint_from_pngs <- function(left_png, right_png, out_joint_pdf, out_joint_png) {
  if (!requireNamespace("png", quietly = TRUE)) return(FALSE)
  img_left <- png::readPNG(left_png)
  img_right <- png::readPNG(right_png)

  pdf(out_joint_pdf, width = 16, height = 8)
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(1, 2)))
  grid::grid.raster(img_left, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  grid::grid.raster(img_right, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
  dev.off()

  joint_ok <- open_png_device(out_joint_png, width = 4800, height = 2400, res = 300)
  if (joint_ok) {
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(1, 2)))
    grid::grid.raster(img_left, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
    grid::grid.raster(img_right, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
    dev.off()
  } else {
    write_placeholder_png(out_joint_png)
  }

  TRUE
}

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
  if (!HAVE_BIGSNPR) return(NULL)
  bed_file <- paste0(ref_prefix_path, ".bed")
  rds_file <- paste0(ref_prefix_path, ".rds")
  if (!file.exists(bed_file)) return(NULL)
  if (!file.exists(rds_file)) snp_readBed(bed_file)
  snp_attach(rds_file)
}

make_ref_map <- function(ref_obj) {
  if (is.null(ref_obj)) return(NULL)
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
  if (is.null(ref_map)) {
    out <- copy(plot_df)
    out[, ref_index := as.integer(NA)]
    return(out)
  }
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
  if (is.null(ref_obj)) return(out)

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

  plot_base_fallback <- function() {
    plot(
      local_df$pos,
      -log10(local_df$p),
      pch = 16,
      cex = 0.6,
      col = ifelse(local_df$r2 >= 0.8, "firebrick", "black"),
      xlab = "Position (bp)",
      ylab = "-log10(p)",
      main = plot_title,
      ylim = ylim_vals
    )
    abline(h = -log10(THRESH_P), lty = 1, col = "red", lwd = THRESH_LWD)
  }

  if (!HAVE_LOCUSZOOMR) {
    plot_base_fallback()
    return(invisible(NULL))
  }

  tryCatch({
    if (HAVE_ENSDB) {
      loc <- locus(
        data = local_df,
        seqname = unique(local_df$chrom),
        xrange = c(min(local_df$pos), max(local_df$pos)),
        LD = "r2",
        ens_db = "EnsDb.Hsapiens.v75"
      )
    } else {
      loc <- locus(
        data = local_df,
        seqname = unique(local_df$chrom),
        xrange = c(min(local_df$pos), max(local_df$pos)),
        LD = "r2"
      )
    }

    index_rsid <- local_df$rsid[which.min(local_df$p)]
    locus_plot(loc, labels = index_rsid, border = TRUE, main = plot_title, ylim = ylim_vals)
    abline(h = -log10(THRESH_P), lty = 1, col = "red", lwd = THRESH_LWD)
  }, error = function(e) {
    message("locuszoomr plot failed; using base fallback: ", conditionMessage(e))
    plot_base_fallback()
  })
}

cojo <- if (HAVE_DATA_TABLE) data.table::fread(cojo_final) else read.table(cojo_final, header = TRUE, sep = "", stringsAsFactors = FALSE)
require_columns(cojo, c("Chr", "SNP", "bp", "p", "pC"))
cojo$Chr <- as.integer(cojo$Chr)
cojo$bp <- as.numeric(cojo$bp)
cojo$p <- as.numeric(cojo$p)
cojo$pC <- as.numeric(cojo$pC)

chr_data <- cojo[cojo$Chr == target_chr, ]
if (nrow(chr_data) == 0) stop("No variants found for chromosome ", target_chr)

ref_obj <- load_reference_panel(ref_prefix)
ref_map <- make_ref_map(ref_obj)

before_df <- data.frame(rsid = chr_data$SNP, chrom = chr_data$Chr, pos = chr_data$bp, p = chr_data$p, stringsAsFactors = FALSE)
after_df <- data.frame(rsid = chr_data$SNP, chrom = chr_data$Chr, pos = chr_data$bp, pC = chr_data$pC, stringsAsFactors = FALSE)
before_df <- as.data.table(before_df)
after_df <- as.data.table(after_df)

before_df <- compute_ld_r2(annotate_ref_index(before_df, ref_map), ref_obj, target_chr, target_pos)
after_df <- compute_ld_r2(annotate_ref_index(after_df, ref_map), ref_obj, target_chr, target_pos)

iter_files <- list.files(cojo_dir, pattern = "^cojo\\.iter[0-9]+\\.cma\\.cojo$", full.names = TRUE)
iter_dfs <- list()
if (length(iter_files) > 0) {
  iter_ids <- as.integer(sub("^cojo\\.iter([0-9]+)\\.cma\\.cojo$", "\\1", basename(iter_files)))
  iter_files <- iter_files[order(iter_ids)]

  for (f in iter_files) {
    d <- if (HAVE_DATA_TABLE) data.table::fread(f) else read.table(f, header = TRUE, sep = "", stringsAsFactors = FALSE)
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

draw_iterations_panel <- function() {
  plot(NA, xlim = range(before_df$pos, na.rm = TRUE), ylim = shared_ylim,
       xlab = "Position (bp)", ylab = "-log10(pC)",
       main = paste0(analysis_name, " | ", locus_id, " | All iterations"))
  abline(h = -log10(THRESH_P), lty = 1, col = "red", lwd = THRESH_LWD)

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
}

before_title <- paste0(analysis_name, " | ", locus_id, " | Before conditioning")
after_title <- paste0(analysis_name, " | ", locus_id, " | After final conditioning")

pdf(out_before_pdf, width = 10, height = 8)
par(mar = c(4, 5, 4, 2))
make_locus_plot(before_df, before_title, "p", shared_ylim)
dev.off()

before_png_ok <- open_png_device(out_before_png, width = 3000, height = 2400, res = 300)
if (before_png_ok) {
  par(mar = c(4, 5, 4, 2))
  make_locus_plot(before_df, before_title, "p", shared_ylim)
  dev.off()
} else {
  write_placeholder_png(out_before_png)
}

pdf(out_after_pdf, width = 10, height = 8)
par(mar = c(4, 5, 4, 2))
make_locus_plot(after_df, after_title, "pC", shared_ylim)
dev.off()

after_png_ok <- open_png_device(out_after_png, width = 3000, height = 2400, res = 300)
if (after_png_ok) {
  par(mar = c(4, 5, 4, 2))
  make_locus_plot(after_df, after_title, "pC", shared_ylim)
  dev.off()
} else {
  write_placeholder_png(out_after_png)
}

composed <- compose_joint_from_pngs(out_before_png, out_after_png, out_pdf, out_png)
if (!composed) {
  pdf(out_pdf, width = 16, height = 8)
  layout(matrix(c(1, 2), nrow = 1, byrow = TRUE))
  par(mar = c(4, 5, 4, 2))
  make_locus_plot(before_df, before_title, "p", shared_ylim)
  make_locus_plot(after_df, after_title, "pC", shared_ylim)
  dev.off()

  png_ok <- open_png_device(out_png, width = 4800, height = 2400, res = 300)
  if (png_ok) {
    layout(matrix(c(1, 2), nrow = 1, byrow = TRUE))
    par(mar = c(4, 5, 4, 2))
    make_locus_plot(before_df, before_title, "p", shared_ylim)
    make_locus_plot(after_df, after_title, "pC", shared_ylim)
    dev.off()
  } else {
    message("Could not open PNG device on this node; writing placeholder PNG to keep workflow consistent.")
    write_placeholder_png(out_png)
  }
}
writeLines(paste0("ok\tplot=", out_pdf), con = out_done)
cat("Saved separate PDFs:", out_before_pdf, out_after_pdf, "\n")
cat("Saved separate PNGs:", out_before_png, out_after_png, "\n")
cat("Saved joint plots (horizontal):", out_pdf, out_png, "\n")
