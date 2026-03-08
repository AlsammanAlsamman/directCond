#!/usr/bin/env Rscript

################################################################################
# COJO Conditional Plot Generator (Before vs After)
# - Before: p
# - After : pC
# - Includes gene tracks via locuszoomr + EnsDb.Hsapiens.v75
################################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(locuszoomr)
  library(bigsnpr)
  library(EnsDb.Hsapiens.v75)
})

################################################################################
# CONFIG
################################################################################

COJO_FILE <- "cojo_chr21_43855067.cma.cojo"
REF_PREFIX <- "ref/g1000_eas_chr21"
OUTPUT_PREFIX <- "cojo_chr21_43855067"
STUDY_LABEL <- "Conditional analysis"
TARGET_CHR <- 21
TARGET_POS <- 43855067

################################################################################
# HELPERS
################################################################################

require_columns <- function(dt, cols) {
  miss <- setdiff(cols, names(dt))
  if (length(miss) > 0) {
    stop("Missing required columns: ", paste(miss, collapse = ", "))
  }
}

pick_first_col <- function(nms, candidates) {
  hit <- intersect(candidates, nms)
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

load_reference_panel <- function(ref_prefix) {
  bed_file <- paste0(ref_prefix, ".bed")
  rds_file <- paste0(ref_prefix, ".rds")

  if (!file.exists(bed_file)) {
    stop("Reference BED not found: ", bed_file)
  }

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
  if (length(id_hit) > 0) {
    out$ref_index[id_hit] <- ref_map$ref_index[id_match[id_hit]]
  }

  na_idx <- which(is.na(out$ref_index))
  if (length(na_idx) > 0) {
    pos_key <- paste0(out$chrom[na_idx], ":", out$pos[na_idx])
    ref_key <- paste0(ref_map$chr, ":", ref_map$bp)
    pos_match <- match(pos_key, ref_key)
    pos_hit <- which(!is.na(pos_match))
    if (length(pos_hit) > 0) {
      out$ref_index[na_idx[pos_hit]] <- ref_map$ref_index[pos_match[pos_hit]]
    }
  }

  out
}

compute_ld_r2 <- function(plot_df, ref_obj, target_chr, target_pos) {
  out <- copy(plot_df)
  out[, r2 := as.numeric(0)]

  valid_idx <- which(!is.na(out$ref_index))
  if (length(valid_idx) == 0) {
    warning("No overlap between COJO SNPs and reference panel.")
    return(out)
  }

  target_candidates <- which(out$chrom == target_chr & out$pos == target_pos & !is.na(out$ref_index))
  if (length(target_candidates) > 0) {
    target_ref_index <- out$ref_index[target_candidates[1]]
  } else {
    mapped <- out[!is.na(ref_index)]
    nearest <- which.min(abs(mapped$pos - target_pos))
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
    if (is.finite(r)) {
      r2_lookup[i] <- r^2
    }
  }

  out$r2[valid_idx] <- r2_lookup[match(out$ref_index[valid_idx], unique_ref)]
  out[ref_index == target_ref_index, r2 := 1]

  out
}

make_locus_plot <- function(plot_df, out_pdf, plot_title, target_pos = NULL, y_limits = NULL) {
  plot_df <- plot_df[is.finite(p) & p > 0]
  if (nrow(plot_df) == 0) {
    stop("No valid p-values available for plotting: ", plot_title)
  }

  seqname <- unique(plot_df$chrom)
  if (length(seqname) != 1) {
    stop("Plot data must contain one chromosome only.")
  }

  x_min <- min(plot_df$pos)
  x_max <- max(plot_df$pos)

  loc <- locus(
    data = plot_df,
    seqname = seqname,
    xrange = c(x_min, x_max),
    LD = "r2",
    ens_db = "EnsDb.Hsapiens.v75"
  )

  index_idx <- which.min(plot_df$p)
  index_rsid <- plot_df$rsid[index_idx]

  pdf(out_pdf, width = 12, height = 10)
  on.exit(dev.off(), add = TRUE)

  locus_plot(
    loc,
    labels = index_rsid,
    border = TRUE,
    main = plot_title,
    ylim = y_limits
  )

  gw_y <- -log10(5e-8)
  abline(h = gw_y, lty = 2, col = "red", lwd = 2.5)
  abline(h = -log10(5e-5), lty = 2, col = "darkgreen", lwd = 1)
  text(
    x = x_min + 0.02 * (x_max - x_min),
    y = gw_y + 0.2,
    labels = expression(5 %*% 10^-8),
    col = "red",
    cex = 0.9,
    pos = 4
  )
  if (!is.null(target_pos)) {
    abline(v = target_pos, lty = 3, col = "blue", lwd = 1)
  }
}

################################################################################
# MAIN
################################################################################

cat("Reading:", COJO_FILE, "\n")
cojo <- fread(COJO_FILE)
require_columns(cojo, c("Chr", "SNP", "bp", "p", "pC"))

cat("Loading reference panel:", REF_PREFIX, "\n")
ref_obj <- load_reference_panel(REF_PREFIX)
ref_map <- make_ref_map(ref_obj)

cojo[, Chr := as.integer(Chr)]
cojo[, bp := as.numeric(bp)]
cojo[, p := as.numeric(p)]
cojo[, pC := as.numeric(pC)]

chr_data <- cojo[Chr == TARGET_CHR]
if (nrow(chr_data) == 0) {
  stop("No variants found for chromosome ", TARGET_CHR, " in ", COJO_FILE)
}

region_start <- min(chr_data$bp, na.rm = TRUE)
region_end <- max(chr_data$bp, na.rm = TRUE)

cat("Variants on chr", TARGET_CHR, ": ", nrow(chr_data), "\n", sep = "")
cat("Region: ", region_start, "-", region_end, "\n", sep = "")

before_df <- chr_data[, .(
  rsid = SNP,
  chrom = Chr,
  pos = bp,
  p = p
)]
before_df <- annotate_ref_index(before_df, ref_map)
before_df <- compute_ld_r2(before_df, ref_obj, TARGET_CHR, TARGET_POS)

after_df <- chr_data[, .(
  rsid = SNP,
  chrom = Chr,
  pos = bp,
  p = pC
)]
after_df <- annotate_ref_index(after_df, ref_map)
after_df <- compute_ld_r2(after_df, ref_obj, TARGET_CHR, TARGET_POS)

cat("LD mapped in before plot: ", sum(before_df$r2 > 0), " SNPs\n", sep = "")
cat("LD mapped in after plot: ", sum(after_df$r2 > 0), " SNPs\n", sep = "")

before_title <- paste0(
  STUDY_LABEL,
  "\nBefore conditioning (p)",
  "\nChr ", TARGET_CHR, ":", format(region_start, big.mark = ","), "-", format(region_end, big.mark = ",")
)

after_title <- paste0(
  STUDY_LABEL,
  "\nAfter conditioning (pC)",
  "\nChr ", TARGET_CHR, ":", format(region_start, big.mark = ","), "-", format(region_end, big.mark = ",")
)

before_pdf <- paste0(OUTPUT_PREFIX, "_before_conditional_ggplot.pdf")
after_pdf <- paste0(OUTPUT_PREFIX, "_after_conditional_ggplot.pdf")

all_p <- c(before_df$p, after_df$p)
all_p <- all_p[is.finite(all_p) & all_p > 0]
shared_ymax <- max(-log10(all_p), na.rm = TRUE)
shared_ylim <- c(0, shared_ymax * 1.05)

make_locus_plot(before_df, before_pdf, before_title, target_pos = TARGET_POS, y_limits = shared_ylim)
make_locus_plot(after_df, after_pdf, after_title, target_pos = TARGET_POS, y_limits = shared_ylim)

cat("Saved:\n")
cat("  ", before_pdf, "\n", sep = "")
cat("  ", after_pdf, "\n", sep = "")
