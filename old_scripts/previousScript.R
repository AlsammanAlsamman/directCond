################################################################################
# LocusZoom Plot Generator for Multiple Loci
# Generates locus plots with LD information from reference panel
################################################################################

# Load required libraries
library(locuszoomr)
library(bigsnpr)
library(data.table)
library(EnsDb.Hsapiens.v75)

################################################################################
# STUDY CONFIGURATION - EDIT THIS SECTION
################################################################################

STUDY_NAME <- "Hispanic GWAS Study"
N_CASES <- 3018      # Number of cases
N_CONTROLS <- 5501   # Number of controls
REF_PANEL_NAME <- "1000 Genomes AMR"
REF_PANEL_PATH <- "ref_panel/g1000_amr"  # Path without .bed extension

# Input files
GWAS_FILE <- "Hisp.txt"       # GWAS summary statistics (rsid, chr, pos, p)
LOCI_FILE <- "Loci.txt"       # Loci to plot (Locus, Chr, Start, End)

# Output directory
OUTPUT_DIR <- "locus_plots"

# Region expansion settings
MIN_REGION_SIZE <- 500000  # Minimum region size in bp (0.5 Mbp)
AUTO_EXPAND_REGIONS <- TRUE  # Automatically expand small regions

# Set working directory
setwd("/s/nath-lab/alsamman/____MyCodes____/locusZoom")

################################################################################
# FUNCTIONS
################################################################################

#' Load or create reference panel in bigsnpr format
#'
#' @param ref_path Path to PLINK files (without extension)
#' @return bigsnpr reference panel object
load_reference_panel <- function(ref_path) {
  rds_file <- paste0(ref_path, ".rds")
  bed_file <- paste0(ref_path, ".bed")
  
  if (!file.exists(rds_file)) {
    cat("Converting PLINK files to bigsnpr format...\n")
    cat("This may take several minutes for large reference panels...\n")
    snp_readBed(bed_file)
  }
  
  cat("Loading reference panel from:", rds_file, "\n")
  ref_panel <- snp_attach(rds_file)
  cat("Reference panel loaded:", nrow(ref_panel$map), "SNPs,", 
      nrow(ref_panel$fam), "samples\n")
  
  return(ref_panel)
}

#' Calculate LD (r2) for a locus
#'
#' @param gwas_subset GWAS data for the locus
#' @param ref_panel bigsnpr reference panel object
#' @param target_chr Chromosome number
#' @param original_start Original region start (for finding index SNP)
#' @param original_end Original region end (for finding index SNP)
#' @return GWAS data with r2 column added
calculate_ld <- function(gwas_subset, ref_panel, target_chr, 
                         original_start = NULL, original_end = NULL) {
  
  cat("\n  Calculating LD for chr", target_chr, "...\n")
  
  # Get genotype matrix and map
  G <- ref_panel$genotypes
  map <- ref_panel$map
  
  # Filter reference panel to target chromosome
  chr_idx <- which(map$chromosome == target_chr)
  if (length(chr_idx) == 0) {
    warning("No SNPs found in reference panel for chromosome ", target_chr)
    gwas_subset$r2 <- 0
    return(gwas_subset)
  }
  
  map_chr <- map[chr_idx, ]
  G_chr <- G[, chr_idx]
  
  # Match positions between GWAS and reference panel
  matched_indices <- match(gwas_subset$pos, map_chr$physical.pos)
  gwas_subset$ref_idx <- matched_indices
  
  # Get indices of SNPs present in both datasets
  valid_snps <- which(!is.na(gwas_subset$ref_idx))
  
  cat("  SNPs in locus:", nrow(gwas_subset), "\n")
  cat("  Overlapping with ref panel:", length(valid_snps), "\n")
  
  # Initialize r2 column
  gwas_subset$r2 <- NA
  
  if (length(valid_snps) > 0) {
    # Find index SNP from original region only (if specified)
    if (!is.null(original_start) && !is.null(original_end)) {
      # Only consider SNPs from original region for index SNP
      original_region_snps <- which(gwas_subset$pos >= original_start & 
                                    gwas_subset$pos <= original_end & 
                                    !is.na(gwas_subset$ref_idx))
      if (length(original_region_snps) > 0) {
        index_idx <- original_region_snps[which.min(gwas_subset$p[original_region_snps])]
        cat("  Index SNP selected from original region\n")
      } else {
        # Fallback: use all valid SNPs if no SNPs in original region
        index_idx <- valid_snps[which.min(gwas_subset$p[valid_snps])]
        cat("  WARNING: No valid SNPs in original region, using expanded region\n")
      }
    } else {
      # Find index SNP (lowest p-value among all valid SNPs)
      index_idx <- valid_snps[which.min(gwas_subset$p[valid_snps])]
    }
    index_pos <- gwas_subset$pos[index_idx]
    index_rsid <- gwas_subset$rsid[index_idx]
    
    cat("  Index SNP:", index_rsid, "at position", index_pos, 
        "with p =", gwas_subset$p[index_idx], "\n")
    
    # Get genotype for index SNP
    index_ref_col <- which(map_chr$physical.pos == index_pos)
    
    if (length(index_ref_col) == 1) {
      index_geno <- G_chr[, index_ref_col]
      
      # Calculate r2 for each overlapping SNP
      for (i in seq_along(valid_snps)) {
        snp_idx <- valid_snps[i]
        ref_col <- which(map_chr$physical.pos == gwas_subset$pos[snp_idx])
        
        if (length(ref_col) == 1) {
          snp_geno <- G_chr[, ref_col]
          # Calculate correlation and square it for r2
          r_val <- cor(index_geno[], snp_geno[], use = "complete.obs")
          gwas_subset$r2[snp_idx] <- r_val^2
        }
      }
      
      # Set index SNP r2 to 1
      gwas_subset$r2[index_idx] <- 1.0
      
      cat("  LD calculated for", sum(!is.na(gwas_subset$r2)), "SNPs\n")
    } else {
      warning("Index SNP not found in reference panel")
    }
  } else {
    warning("No overlapping SNPs found between GWAS and reference panel")
  }
  
  # Remove temporary column
  gwas_subset$ref_idx <- NULL
  
  return(gwas_subset)
}

#' Create locus plot
#'
#' @param locus_data GWAS data for locus with LD information
#' @param locus_name Name of the locus
#' @param chr_num Chromosome number
#' @param start_pos Start position
#' @param end_pos End position
#' @param output_prefix Output file prefix
#' @param study_info Study information for title
create_locus_plot <- function(locus_data, locus_name, chr_num, start_pos, 
                               end_pos, output_prefix, study_info) {
  
  cat("\n  Creating locus plot...\n")
  
  # Create locus object
  loc <- locus(
    data = locus_data,
    seqname = chr_num,
    xrange = c(start_pos, end_pos),
    LD = "r2",
    ens_db = "EnsDb.Hsapiens.v75"
  )
  
  # Get index SNP rsid for label
  index_idx <- which.min(locus_data$p)
  index_rsid <- locus_data$rsid[index_idx]
  
  # Create plot title
  plot_title <- paste0(study_info$study_name, " - ", locus_name, "\n",
                      "Chr ", chr_num, ":", format(start_pos, big.mark = ","), 
                      "-", format(end_pos, big.mark = ","), "\n",
                      "Cases: ", study_info$n_cases, ", Controls: ", 
                      study_info$n_controls, " | LD: ", study_info$ref_panel_name)
  
  # Base graphics plot
  pdf(paste0(output_prefix, "_base.pdf"), width = 12, height = 7)
  locus_plot(loc, 
             border = TRUE,
             labels = index_rsid,
             main = plot_title)
  dev.off()
  
  # ggplot2 version
  pdf(paste0(output_prefix, "_ggplot.pdf"), width = 12, height = 7)
  
  # Create scatter plot with threshold lines
  scatter <- gg_scatter(loc, labels = index_rsid) +
    ggplot2::geom_hline(yintercept = -log10(5e-8), linetype = "dashed", 
                       color = "red", linewidth = 0.5) +
    ggplot2::geom_hline(yintercept = -log10(5e-5), linetype = "dashed", 
                       color = "green", linewidth = 0.5)
  
  # Add gene tracks below
  p <- gg_addgenes(scatter, loc) +
    ggplot2::ggtitle(plot_title) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 10, hjust = 0.5))
  
  print(p)
  dev.off()
  
  cat("  Plots saved:\n")
  cat("    ", paste0(output_prefix, "_base.pdf\n"))
  cat("    ", paste0(output_prefix, "_ggplot.pdf\n"))
  
  return(loc)
}

################################################################################
# MAIN SCRIPT
################################################################################

cat("################################################################################\n")
cat("LocusZoom Plot Generator\n")
cat("Study:", STUDY_NAME, "\n")
cat("Cases:", N_CASES, "| Controls:", N_CONTROLS, "\n")
cat("Reference Panel:", REF_PANEL_NAME, "\n")
cat("################################################################################\n\n")

# Create output directory
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
  cat("Created output directory:", OUTPUT_DIR, "\n")
}

# Load GWAS data
cat("\n[1/4] Loading GWAS summary statistics...\n")
gwas_data <- fread(GWAS_FILE)
cat("Loaded", nrow(gwas_data), "SNPs from", GWAS_FILE, "\n")

# Ensure column names are correct
if ("rsid" %in% colnames(gwas_data)) {
  setnames(gwas_data, old = c("rsid", "chr", "pos", "p"), 
           new = c("rsid", "chrom", "pos", "p"), skip_absent = TRUE)
} else {
  stop("GWAS file must contain columns: rsid, chr, pos, p")
}

# Convert p-values to numeric
gwas_data$p <- as.numeric(gwas_data$p)

# Load loci file
cat("\n[2/4] Loading loci to plot...\n")
loci <- fread(LOCI_FILE)
cat("Loaded", nrow(loci), "loci from", LOCI_FILE, "\n")
print(loci)

# Load reference panel
cat("\n[3/4] Loading reference panel...\n")
ref_panel <- load_reference_panel(REF_PANEL_PATH)

# Store study info
study_info <- list(
  study_name = STUDY_NAME,
  n_cases = N_CASES,
  n_controls = N_CONTROLS,
  ref_panel_name = REF_PANEL_NAME
)

# Process each locus
cat("\n[4/4] Processing loci and creating plots...\n")
cat("================================================================================\n")

for (i in 1:nrow(loci)) {
  locus_name <- loci$Locus[i]
  chr_num <- loci$Chr[i]
  original_start <- loci$Start[i]
  original_end <- loci$End[i]
  
  cat("\nProcessing locus", i, "of", nrow(loci), ":", locus_name, "\n")
  cat("  Original region: chr", chr_num, ":", original_start, "-", original_end, "\n")
  
  # Check if region needs expansion
  region_size <- original_end - original_start
  cat("  Region size:", format(region_size, big.mark = ","), "bp\n")
  
  start_pos <- original_start
  end_pos <- original_end
  
  if (AUTO_EXPAND_REGIONS && region_size < MIN_REGION_SIZE) {
    # Calculate how much to expand
    expansion_needed <- MIN_REGION_SIZE - region_size
    flank <- ceiling(expansion_needed / 2)
    
    start_pos <- max(1, original_start - flank)
    end_pos <- original_end + flank
    
    cat("  Region expanded to:", format(MIN_REGION_SIZE, big.mark = ","), 
        "bp (adding", format(flank, big.mark = ","), "bp flanks)\n")
    cat("  Expanded region: chr", chr_num, ":", start_pos, "-", end_pos, "\n")
    cat("  NOTE: Index SNP will be selected only from original region\n")
  }
  
  # Extract GWAS data for this locus (use expanded region)
  locus_gwas <- gwas_data[chrom == chr_num & pos >= start_pos & pos <= end_pos]
  
  if (nrow(locus_gwas) == 0) {
    cat("  WARNING: No SNPs found in this region. Skipping...\n")
    next
  }
  
  cat("  Found", nrow(locus_gwas), "SNPs in this region\n")
  
  # Calculate LD (pass original region boundaries if expanded)
  if (AUTO_EXPAND_REGIONS && region_size < MIN_REGION_SIZE) {
    locus_gwas <- calculate_ld(locus_gwas, ref_panel, chr_num, 
                               original_start, original_end)
  } else {
    locus_gwas <- calculate_ld(locus_gwas, ref_panel, chr_num)
  }
  
  # Create output file prefix
  output_prefix <- file.path(OUTPUT_DIR, 
                             paste0(locus_name, "_chr", chr_num, "_", 
                                   start_pos, "-", end_pos))
  
  # Create plot
  tryCatch({
    create_locus_plot(locus_gwas, locus_name, chr_num, start_pos, end_pos,
                     output_prefix, study_info)
  }, error = function(e) {
    cat("  ERROR creating plot:", e$message, "\n")
  })
  
  cat("  Completed:", locus_name, "\n")
  cat("--------------------------------------------------------------------------------\n")
}

cat("\n################################################################################\n")
cat("All loci processed!\n")
cat("Output directory:", OUTPUT_DIR, "\n")
cat("################################################################################\n")
