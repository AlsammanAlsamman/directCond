#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

log_msg <- function(...) {
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(sprintf("[%s] ", ts), ..., "\n", sep = "")
  flush(stdout())
}

get_arg <- function(flag, required = TRUE, default = NA_character_) {
  idx <- which(args == flag)
  if (length(idx) == 0) {
    if (required) stop(paste0("Missing required argument: ", flag))
    return(default)
  }
  if (idx[1] == length(args)) stop(paste0("Missing value for argument: ", flag))
  args[idx[1] + 1]
}

analysis_name <- get_arg("--analysis")
gwas_file <- get_arg("--gwas-file")
cojo_dir <- get_arg("--cojo-dir")
out_xlsx <- get_arg("--out-xlsx")

if (!requireNamespace("openxlsx", quietly = TRUE)) {
  stop("Package openxlsx is required to generate the Excel workbook.")
}
if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("Package data.table is required (fread-based GWAS extraction).")
}

trim_text <- function(x) {
  gsub("^[[:space:]]+|[[:space:]]+$", "", x)
}

read_nonempty_lines <- function(path) {
  if (!file.exists(path)) return(character())
  lines <- readLines(path, warn = FALSE)
  lines <- trim_text(lines)
  lines[nzchar(lines)]
}

parse_done_file <- function(path) {
  fields <- list()
  if (!file.exists(path)) return(fields)
  line <- readLines(path, n = 1, warn = FALSE)
  if (!length(line)) return(fields)
  tokens <- strsplit(line, "\t", fixed = TRUE)[[1]]
  if (length(tokens) == 0) return(fields)
  fields$status <- tokens[1]
  if (length(tokens) > 1) {
    for (token in tokens[-1]) {
      token <- trim_text(token)
      if (!nzchar(token) || !grepl("=", token, fixed = TRUE)) next
      parts <- strsplit(token, "=", fixed = TRUE)[[1]]
      key <- parts[1]
      value <- paste(parts[-1], collapse = "=")
      fields[[key]] <- value
    }
  }
  fields
}

find_snp_column <- function(headers) {
  lowered <- tolower(trim_text(headers))
  candidates <- c("snp", "rsid", "varid", "marker", "id")
  idx <- match(candidates, lowered, nomatch = 0)
  idx <- idx[idx > 0]
  if (!length(idx)) stop("Could not find SNP identifier column in original GWAS file.")
  idx[1]
}

read_selected_gwas_rows <- function(path, target_snps) {
  target_snps <- unique(target_snps[nzchar(target_snps)])
  if (!length(target_snps)) return(list(headers = character(), rows = list()))

  log_msg("Reading GWAS via fread and filtering conditioned rsIDs: ", length(target_snps), " targets")

  dt <- data.table::fread(path, sep = "\t", quote = "", data.table = TRUE, showProgress = FALSE)
  headers <- names(dt)
  snp_col_idx <- find_snp_column(headers)
  snp_col <- headers[snp_col_idx]

  dt[, .snp_key := trim_text(as.character(get(snp_col)))]
  filtered <- dt[.snp_key %in% target_snps]
  if (nrow(filtered)) {
    filtered <- filtered[!duplicated(.snp_key)]
  }

  rows <- list()
  if (nrow(filtered)) {
    for (i in seq_len(nrow(filtered))) {
      key <- filtered$.snp_key[i]
      row_vals <- as.list(filtered[i, ..headers])
      rows[[key]] <- row_vals
    }
  }

  log_msg("GWAS rsID filter complete: found=", length(rows), "/", length(target_snps))
  list(headers = headers, rows = rows)
}

read_cojo_table <- function(path) {
  if (!file.exists(path)) return(data.frame(stringsAsFactors = FALSE))
  read.table(path, header = TRUE, sep = "", stringsAsFactors = FALSE, check.names = FALSE, quote = "", comment.char = "")
}

build_summary_row <- function(analysis_name, locus_id, cond_snps, done_fields) {
  conditioned_text <- if (length(cond_snps)) paste(cond_snps, collapse = ", ") else ""
  data.frame(
    analysis = analysis_name,
    locus_id = locus_id,
    n_conditioned = length(cond_snps),
    conditioned_snps = conditioned_text,
    mode = if (!is.null(done_fields$mode)) done_fields$mode else "",
    last_cond_snp = if (!is.null(done_fields$last_cond_snp)) done_fields$last_cond_snp else "",
    initial_candidates = if (!is.null(done_fields$initial_candidates)) done_fields$initial_candidates else "",
    remaining_candidates = if (!is.null(done_fields$remaining_candidates)) done_fields$remaining_candidates else "",
    skipped_collinear = if (!is.null(done_fields$skipped_collinear)) done_fields$skipped_collinear else "",
    status = if (!is.null(done_fields$status)) done_fields$status else "",
    stringsAsFactors = FALSE
  )
}

if (!dir.exists(cojo_dir)) stop("COJO directory not found: ", cojo_dir)
if (!file.exists(gwas_file)) stop("Original GWAS file not found: ", gwas_file)

log_msg("Starting conditioned SNP Excel export for analysis: ", analysis_name)
log_msg("GWAS file: ", gwas_file)
log_msg("COJO directory: ", cojo_dir)

locus_dirs <- list.dirs(cojo_dir, recursive = FALSE, full.names = TRUE)
locus_dirs <- locus_dirs[dir.exists(locus_dirs)]
locus_dirs <- locus_dirs[order(basename(locus_dirs))]

summary_rows <- list()
condition_rows <- list()
all_conditioned_snps <- character()
cojo_tables <- list()

for (locus_path in locus_dirs) {
  locus_id <- basename(locus_path)
  cond_file <- file.path(locus_path, "cojo.cond.snp")
  done_file <- file.path(locus_path, "cojo.done")
  final_file <- file.path(locus_path, "cojo.cma.cojo")

  cond_snps <- read_nonempty_lines(cond_file)
  done_fields <- parse_done_file(done_file)
  summary_rows[[length(summary_rows) + 1]] <- build_summary_row(analysis_name, locus_id, cond_snps, done_fields)

  if (file.exists(final_file)) {
    cojo_tables[[locus_id]] <- read_cojo_table(final_file)
  } else {
    cojo_tables[[locus_id]] <- data.frame(stringsAsFactors = FALSE)
  }

  if (!length(cond_snps)) next

  all_conditioned_snps <- c(all_conditioned_snps, cond_snps)
  for (idx in seq_along(cond_snps)) {
    condition_rows[[length(condition_rows) + 1]] <- data.frame(
      analysis = analysis_name,
      locus_id = locus_id,
      condition_order = idx,
      conditioned_snp = cond_snps[idx],
      stringsAsFactors = FALSE
    )
  }
}

log_msg("Loci discovered: ", length(locus_dirs), ", loci with conditioned SNPs: ", length(unique(vapply(condition_rows, function(x) x$locus_id, character(1)))))

summary_df <- if (length(summary_rows)) do.call(rbind, summary_rows) else data.frame(stringsAsFactors = FALSE)
detail_base <- if (length(condition_rows)) do.call(rbind, condition_rows) else data.frame(
  analysis = character(),
  locus_id = character(),
  condition_order = integer(),
  conditioned_snp = character(),
  stringsAsFactors = FALSE
)

gwas_lookup <- read_selected_gwas_rows(gwas_file, all_conditioned_snps)
gwas_headers <- gwas_lookup$headers
gwas_rows <- gwas_lookup$rows

detail_records <- list()
if (nrow(detail_base)) {
  for (row_idx in seq_len(nrow(detail_base))) {
    locus_id <- detail_base$locus_id[row_idx]
    snp <- detail_base$conditioned_snp[row_idx]
    gwas_row <- gwas_rows[[snp]]
    if (is.null(gwas_row)) {
      gwas_row <- as.list(setNames(rep("", length(gwas_headers)), gwas_headers))
    }

    cojo_tbl <- cojo_tables[[locus_id]]
    cojo_match <- if (nrow(cojo_tbl) && "SNP" %in% names(cojo_tbl)) {
      cojo_tbl[cojo_tbl$SNP == snp, , drop = FALSE]
    } else {
      data.frame(stringsAsFactors = FALSE)
    }

    cojo_fields <- list(
      cojo_bp = if (nrow(cojo_match) && "bp" %in% names(cojo_match)) as.character(cojo_match$bp[1]) else "",
      cojo_p = if (nrow(cojo_match) && "p" %in% names(cojo_match)) as.character(cojo_match$p[1]) else "",
      cojo_pC = if (nrow(cojo_match) && "pC" %in% names(cojo_match)) as.character(cojo_match$pC[1]) else "",
      cojo_b = if (nrow(cojo_match) && "b" %in% names(cojo_match)) as.character(cojo_match$b[1]) else "",
      cojo_bC = if (nrow(cojo_match) && "bC" %in% names(cojo_match)) as.character(cojo_match$bC[1]) else ""
    )

    detail_records[[length(detail_records) + 1]] <- c(
      list(
        analysis = detail_base$analysis[row_idx],
        locus_id = locus_id,
        condition_order = detail_base$condition_order[row_idx],
        conditioned_snp = snp
      ),
      gwas_row,
      cojo_fields
    )
  }
}

detail_df <- if (length(detail_records)) {
  as.data.frame(do.call(rbind, lapply(detail_records, function(record) {
    as.data.frame(record, stringsAsFactors = FALSE)
  })), stringsAsFactors = FALSE)
} else {
  data.frame(
    analysis = character(),
    locus_id = character(),
    condition_order = integer(),
    conditioned_snp = character(),
    stringsAsFactors = FALSE
  )
}

numeric_columns <- intersect(c("condition_order", "P", "N", "cojo_bp", "cojo_p", "cojo_pC", "cojo_b", "cojo_bC", "BETA", "SE", "OR", "POS", "CHR"), names(detail_df))
for (col_name in numeric_columns) {
  suppressWarnings(detail_df[[col_name]] <- as.numeric(detail_df[[col_name]]))
}

dir.create(dirname(out_xlsx), recursive = TRUE, showWarnings = FALSE)

wb <- openxlsx::createWorkbook(creator = "GitHub Copilot")

header_style <- openxlsx::createStyle(
  fgFill = "#1F4E78",
  fontColour = "#FFFFFF",
  textDecoration = "bold",
  halign = "center",
  border = "Bottom"
)
summary_style <- openxlsx::createStyle(fgFill = "#D9EAF7")
sig_style <- openxlsx::createStyle(fontColour = "#9C0006", fgFill = "#FFC7CE")
cond_style <- openxlsx::createStyle(fgFill = "#FFF2CC")

openxlsx::addWorksheet(wb, "Summary", tabColour = "#4F81BD")
openxlsx::writeDataTable(wb, "Summary", summary_df, tableStyle = "TableStyleMedium9", withFilter = TRUE)
openxlsx::addStyle(wb, "Summary", header_style, rows = 1, cols = seq_len(max(1, ncol(summary_df))), gridExpand = TRUE)
if (nrow(summary_df)) {
  openxlsx::conditionalFormatting(wb, "Summary", cols = which(names(summary_df) == "n_conditioned"), rows = 2:(nrow(summary_df) + 1), rule = ">0", style = summary_style)
}
openxlsx::freezePane(wb, "Summary", firstRow = TRUE)
openxlsx::setColWidths(wb, "Summary", cols = seq_len(max(1, ncol(summary_df))), widths = "auto")

openxlsx::addWorksheet(wb, "Conditioned_SNPs", tabColour = "#C0504D")
openxlsx::writeDataTable(wb, "Conditioned_SNPs", detail_df, tableStyle = "TableStyleMedium2", withFilter = TRUE)
openxlsx::addStyle(wb, "Conditioned_SNPs", header_style, rows = 1, cols = seq_len(max(1, ncol(detail_df))), gridExpand = TRUE)
openxlsx::freezePane(wb, "Conditioned_SNPs", firstRow = TRUE)
openxlsx::setColWidths(wb, "Conditioned_SNPs", cols = seq_len(max(1, ncol(detail_df))), widths = "auto")

if (nrow(detail_df)) {
  openxlsx::addStyle(wb, "Conditioned_SNPs", cond_style, rows = 2:(nrow(detail_df) + 1), cols = match("conditioned_snp", names(detail_df)), gridExpand = TRUE, stack = TRUE)

  if ("P" %in% names(detail_df)) {
    openxlsx::conditionalFormatting(wb, "Conditioned_SNPs", cols = match("P", names(detail_df)), rows = 2:(nrow(detail_df) + 1), rule = "<5E-8", style = sig_style)
  }
  if ("cojo_pC" %in% names(detail_df)) {
    openxlsx::conditionalFormatting(wb, "Conditioned_SNPs", cols = match("cojo_pC", names(detail_df)), rows = 2:(nrow(detail_df) + 1), rule = "<5E-8", style = sig_style)
  }
}

openxlsx::saveWorkbook(wb, out_xlsx, overwrite = TRUE)
log_msg("Wrote conditioned SNP Excel report: ", out_xlsx)