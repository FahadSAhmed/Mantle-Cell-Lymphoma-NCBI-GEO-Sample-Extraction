# ================================
# GEO MZL → GPL-aware normalize → Single Excel (streamed) → ROR1 plot
# Base dir: C:/Users/NAME/Downloads/MZL
# Sheets: one per MZL GSE block + ROR1_MZL summary
# Author: Fahad S. Ahmed, MD
# ================================

# ---- Config paths ----
base_dir  <- "C:/Users/NAME/Downloads/MZL"
cache_dir <- file.path(base_dir, "cache")   # GEO download cache
rds_dir   <- file.path(base_dir, "rds")     # per-GSE normalized matrices on disk
dir.create(base_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(rds_dir,   showWarnings = FALSE, recursive = TRUE)

infile_mzl <- file.path(base_dir, "gds_result_MZL.txt")

excel_path <- file.path(base_dir, "MZL_merged.xlsx")
plot_mzl   <- file.path(base_dir, "ROR1_expression_MZL.png")
log_path   <- file.path(base_dir, "run_log_MZL.txt")

# ---- Optional: reduce per-sheet size by keeping only selected genes (NULL = keep all)
GENE_WHITELIST <- NULL  # e.g., c("ROR1","MS4A1","CD19","CD79B","BTK","PAX5")

# ---- Packages ----
cran_pkgs <- c("stringr","readr","dplyr","tibble","purrr","janitor","ggplot2","tidyr","openxlsx")
to_install_cran <- setdiff(cran_pkgs, rownames(installed.packages()))
if (length(to_install_cran)) install.packages(to_install_cran, repos = "https://cloud.r-project.org")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
bioc_pkgs <- c("GEOquery","Biobase","limma","edgeR","annotate")
to_install_bioc <- setdiff(bioc_pkgs, rownames(installed.packages()))
if (length(to_install_bioc)) BiocManager::install(to_install_bioc, ask = FALSE, update = TRUE)

suppressPackageStartupMessages({
  library(GEOquery); library(Biobase); library(limma); library(edgeR); library(annotate)
  library(stringr);  library(readr);   library(dplyr); library(tibble); library(purrr)
  library(janitor);  library(openxlsx); library(ggplot2); library(tidyr)
})

# ---- Logging ----
log_con <- file(log_path, open = "wt")
on.exit(try(close(log_con), silent = TRUE), add = TRUE)
log_msg <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), sprintf(...)),
      file = log_con)
  flush(log_con)
}

# ---- Utilities ----
extract_gse <- function(txt) unique(str_extract_all(txt, "GSE\\d+") |> unlist())

is_cell_line <- function(x) {
  any(str_detect(tolower(paste(x, collapse = " | ")), "cell\\s*line|cell-line|immortalized"))
}
has_disease  <- function(x, keys) {
  s <- tolower(paste(x, collapse = " | "))
  any(str_detect(s, paste0("\\b(", paste(keys, collapse="|"), ")\\b")))
}

pick_patient_id <- function(pdat_row) {
  cand_names <- names(pdat_row); vals <- as.list(pdat_row)
  keys <- c("patient","patient_id","subject id","subjectid","sample id",
            "individual","case id","donor id","participant id")
  for (k in keys) {
    hit <- which(tolower(cand_names) == k |
                   str_detect(tolower(cand_names), paste0("\\b", k, "\\b")))
    if (length(hit)) {
      v <- as.character(vals[[hit[1]]])
      if (!is.na(v) && nzchar(v)) return(v)
    }
  }
  NA_character_
}

collapse_to_genes <- function(expr, fdat,
                              gene_cols = c("Gene.symbol","Gene Symbol","Symbol","SYMBOL",
                                            "GENE_SYMBOL","GeneID","Gene name","GENE","Entrez_Gene_ID")) {
  genes <- NULL
  if (!is.null(fdat) && nrow(fdat) == nrow(expr)) {
    for (gc in gene_cols) {
      if (gc %in% colnames(fdat)) { genes <- fdat[[gc]]; break }
    }
  }
  if (is.null(genes)) {
    rn <- rownames(expr)
    genes <- if (!is.null(rn) && mean(str_detect(rn, "^[A-Za-z0-9._-]+$")) > 0.9)
      rn else paste0("FEATURE_", seq_len(nrow(expr)))
  }
  genes <- make.names(genes)
  keep <- !is.na(genes) & genes != "" & genes != "NA"
  expr <- expr[keep, , drop=FALSE]; genes <- genes[keep]
  ag <- rowsum(expr, group = genes, reorder = FALSE) / as.vector(table(genes))
  as.matrix(ag)
}

pheno_strings <- function(pheno_row) {
  nms <- names(pheno_row)
  sel <- str_detect(tolower(nms),
                    "title|source|characteristics|description|disease|diagnos|sample")
  as.character(unlist(pheno_row[sel]))
}
is_MZL      <- function(pheno_row) has_disease(pheno_strings(pheno_row),
                                               c("mzl","marginal zone","splenic marginal",
                                                 "nodal marginal","malt"))
is_CellLine <- function(pheno_row) is_cell_line(pheno_strings(pheno_row))

# ---- Hardened numeric helpers ----
finite_numeric_matrix <- function(m) {
  m <- suppressWarnings(as.matrix(m))
  storage.mode(m) <- "numeric"
  m[!is.finite(m)] <- NA_real_
  m
}
safe_q75 <- function(m) {
  v <- as.numeric(m[is.finite(m)])
  if (!length(v)) return(NA_real_)
  suppressWarnings(stats::quantile(v, 0.75, na.rm = TRUE, names = FALSE))
}
to_log2_if_needed_safe <- function(m) {
  q75 <- safe_q75(m)
  if (is.na(q75)) return(m)
  if (q75 > 100) log2(m + 1) else m
}
is_counts_like <- function(mat) {
  mat <- finite_numeric_matrix(mat)
  if (sum(is.finite(mat)) < 2L) return(FALSE)
  is_integerish <- function(v){
    v <- v[is.finite(v)]
    if (!length(v)) return(FALSE)
    mean(abs(v - round(v)) < 1e-6) > 0.9
  }
  pct_int_cols <- mean(apply(mat, 2, is_integerish))
  big_range <- {
    r <- apply(mat, 2, function(x){
      x <- x[is.finite(x)]
      if (!length(x)) return(0)
      diff(range(x))
    })
    median(r, na.rm=TRUE) > 1000
  }
  many_zeros <- mean(mat == 0, na.rm=TRUE) > 0.05
  isTRUE(pct_int_cols > 0.6 & (big_range | many_zeros))
}

# ---- GPL-aware platform detection ----
detect_gpl <- function(eset) {
  gpl <- tryCatch(Meta(eset)$platform_id, error = function(e) NA_character_)
  if (is.na(gpl)) {
    pd <- tryCatch(pData(eset), error = function(e) NULL)
    if (!is.null(pd)) {
      cand <- intersect(colnames(pd), c("platform_id","Platform_geo_accession"))
      if (length(cand)) gpl <- as.character(pd[[cand[1]]][1])
    }
  }
  if (is.na(gpl)) gpl <- tryCatch(eset@experimentData@other$platform_id, error = function(e) NA_character_)
  gpl
}
platform_class <- function(gpl) {
  if (is.na(gpl) || !nzchar(gpl))
    return(list(class="Unknown", channels=NA_character_, vendor=NA_character_))
  gpl_meta <- tryCatch({ g <- getGPL(gpl, destdir = cache_dir); Meta(g) }, error = function(e) NULL)
  title <- if (!is.null(gpl_meta) && !is.null(gpl_meta$title)) gpl_meta$title else ""
  manuf <- if (!is.null(gpl_meta) && !is.null(gpl_meta$manufacturer)) gpl_meta$manufacturer else ""
  channels <- if (!is.null(gpl_meta) && !is.null(gpl_meta$channel_count))
    as.character(gpl_meta$channel_count) else NA_character_
  txt <- paste(title, manuf, gpl)
  vendor <- dplyr::case_when(
    grepl("Affymetrix|GeneChip|HTA|Clariom", txt, ignore.case = TRUE) ~ "Affymetrix",
    grepl("Agilent", txt, ignore.case = TRUE)                         ~ "Agilent",
    grepl("Illumina|BeadChip|HT-12|WG-6", txt, ignore.case = TRUE)    ~ "Illumina",
    grepl("NanoString", txt, ignore.case = TRUE)                      ~ "NanoString",
    TRUE                                                              ~ "Unknown"
  )
  class <- vendor
  if (vendor == "Agilent" && !is.na(channels)) {
    if (channels %in% c("2","two","Two-color","two-color")) class <- "Agilent_2color"
    if (channels %in% c("1","one","One-color","one-color")) class <- "Agilent_1color"
  }
  list(class = class, channels = channels, vendor = vendor)
}

# ---- Platform-aware normalization ----
normalize_platform_matrix <- function(eset, expr) {
  expr <- finite_numeric_matrix(expr)
  if (sum(is.finite(expr)) < 2L || ncol(expr) == 0L || nrow(expr) == 0L) {
    attr(expr, "norm_method") <- "pass_through_insufficient_data"
    if (exists("log_msg")) log_msg("Insufficient finite values; passing through.")
    return(expr)
  }
  gpl <- detect_gpl(eset); pf <- platform_class(gpl)
  if (exists("log_msg"))
    log_msg("GPL %s detected: class=%s, vendor=%s, channels=%s",
            gpl, pf$class, pf$vendor, pf$channels)
  
  if (is_counts_like(expr)) {
    dge <- edgeR::DGEList(counts = expr)
    dge <- edgeR::calcNormFactors(dge, method = "TMM")
    out <- edgeR::cpm(dge, log = TRUE, prior.count = 1)
    attr(out, "norm_method") <- "edgeR_TMM_logCPM"
    return(out)
  }
  
  cls <- pf$class; tolog <- to_log2_if_needed_safe
  if (cls == "Affymetrix") {
    m <- tolog(expr); out <- limma::normalizeBetweenArrays(m, "quantile")
    attr(out,"norm_method") <- "Affymetrix_series_quantile"; return(out)
  }
  if (cls == "Illumina") {
    m <- tolog(expr); out <- limma::normalizeBetweenArrays(m, "quantile")
    attr(out,"norm_method") <- "Illumina_series_quantile";   return(out)
  }
  if (cls == "Agilent_1color") {
    m <- tolog(expr); out <- limma::normalizeBetweenArrays(m, "quantile")
    attr(out,"norm_method") <- "Agilent_1color_series_quantile"; return(out)
  }
  if (cls == "Agilent_2color") {
    has_negs <- mean(expr < 0, na.rm = TRUE) > 0.1
    if (has_negs) {
      out <- limma::normalizeBetweenArrays(expr, "Aquantile")
      attr(out,"norm_method") <- "Agilent_2color_Aquantile"; return(out)
    }
    m <- tolog(expr); out <- limma::normalizeBetweenArrays(m, "quantile")
    attr(out,"norm_method") <- "Agilent_2color_series_quantile_fallback"; return(out)
  }
  if (cls == "NanoString") {
    m <- tolog(expr)
    if (sum(is.finite(m)) < 2L) {
      attr(m,"norm_method") <- "NanoString_pass_through"; return(m)
    }
    m <- sweep(m, 2, apply(m, 2, function(x) median(x, na.rm=TRUE)), "-")
    attr(m,"norm_method") <- "NanoString_series_log2_medianCenter"; return(m)
  }
  m <- tolog(expr)
  out <- limma::normalizeBetweenArrays(m, "quantile")
  attr(out,"norm_method") <- "Default_series_quantile"
  return(out)
}

# ---- Cache check ----
series_matrix_present <- function(gse_id) {
  any(grepl(paste0("^", gse_id, ".*series_matrix.*gz$"),
            list.files(cache_dir, full.names = FALSE)))
}

# ---- Process one GSE → write RDS, return minimal metadata + ROR1 ----
process_one_GSE <- function(gse_id, disease_flag_fun) {
  log_msg("Fetching %s ...", gse_id)
  if (series_matrix_present(gse_id))
    log_msg("Series Matrix cached for %s → skip re-download.", gse_id)
  
  esets <- tryCatch(getGEO(gse_id, GSEMatrix = TRUE, AnnotGPL = TRUE,
                           destdir = cache_dir), error = function(e) NULL)
  if (is.null(esets)) { log_msg("No Series Matrix for %s, skipping.", gse_id); return(NULL) }
  if (inherits(esets, "ExpressionSet")) esets <- list(esets)
  
  out <- list()
  for (i in seq_along(esets)) {
    eset <- esets[[i]]
    pdat <- pData(eset) |> tibble::as_tibble(rownames = "GSM")
    keep <- purrr::map_lgl(seq_len(nrow(pdat)), ~ disease_flag_fun(pdat[.x,])) &
      !purrr::map_lgl(seq_len(nrow(pdat)), ~ is_CellLine(pdat[.x,]))
    if (!any(keep)) next
    
    psub <- pdat[keep, , drop = FALSE]
    expr <- exprs(eset)[, psub$GSM, drop = FALSE]
    expr <- finite_numeric_matrix(expr)
    expr <- expr[, colSums(is.finite(expr)) > 0, drop = FALSE]
    if (ncol(expr) == 0L || nrow(expr) == 0L) {
      log_msg("Empty matrix after cleanup; skipping %s[%d].", gse_id, i)
      next
    }
    
    expr_n <- normalize_platform_matrix(eset, expr)
    fdat   <- tryCatch(fData(eset), error = function(e) NULL)
    expr_g <- collapse_to_genes(expr_n, fdat)
    
    if (!is.null(GENE_WHITELIST)) {
      rn <- toupper(rownames(expr_g))
      keep_rows <- rn %in% toupper(GENE_WHITELIST)
      expr_g <- expr_g[keep_rows, , drop = FALSE]
    }
    
    expr_g <- expr_g[, colnames(expr_g) %in% psub$GSM, drop=FALSE]
    psub   <- psub |> dplyr::filter(GSM %in% colnames(expr_g))
    if (!nrow(psub) || !ncol(expr_g)) {
      log_msg("No overlap after collapse; skipping %s[%d].", gse_id, i)
      next
    }
    
    psub <- psub |>
      dplyr::mutate(
        patient_id = purrr::map_chr(seq_len(dplyr::n()), ~ pick_patient_id(psub[.x,])),
        patient_id = ifelse(is.na(patient_id) | !nzchar(patient_id), GSM, patient_id)
      ) |>
      dplyr::distinct(patient_id, .keep_all = TRUE)
    
    expr_g <- expr_g[, psub$GSM, drop=FALSE]
    if (!nrow(psub) || !ncol(expr_g)) next
    
    block_id <- sprintf("%s_%d", gse_id, i)
    saveRDS(list(expr = expr_g, meta = psub),
            file = file.path(rds_dir, paste0(block_id, ".rds")))
    
    ror1_value <- NA_real_
    rn_up <- toupper(rownames(expr_g))
    idx <- which(rn_up == "ROR1")
    if (length(idx)) ror1_value <- as.numeric(expr_g[idx[1], ])
    out[[block_id]] <- list(
      block_id = block_id,
      gse_id   = gse_id,
      gsms     = colnames(expr_g),
      ror1     = ror1_value
    )
    
    rm(expr, expr_n, fdat, expr_g, psub); invisible(gc())
  }
  rm(esets); invisible(gc())
  if (!length(out)) return(NULL)
  out
}

# ---- Read MZL list ----
if (!file.exists(infile_mzl)) stop("MZL list not found: ", infile_mzl)
txt_mzl <- read_file(infile_mzl)
gse_mzl <- extract_gse(txt_mzl)
log_msg("Candidate series — MZL: %d", length(gse_mzl))

# ---- Harvest MZL (streamed, with periodic gc) ----
harvest_disease <- function(gse_ids, flag_fun, clean_every = 3L) {
  out <- list(); k <- 0L
  for (i in seq_along(gse_ids)) {
    gid <- gse_ids[[i]]
    res <- process_one_GSE(gid, flag_fun)
    if (!is.null(res)) out <- c(out, res)
    k <- k + 1L
    if (k %% clean_every == 0L) {
      log_msg("Cleanup checkpoint after %d series", k)
      invisible(gc())
    }
  }
  if (!length(out)) return(NULL)
  out
}
mzl_blocks <- harvest_disease(gse_mzl, is_MZL, clean_every = 3L)

# ---- Build Excel workbook (MZL only, sheet-by-sheet) ----
wb <- createWorkbook()

write_blocks_to_excel <- function(blocks, disease_prefix) {
  if (is.null(blocks) || !length(blocks)) return(invisible(NULL))
  for (b in blocks) {
    rds_path <- file.path(rds_dir, paste0(b$block_id, ".rds"))
    if (!file.exists(rds_path)) next
    obj <- readRDS(rds_path)
    df <- obj$expr |> as.data.frame() |> tibble::rownames_to_column("Gene")
    sheet_name <- substr(paste0(disease_prefix, "_", b$block_id), 1, 31)  # Excel sheet limit
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet = sheet_name, x = df, rowNames = FALSE)
    rm(obj, df); invisible(gc())
  }
}

write_blocks_to_excel(mzl_blocks, "MZL")

# ---- Compact ROR1 sheet (MZL only) ----
ror_table <- function(blocks) {
  empty <- tibble::tibble(
    Block = character(), GSE = character(), GSM = character(), ROR1 = double()
  )
  if (is.null(blocks) || !length(blocks)) return(empty)
  
  out <- purrr::map_dfr(blocks, function(b) {
    v <- b$ror1
    if (is.null(v) || length(v) == 0L || all(is.na(v))) return(empty)
    n_gsm <- length(b$gsms)
    n_v   <- length(v)
    n     <- max(0L, min(n_gsm, n_v))
    if (n == 0L) return(empty)
    tibble::tibble(
      Block = rep(b$block_id, n),
      GSE   = rep(b$gse_id,   n),
      GSM   = b$gsms[seq_len(n)],
      ROR1  = as.numeric(v[seq_len(n)])
    )
  })
  
  if (!nrow(out)) return(empty)
  dplyr::filter(out, is.finite(ROR1))
}

mzl_ror1 <- ror_table(mzl_blocks)

addWorksheet(wb, "ROR1_MZL")
if (nrow(mzl_ror1))
  writeData(wb, "ROR1_MZL", mzl_ror1, rowNames = FALSE)

# ---- Save Excel (single workbook) ----
saveWorkbook(wb, excel_path, overwrite = TRUE)
log_msg("Excel written: %s", excel_path)

# ---- ROR1 plot (MZL only) ----
plot_ror1 <- function(df, disease_name, outfile) {
  if (is.null(df) || !nrow(df)) {
    log_msg("No ROR1 values for %s", disease_name); return(invisible(NULL))
  }
  p <- ggplot(df, aes(x = "", y = ROR1)) +
    geom_violin(trim = TRUE) +
    geom_boxplot(width = 0.2, outlier.shape = NA) +
    geom_jitter(width = 0.1, alpha = 0.5) +
    labs(title = paste0("ROR1 expression: ", disease_name),
         x = NULL, y = "Normalized expression (log-scale where applicable)") +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_blank(),
          plot.title = element_text(hjust = 0.5))
  ggsave(outfile, p, width = 5, height = 4, dpi = 300)
  log_msg("Plot written: %s", outfile)
}

plot_ror1(mzl_ror1, "MZL", plot_mzl)

# ---- Final cleanup ----
rm(list = setdiff(ls(), c("excel_path","plot_mzl","log_path","base_dir"))); invisible(gc())
log_msg("DONE.")
message("Outputs:\n",
        " - Excel workbook (MZL): ", excel_path, "\n",
        " - Plot: ", plot_mzl, "\n",
        " - Log: ", log_path, "\n",
        " - Per-GSE RDS (normalized): ", file.path(base_dir, "rds"))

