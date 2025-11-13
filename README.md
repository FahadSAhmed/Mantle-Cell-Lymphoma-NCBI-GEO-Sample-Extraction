# Mantle-Cell-Lymphoma-NCBI-GEO-Sample-Extraction
R pipeline to harvest, normalize, and analyze Marginal Zone Lymphoma (MZL) datasets from NCBI GEO. Platform-aware processing, gene-level aggregation, and streaming I/O generate a unified Excel workbook and ROR1 expression summaries, including per-GSE sheets, a compact ROR1 table, and diagnostic plots. 

# MZL ROR1 Meta-Analysis Pipeline

This repository contains an R pipeline to harvest Marginal Zone Lymphoma (MZL) datasets from NCBI GEO, normalize expression data in a platform-aware manner, and summarize ROR1 expression across studies into a single Excel workbook and diagnostic plots.

## Overview

Given a text file listing GEO series (e.g. `gds_result_MZL.txt` exported from GEO DataSets), the script:

1. **Parses candidate GSE IDs** from the text input.
2. **Downloads and caches** GEO Series Matrix files (GSEMatrix) to avoid redundant network calls.
3. **Filters samples to MZL only** based on phenotype annotations and excludes obvious cell lines.
4. **Per-platform normalization**:
   - Infers GPL/platform and vendor (Affymetrix, Illumina, Agilent 1-color / 2-color, NanoString, counts).
   - Applies appropriate normalization:
     - Microarray: log2 transform (if needed) + `limma::normalizeBetweenArrays` (quantile/Aquantile).
     - Count-like data: `edgeR::calcNormFactors` (TMM) + log-CPM.
     - NanoString: log2 + median centering.
5. **Collapses probes to gene-level** using available GPL annotations (e.g. `Gene.symbol`), aggregating probes by gene.
6. **Streams per-GSE blocks to disk** as RDS objects instead of constructing one monolithic matrix in memory (to avoid multi-GB allocations).
7. **Compiles a single Excel workbook**:
   - One sheet per MZL GSE block (gene × sample expression matrices).
   - A compact summary sheet (`ROR1_MZL`) with per-sample ROR1 values and metadata (GSE, GSM, block).
8. **Generates a violin/box/jitter plot** of ROR1 expression across all MZL samples.

The pipeline is designed to be conservative with memory and robust to heterogeneous platforms and incomplete annotations.

## Key Features

- **Disease-specific filtering**: Retains only MZL cases based on text mining of GEO phenotype annotations.
- **Cell line exclusion**: Removes samples likely annotated as cell lines or immortalized lines.
- **Platform-aware normalization** via GPL introspection:
  - Affymetrix, Illumina, Agilent 1-color and 2-color arrays
  - NanoString
  - Count-like RNA-seq / expression tables
- **Robust numeric handling**:
  - Safeguards around NA/Inf values
  - Adaptive log2 conversion (based on 75th percentile)
  - Detection of count-like matrices
- **Streaming I/O**:
  - Per-GSE gene-level matrices written to RDS (`rds/`)
  - Excel workbook constructed sheet-by-sheet to avoid large in-memory unions
- **ROR1-centric outputs**:
  - Per-sample ROR1 values for downstream exploratory analysis
  - Simple summary plot of ROR1 distribution in MZL

## Requirements

- **R** ≥ 4.2
- CRAN packages:
  - `stringr`, `readr`, `dplyr`, `tibble`, `purrr`, `janitor`, `ggplot2`, `tidyr`, `openxlsx`
- Bioconductor packages:
  - `GEOquery`, `Biobase`, `limma`, `edgeR`, `annotate`

The script automatically installs missing packages from CRAN and Bioconductor.

## Input

- `gds_result_MZL.txt`  
  Text file containing GEO DataSets search results for MZL, including GSE identifiers (e.g. “GSE12345”).

Place this file in the configured base directory, e.g.:

```r
base_dir  <- "C:/Users/USERNAME/Downloads/ROR1_CLL_MZL"
infile_mzl <- file.path(base_dir, "gds_result_MZL.txt")

Outputs

All outputs are written under base_dir:

Excel workbook

MZL_merged.xlsx

One sheet per MZL GSE block: MZL_<GSE>_<index>

Summary sheet: ROR1_MZL with columns:

Block (GSE block ID)

GSE (series accession)

GSM (sample accession)

ROR1 (normalized expression)

Plots

ROR1_expression_MZL.png
Violin + box + jitter plot of ROR1 across all retained MZL samples.

Log

run_log_MZL.txt
Text log capturing GSE processing, platform detection, filtering decisions, and any warnings.

Per-GSE normalized matrices (RDS)

rds/<GSE>_<index>.rds
Each RDS contains a list with:

expr: gene-level normalized expression matrix (genes × samples)

meta: corresponding sample-level metadata (phenotype)

Usage

Export or construct gds_result_MZL.txt from NCBI GEO (searching for MZL / marginal zone lymphoma).

Adjust base_dir in the script to your desired working directory.

Run the R script:

source("run_MZL_pipeline.R")


Inspect outputs:

Open MZL_merged.xlsx in Excel or R for cross-study expression review.

Review ROR1_expression_MZL.png for a global view of ROR1 expression.

Check run_log_MZL.txt for any datasets skipped due to missing data or annotation.

Customization

Gene whitelist

To restrict output matrices (e.g. to a focused biomarker panel), set:

GENE_WHITELIST <- c("ROR1","MS4A1","CD19","CD79B","BTK","PAX5")


Only these genes (if present) will be retained per GSE block and written to Excel.

Phenotype filters

The MZL classifier is implemented in is_MZL() and uses free-text pattern matching. This function can be refined if specific subtypes or exclusion criteria are needed.

Normalization behavior

Platform→normalization mapping is encapsulated in normalize_platform_matrix(). This can be modified if a different normalization strategy is preferred for a given platform class.

Caveats

The disease classification relies on GEO phenotype text; misannotated series may be included or excluded.

ROR1 detection depends on correct gene symbol mapping through GPL annotation; incomplete or non-standard annotations may lead to missing ROR1 values for some platforms.

The script assumes that series matrix files contain expression data (not raw CEL/IDAT files).
