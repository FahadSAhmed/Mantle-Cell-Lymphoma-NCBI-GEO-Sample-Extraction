
# **MZL ROR1 Analysis Pipeline**

This repository contains a complete R pipeline for **systematic retrieval, normalization, and aggregation of Marginal Zone Lymphoma (MZL) gene-expression datasets** from NCBI GEO, with a specific focus on **ROR1** expression across studies.  
The workflow is designed to handle heterogeneous data platforms, avoid high-memory operations, and generate unified outputs suitable for downstream analysis.

---

## **📌 Pipeline Summary**

Given a text file listing GEO series identifiers (e.g., *gds_result_MZL.txt* produced from GEO DataSets; step 0 is a manual process), the script performs:

0. **Producing the gsd_list_MZL.txt**
   - Go to "https://www.ncbi.nlm.nih.gov/gds" and search for "Marginal Zone lymphoma" in the search bar
   - Saving the list to file "gds_list_MZL.txt"
      
2. **GSE Extraction**
   - Reads GEO search results and extracts all valid GSE accessions.

3. **Automated GEO Downloading + Caching**
   - Downloads *Series Matrix* files only once.
   - Reuses previously downloaded versions to reduce run-time and bandwidth.

4. **MZL Sample Filtering**
   - Identifies MZL disease samples using robust text-based phenotype parsing.
   - Removes cell-line or immortalized samples.

5. **Platform-Aware Normalization**
   - Automatically detects GPL platform → maps to appropriate normalization:
     - Affymetrix → log2 + quantile normalization  
     - Illumina → log2 + quantile normalization  
     - Agilent 1-color / 2-color → quantile or Aquantile  
     - NanoString → log2 + median centering  
     - Count-like data → edgeR TMM + logCPM  
     - Unknown → safe fallback normalization

6. **Gene-Level Aggregation**
   - Collapses probes → genes using GPL annotations (`Gene.symbol`, `SYMBOL`, etc.).
   - Averages probes mapping to the same gene.

7. **Streaming-Based Data Handling**
   - Each normalized GSE block is written as a separate **RDS file**.
   - Excel workbook is assembled sheet-by-sheet → avoids multi-GB memory use.

8. **ROR1 Summary Extraction**
   - Extracts ROR1 expression for each sample.
   - Builds a compact summary table in Excel.

9. **Plot Generation**
   - Produces a violin + box + jitter plot of ROR1 expression across all MZL samples.

---

## **📁 Directory Structure**

```
base_dir/
│
├── gds_result_MZL.txt               # Input list of GSE datasets
├── MZL_merged.xlsx                  # Output Excel workbook
├── ROR1_expression_MZL.png          # ROR1 distribution plot
├── run_log_MZL.txt                  # Full processing log
│
└── rds/                              # Normalized per-GSE blocks
      ├── GSEXXXXX_1.rds
      ├── GSEXXXXX_2.rds
      └── ...
```

---

## **📘 Input Requirements**

### **1. MZL Input File**
A text file containing GEO DataSets export results:

```
gds_result_MZL.txt
```

Must contain strings like:  
```
GSE12345
GSE67890
...
```

The script automatically detects and extracts them.

### **2. Directory Setup**
Set the working directory in the script:

```r
base_dir <- "C:/Users/USERNAME/Downloads/ROR1_CLL_MZL"
```

---

## **📦 Dependencies**

The script automatically installs missing packages.

### **CRAN**
- `stringr`
- `readr`
- `dplyr`
- `tibble`
- `purrr`
- `janitor`
- `tidyr`
- `openxlsx`
- `ggplot2`

### **Bioconductor**
- `GEOquery`
- `Biobase`
- `limma`
- `edgeR`
- `annotate`

---

## **🧠 Core Features**

### **✔ Robust Disease Filtering**
- Recognizes MZL, marginal zone lymphoma, SMZL, NMZL, MALT, etc.
- Rejects cell-line samples.

### **✔ Platform-Aware Normalization**
The pipeline recognizes platform vendors through the GPL metadata:

| Vendor        | Class                | Method |
|---------------|----------------------|--------|
| Affymetrix    | Microarray           | log2 + quantile |
| Illumina      | BeadChip             | log2 + quantile |
| Agilent 1-color | Microarray         | quantile |
| Agilent 2-color | Microarray         | Aquantile or fallback |
| NanoString    | nCounter             | log2 + median centering |
| RNA-seq like  | Count-like           | edgeR TMM + logCPM |

### **✔ Memory-Safe Handling**
Avoids massive `rbind`/`cbind` operations.  
All blocks are saved to RDS and appended individually into Excel.

### **✔ Gene Collapsing**
Handles datasets that provide probes instead of gene symbols.

### **✔ Clean and Reproducible Outputs**
- One Excel sheet per normalized GSE block  
- One summary sheet: `ROR1_MZL`  
- One high-resolution plot: `ROR1_expression_MZL.png`

---

## **📊 Output Files**

### **1. Excel Workbook: `MZL_merged.xlsx`**

Contains:

#### **Per-GSE Sheets**
Example names:
```
MZL_GSE12345_1
MZL_GSE67890_1
MZL_GSE67890_2
...
```

Each sheet has:

| Gene | GSM1 | GSM2 | ... |
|------|------|------|-----|
| ROR1 | ...  | ...  | ... |
| CD79A | ... | ...  | ... |

#### **Summary Sheet**
```
ROR1_MZL
```

| Block | GSE | GSM | ROR1 |
|-------|-----|-----|------|
| GSE12345_1 | GSE12345 | GSM10001 | 8.23 |
| ... | ... | ... | ... |

---

## **📈 Plot: ROR1 Expression**

Generated file:

```
ROR1_expression_MZL.png
```

This plot displays:
- Distribution of ROR1 across all MZL samples  
- Violin plot (density)
- Boxplot (median, IQR)
- Jittered sample-level points  

---

## **🚀 How to Run**

1. Make sure `gds_result_MZL.txt` exists in the working directory.
2. Open R or RStudio.
3. Run the script:

```r
source("run_MZL_pipeline.R")
```

Note: alternatively, the code can also be run piecemeal in the RStudios

4. View outputs in the working directory.

---

## **⚙ Optional Settings**

### **Whitelist Genes**
To keep only specific genes in each Excel sheet:

```r
GENE_WHITELIST <- c("ROR1","MS4A1","CD19")
```

### **Adjust Cleanup Frequency**

```r
clean_every = 1L
```

---

## **❗ Caveats & Notes**

- Disease classification depends on GEO annotation text; mislabeled samples may slip through.
- Some GPLs lack proper gene annotations; ROR1 may be missing for specific platforms.
- Series Matrix files are used; raw-level preprocessing is not performed.

---

## **📬 Contact**

Open a GitHub issue for feedback or enhancements.

---

## **✔ Summary**

This pipeline provides a **fully automated**, **platform-aware**, and **memory-efficient** approach to:

- Harvesting MZL samples from GEO  
- Normalizing heterogeneous datasets  
- Extracting and comparing ROR1 expression  
- Producing interpretable, publication-ready outputs  

Perfect for biomarker discovery and cross-study expression meta-analysis.

---
