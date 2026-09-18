# Phosphoproteomics analysis of HCT-116 cells treated with tumor necrosis factor (TNF)

R/Quarto analysis accompanying:

**LUBAC PUB domain interactions restrict Met1-linked ubiquitination to prevent embryonic lethality and immune pathology in mice**

This repository contains the phosphoproteomics analysis of: HOIP KO HCT-116 human colon carcinoma cells, reconstituted with HOIP WT or with HOIP(N102D), following TNF stimulation for 5 and 15 minutes.

## Data

The phosphoproteomics data are deposited in PRIDE.

**PRIDE accession: PXD060649**

The analysis starts from:

`20241126_Report_PTM_pivot_GF (Pivot).tsv`

Place this file in the `data/` folder before running the analysis.

## Analysis

The complete analysis is contained in `analysis.qmd` and includes:

- data filtering and quality control
- log2 transformation and median centering
- differential phosphorylation analysis with limma
- regulated phosphosite summaries and overlap analysis
- phosphosite heatmaps
- KEGG TNF and NF-kB pathway heatmap

Helper functions used by the analysis are stored in `R/functions.R`.

## Running the analysis

Open `TNF_project.Rproj` in RStudio and render `analysis.qmd`.

Required R packages are listed at the beginning of the analysis.

Input data and generated outputs are kept locally and are excluded from the Git repository.

## Repository structure

``` text
analysis.qmd        Main phosphoproteomics analysis
R/functions.R       Helper functions
TNF_project.Rproj   RStudio project
data/               Input data (not tracked by Git)
outputs/            Generated tables and figures (not tracked by Git)
```
