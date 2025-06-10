![LGG_mIHC](https://img.shields.io/badge/Project-LGG_mIHC-green?style=for-the-badge&logo=R&logoColor=white)

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)  
[![Last Commit](https://img.shields.io/github/last-commit/3DBMandNE-Lab/LGG_mIHC?color=blue)](https://github.com/3DBMandNE-Lab/LGG_mIHC/commits/main)  
[![R](https://img.shields.io/badge/Made%20with-R-1f425f.svg?logo=R&logoColor=white)](https://cran.r-project.org/)  
[![Issues](https://img.shields.io/github/issues/3DBMandNE-Lab/LGG_mIHC)](https://github.com/3DBMandNE-Lab/LGG_mIHC/issues)

# LGG_mIHC Pipeline

[![Email](https://img.shields.io/badge/Email-kevin.joseph@uniklinik--freiburg.de-blue?logo=gmail)](mailto:kevin.joseph@uniklinik-freiburg.de)  
[![ORCID](https://img.shields.io/badge/ORCID-0000--0003--0062--2099-a6ce39?logo=orcid)](https://orcid.org/0000-0003-0062-2099)  

---

## Overview

**LGG_mIHC** is an automated R-based pipeline to quantify and dissect **spatial network interactions** of immune and tumor cells in **low-grade gliomas (LGG)** using **multiplex immunohistochemistry (mIHC)** data.

The pipeline processes raw CSV exports of **cell-level marker intensities** to produce:

✅ Cell data annotated by **hypoxia** and **IDH status**  
✅ **Sub-region (ROI) analyses** — spatial distance metrics & networks  
✅ **Whole-image global network analyses**  
✅ Ready-to-use **RDS objects**, **publication-ready PDFs**, and **summary tables**

---

## Repository Structure

```text
LGG_mIHC/
├── DataPrep.R         # Data loading, annotation, merging, RDS export
├── ROIAnalysis.R      # Sub-region (ROI) spatial distance & network analysis
├── WholeImage.R       # Global network analysis on entire TMA images
├── outputs/           # Generated RDS, PDFs, figures, and summary tables
├── LICENSE            # MIT License terms
└── README.md          # This landing page
```

---

## Prerequisites

- **R** version ≥ 4.0  
- Required R packages:

```r
install.packages(c(
  "tidyverse", "readxl", "reshape2", "ggnetwork", "igraph",
  "ggraph", "tidygraph", "circlize", "ggforce", "pheatmap",
  "spatstat", "EBImage", "rmarkdown"
))
```

---

## Quick Start

### Clone this repository

```bash
git clone https://github.com/3DBMandNE-Lab/LGG_mIHC.git
cd LGG_mIHC
```

### Configure file paths

Edit the first lines of each R script to point to your input `data/` and desired `outputs/` directories.

### Run the pipeline

```bash
Rscript DataPrep.R
Rscript ROIAnalysis.R
Rscript WholeImage.R
```

### Explore results

Results will appear in the `outputs/` folder:

- `processed_data.rds` — merged and annotated cell-level data  
- `ROI_networks.pdf`, `ROI_violin_plots.pdf`  
- `WholeImage_networks.pdf`, `WholeImage_violin_plots.pdf`  
- Summary CSVs: `roi_metrics.csv`, `whole_image_metrics.csv`

---

## Script Details

### DataPrep.R

**Purpose:**  
- Import raw CSVs of mIHC marker intensities  
- Annotate with hypoxia status and IDH status  
- Save merged dataset as `processed_data.rds`

**Key steps:**  
- `list.files()` → `read_csv()` → `bind_rows()`  
- Add `Parent` column labels → Hypoxic / Normoxic  
- Save RDS for downstream use.

---

### ROIAnalysis.R

**Purpose:**  
- Divide each image into **ROIs** (default: 150 × 150 µm)  
- Compute **pairwise distances** between selected cell types  
- Aggregate statistics across ROIs  
- Build **ROI-level spatial networks**

**Outputs:**  
- `ROI_networks.pdf` → network plots for each condition  
- `ROI_violin_plots.pdf` → violin/boxplots of distance distributions  
- `roi_metrics.csv` → distance metrics per ROI

---

### WholeImage.R

**Purpose:**  
- Compute **global cell–cell distance networks** for whole TMA images  
- Analyze **global interaction patterns** across hypoxia and IDH groups

**Outputs:**  
- `WholeImage_networks.pdf` → whole-image network layouts  
- `WholeImage_violin_plots.pdf` → violin/boxplots of distances  
- `whole_image_metrics.csv` → global summary statistics

---

## Outputs & Visualization

The pipeline produces:

- **RDS** → Load in R or Shiny for further analysis  
- **PDFs** → Network graphs and violin/boxplots suitable for publication  
- **CSVs** → Numeric summary tables for statistical modeling or reporting

---

**Note:** Always re-run `DataPrep.R` before submitting a PR to ensure consistent RDS output.

---

## 📜 License

Distributed under the **MIT License**. See [LICENSE](LICENSE) for details.

---

## 📫 Contact

**Dr.-Ing. Kevin Joseph**
Department of Neurosurgery, University Hospital Freiburg  
[kevin.joseph@uniklinik-freiburg.de](mailto:kevin.joseph@uniklinik-freiburg.de)
