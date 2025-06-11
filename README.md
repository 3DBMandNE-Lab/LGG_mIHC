![R](https://img.shields.io/badge/R-4.4.0-blue)
![Seurat](https://img.shields.io/badge/Seurat-5.0.3-blueviolet)
![License](https://img.shields.io/badge/license-MIT-green)
![Status](https://img.shields.io/badge/status-Research--Collaboration-yellow)

# Spatial Network Interactions in GBM (mIHC imaging based)

---

## Overview

**LGG_mIHC** is an automated R-based pipeline to quantify and dissect **spatial network interactions** of immune and tumor cells in **low-grade** and **high grade** gliomas using **multiplex immunohistochemistry (mIHC)** data obtained from tissue microarrays.

The pipeline processes raw CSV exports of **cell-level marker intensities** to produce:

- Cell data annotated by **hypoxia**, and **cannonical cellular markers**  
- **ROI analyses** — Spatial Distance Metrics & Networks  
- **Whole-image global network analyses**  


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

## Quick Start

### Clone this repository

```bash
git clone https://github.com/3DBMandNE-Lab/LGG_mIHC.git
cd LGG_mIHC
```

### Setup

Modify each script to point to your input `data/` and desired `outputs/` directories.

### Explore results

Results will appear in the `outputs/` folder:

- `processed_data.rds` — merged and annotated cell-level data  
- `ROI_networks.pdf`, `ROI_violin_plots.pdf`  
- `WholeImage_networks.pdf`, `WholeImage_violin_plots.pdf`  
- Summary CSVs: `roi_metrics.csv`, `whole_image_metrics.csv`

---

## Script Details

### DataPrep.R

**Features:**  
- Import raw CSVs of mIHC marker intensities  
- Annotate with hypoxia status and IDH status  
- Save merged dataset as `processed_data.rds`

**Key steps:**  
- `list.files()` → `read_csv()` → `bind_rows()`  
- Add `Parent` column labels → Hypoxic / Normoxic  
- Save RDS for downstream use.

---

### ROIAnalysis.R

**Features:**  
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

**Features:**  
- Compute **global cell–cell distance networks** for whole TMA images  
- Analyze **global interaction patterns** across hypoxia and IDH groups

**Outputs:**  
- `WholeImage_networks.pdf` → whole-image network layouts  
- `WholeImage_violin_plots.pdf` → violin/boxplots of distances  
- `whole_image_metrics.csv` → global summary statistics

---

## Author

**Kevin Joseph**  
3DBM & NE, Neurosurgery, University Hospital Freiburg  
Germany  

---

