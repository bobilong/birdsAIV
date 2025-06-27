# GBIF Data Processing and Monthly Species Distribution Modeling (SDMs)

The final vertion is in file 'WAE_bird'

## Overview
This repository contains scripts for processing GBIF (Global Biodiversity Information Facility) data, performing species distribution modeling (SDMs), and analyzing ecological patterns. The workflow includes data cleaning, spatial filtering, model training/validation, and visualization. Key steps include:
- **Data Preprocessing**: Filtering species records and removing outliers.
- **Spatial Analysis**: Aligning data with IUCN boundaries and environmental variables.
- **Modeling**: Training machine learning models (e.g., Random Forest, XGBoost) for SDMs.
- **Validation**: Evaluating model performance using AUC-ROC and waterbird activity entropy (WAE).
- **Visualization**: Generating maps for species richness, AE, and risk hotspots.

---

## Requirements
### R Libraries
- `data.table`, `terra`, `tidyterra`, `stringr`
- `ggplot2`, `dplyr`, `rnaturalearth`, `paletteer`
- `vegan`, `linkET` (for Mantel tests)
- Full list in `library()` calls in the code.


### External Data
- **GBIF Data**: CSV files stored in `[https://figshare.com/articles/dataset/_b_Enhancing_Global_Avian_Influenza_Exposure_Risk_Assessment_with_Waterbird_Activity_Entropy_b_/28504868]`.
- **IUCN Boundaries**: Shapefiles in `[/root/autodl-tmp/IUCN_Data/](https://datazone.birdlife.org/species/requestdis]`.
- **
- **Environmental Layers**: 
  - Climate data (e.g., monthly climate variables from WorldClim).
  - Water percentage:`https://land.copernicus.eu/en`.
  - Human/poultry/cattle density rasters:`https://data.apps.fao.org/catalog/iso/15f8c56c-5499-45d5-bd89-59ef6c026704`.
- **Validation Data**: Outbreak records (e.g., avian influenza cases,`https://empres-i.apps.fao.org/`).

---

## Workflow
### 1. Data Preparation
- **Input**: GBIF CSV files, IUCN boundaries, climate rasters.
- **Steps**:
  - Filter records by observation type (`HUMAN_OBSERVATION`, `MACHINE_OBSERVATION`).
  - Remove species with fewer than 108 records.
  - Clip points to IUCN species ranges and thin spatially to reduce autocorrelation.

### 2. Model Training
- **Features**: Climate variables (temperature, precipitation), elevation, water percentage.
- **Models**: Random Forest, XGBoost, Logistic Regression (MaxEnt equivalent).
- **Process**:
  - Sample background points for pseudo-absences.
  - Split data into training/testing sets.
  - Train models and save predictions (probability rasters).

### 3. Validation
- **AUC-ROC**: Evaluate model discrimination ability.
- **Waterbird Activity Entropy (WAE)**: Quantify spatiotemporal aggregation of species records.
- **Hotspot Analysis**: Identify high-risk regions using AE thresholds and environmental covariates.

### 4. Visualization
- **Maps**: Species richness, WAE, and risk hotspots.
- **Statistical Plots**: Boxplots for model comparisons, ROC curves, and Mantel tests for correlation analysis.


---

## Usage
1. **Environment Setup**:
   - Install required R/Python packages.
   - Update file paths in the code (e.g., `/root/autodl-tmp/` to your local paths).

2. **Run Workflow**:
   - Execute R sections for data preprocessing, modeling, and visualization.
   - Run Python blocks for SDM simulations (ensure `pyimpute` and `sklearnex` are installed).

3. **Outputs**:
   - Model predictions.
   - Validation metrics (AUC scores, WAE maps).
   - Figures (PDF/PNG) for ecological patterns and risk hotspots.

---

## Notes
- **Data Paths**: Modify hardcoded paths (e.g., `/root/autodl-tmp/`) to match your local environment.
- **Parallelization**: Use `n_jobs=-1` in Python sections for multi-core processing.
- **Citation**: Ensure proper attribution for GBIF, IUCN, and WorldClim data sources.



