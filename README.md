# GBIF Data Processing and Species Distribution Modeling (SDM)

## Overview
This repository contains scripts for processing GBIF (Global Biodiversity Information Facility) data, performing species distribution modeling (SDM), and analyzing ecological patterns. The workflow includes data cleaning, spatial filtering, model training/validation, and visualization. Key steps include:
- **Data Preprocessing**: Filtering species records and removing outliers.
- **Spatial Analysis**: Aligning data with IUCN boundaries and environmental variables.
- **Modeling**: Training machine learning models (e.g., Random Forest, XGBoost) for SDM.
- **Validation**: Evaluating model performance using AUC-ROC and activity entropy (AE).
- **Visualization**: Generating maps for species richness, AE, and risk hotspots.

---

## Requirements
### R Libraries
- `data.table`, `terra`, `tidyterra`, `stringr`
- `ggplot2`, `dplyr`, `rnaturalearth`, `paletteer`
- `vegan`, `linkET` (for Mantel tests)
- Full list in `library()` calls in the code.

### Python Libraries (for SDM simulation)
- `elapid`, `sklearn`, `xgboost`, `lightgbm`, `geopandas`

### External Data
- **GBIF Data**: CSV files stored in `/root/autodl-tmp/GBIF_Data/`.
- **IUCN Boundaries**: Shapefiles in `/root/autodl-tmp/IUCN_Data/`.
- **Environmental Layers**: 
  - Climate data (e.g., BIO variables from WorldClim).
  - Elevation, water percentage, and human/poultry density rasters.
- **Validation Data**: Outbreak records (e.g., avian influenza cases).

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
- **Activity Entropy (AE)**: Quantify spatiotemporal aggregation of species records.
- **Hotspot Analysis**: Identify high-risk regions using AE thresholds and environmental covariates.

### 4. Visualization
- **Maps**: Species richness, AE, and risk hotspots.
- **Statistical Plots**: Boxplots for model comparisons, ROC curves, and Mantel tests for correlation analysis.

---

## File Structure
