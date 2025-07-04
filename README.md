# Global Avian Influenza Risk: Waterbird entropy for targeted surveillance

The updated vertion is file 'Code of WAE (Version 202507).R'

## Overview
This repository presents a spatial modeling framework that evaluates the ecological drivers and predictive performance of Waterbird Activity Entropy (WAE) as a proxy for avian influenza virus (AIV) outbreak risk. Leveraging species distribution models (SDMs), global surveillance data, and human/cattle/poultry density, the project quantifies transmission potential across different regions/countries and functional bird groups.
- **1.SDMs modeling**: Constructs monthly niche models using filtered occurrence data and climate layers
- **2.SDMs validation**: Validation via True Skill Statistic (TSS)
- **3.Statistic for Main Figures**: 
  -**(1)Waterbird activity entropy (WAE) calculations**
  -**(2)Performance of WAE for predicting reported cases of AIV**
  -**(3)Hotspot identification**
  -**(4)Dominant waterbird functional groups identification**
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
- **BirdLife data zone**: BirdLife (https://datazone.birdlife.org/species/requestdis)
- **Environmental Layers**: 
  - Climate data (e.g., monthly climate variables from WorldClim).
  - Water percentage:`https://land.copernicus.eu/en`.
  - Human/poultry/cattle density rasters:`https://data.apps.fao.org/catalog/iso/15f8c56c-5499-45d5-bd89-59ef6c026704`.
- **AIVs Data**: Outbreak records (e.g., avian influenza cases,`https://empres-i.apps.fao.org/`).

---

## Workflow
### 1. SDMs modeling
- **Input**: Monthly Climate variables (temperature, precipitation), elevation, water percentage.
- **Models**: random forest (RF), maximum entropy (Maxentropy), and eXtreme Gradient Boosting (XGBoost)
- **Process**:
  - Sample background points for pseudo-absences.
  - Split data into training/testing sets.
  - Train models and save predictions (probability rasters).

### 2. SDMs validation
- **Process**: calibration scores .tif projections per species per month

### 3. Statistic for Main Figures
#### (1) Waterbird activity entropy (WAE) calculations
- **WAE calculation**: Entropy calculated using Shannon index across monthly species presence
- **Species richness and CV calculation**: Richness and variability assessed via species CV and latitudinal shifts
- **Visualazion**: Spatial patterns visualized in Figure 1a–1b
- *Output*: spNumTif.tif, sp_cv.tif, WAE.tif Fig.1a, Fig.1b
- 
#### (2) Predictive Performance of WAE
- **AIVs reports Filtering**: Filters outbreaks to exclude poultry–poultry transmission
- **AUC-ROC**: ROC curves assess entropy's predictive power (AUC ≥ 0.8)
- **Country-level statistic**: Country-level precision, recall, and F1 scores
- **Serotype statistic**: Serotype-specific sensitivity analysis
- *Output*: result_country.csv, Serotype_edit.csv Fig.2a–2d

#### (3) Hotspot Identification
- **Identification**: Identifies areas with co-occurrence of high-risk (using WAE thresholds) and high-density of human/cattle/poultry
- **Classification**: Categorical raster assigns composite hotspot typologies
- **Regional Statistic**: Regional summaries across regions/countries
- *Output*: hotentrpoppoulcat.tif, hot_statistic.csv, country_hot.csv Fig.3a, Fig.3b

#### (4) Functional Group Analysis
- **Activoty Entropy calculation**: Entropy calculated separately for bird functional groups (Waterfowl, Shorebirds, Seabirds, Waders, Others) and Two host definitions (Confirmed, All Suspected)
- **Person correlation**: Results correlated with human/cattle/poultry density using Person r
- **Person correlation**: Results correlated with human/cattle/poultry density using Person r
- 
- 
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



