# Global Avian Influenza Risk: Waterbird entropy for targeted surveillance

**Version**: Code of WAE (Version 202507).R 

## Overview
This repository presents a spatial modeling framework that evaluates the ecological drivers and predictive performance of Waterbird Activity Entropy (WAE) as a proxy for avian influenza virus (AIV) outbreak risk. Leveraging species distribution models (SDMs), global surveillance data, and human/cattle/poultry density, the project quantifies transmission potential across different regions/countries and functional bird groups.
- **1.SDMs Construction**: Monthly niche modeling with occurrence and environment data
- **2.SDMs Validation**: Evaluation using True Skill Statistic (TSS)
- **3.Statistical Framework for Core Figures**: 
  -**(1)Waterbird activity entropy (WAE) raster generation**
  -**(2)AIV cases prediction performance**
  -**(3)Hotspot identification**
  -**(4)Functional group prioritization and correlation analysis**
---

## Requirements
### R Libraries
- `data.table`, `terra`, `tidyterra`, `stringr`
- `ggplot2`, `dplyr`, `rnaturalearth`, `paletteer`
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
### 1. Species Distribution Modeling
- **Input**: Filtered monthly occurrences; Monthly Climate variables (temperature, precipitation), elevation, water percentage.
- **Models**: random forest (RF), maximum entropy (Maxentropy), and eXtreme Gradient Boosting (XGBoost)
- **Process**:
  - Sample background points for pseudo-absences.
  - Split data into training/testing sets.
  - Train models and save predictions (probability rasters).

### 2. SDMs Validation
- **Process**: Model calibration assessed via True Skill Statistic (TSS); Summary metrics compiled across species/months
- *Output*: SDMs_TSS.csv

### 3. Statistical Framework for Core Figures
#### (1) WAE raster generation
- **WAE calculation**: Entropy calculated using Shannon index across monthly species presence
- **Species richness and CV calculation**: Richness and variability assessed via species CV and latitudinal shifts
- **Visualazion**: Spatial patterns visualized in Figure 1a–1b
- *Output*: spNumTif.tif, sp_cv.tif, WAE.tif Fig.1a, Fig.1b
- 
#### (2) AIV cases prediction performance
- **AIVs reports Filtering**: Filters outbreaks to exclude poultry–poultry transmission
- **AUC-ROC**: ROC curves assess entropy's predictive power (AUC ≥ 0.8)
- **Country-level statistic**: Country-level precision, recall, and F1 scores
- **Serotype statistic**: Serotype-specific sensitivity analysis
- *Output*: result_country.csv, Serotype_edit.csv Fig.2a–2d

#### (3) Hotspot identification
- **Identification**: Identifies areas with co-occurrence of high-risk (using WAE thresholds) and high-density of human/cattle/poultry
- **Classification**: Categorical raster assigns composite hotspot typologies
- **Regional Statistic**: Regional summaries across regions/countries
- *Output*: hotentrpoppoulcat.tif, hot_statistic.csv, country_hot.csv Fig.3a, Fig.3b

#### (4) Functional Group Analysis
- **Activoty Entropy calculation**: Entropy calculated separately for bird functional groups (Waterfowl, Shorebirds, Seabirds, Waders, Others) and Two host definitions (Confirmed Hosts, All Suspected Hosts)
- **Person correlation**: - WAE vs. human/cattle/poultry densities (among key exposure regions/countries)
- **Spatial autocorrelation**: Bivariate Moran’s I evaluates spatial autocorrelation (among key exposure regions/countries)
- *Output*: Global_CHEntropy_df.csv, Global_birdEntropy_df.csv All5_poppoulcatresults_df.csv, Moran_Func.csv

  
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



