# The Arctic warming-driven host convergence bridges continents for avian influenza spreading

**Version**: Code of Arctic_AIVhost / CHEC Analysis (Version 202502)

## Overview
This repository implements a spatial and trait-based modeling framework to evaluate how climate change reshapes cross-continental host encounter combinations (CHEC) and the potential for avian influenza (AIV) viral sharing across the Arctic.
We project present and future summer distributions (May–August) of  confirmed AIV host species and quantify how both species richness and cross‑continental encounter opportunities change under multiple climate scenarios. Specifically, the workflow includes:

- **1.SDMs Construction**: Monthly niche modeling with occurrence and environment data and Evaluation using True Skill Statistic (TSS)
- **2.Arctic host distribution projections**: Present and future (to end of century) summer distributions of 415 confirmed AIV hosts 
- **3.Statistical Framework for Core Figures**: 
  - **(1)Species richness trends under climate scenarios**
  - **(2)Cross‑continental host encounter combinations (CHEC) and Future changes in CHEC (ΔCHEC)**
  - **(3)Sector-based trans‑Arctic comparison**
  - **(4)Trait-based drivers of CHEC**
---

## Requirements
### Core R packages used in this workflow include:
terra, data.table, dplyr, stringr
ggplot2, paletteer, patchwork
rnaturalearth, grfxtools
ggtree, ggtreeExtra, ggnewscale
brms / rstanarm / tidybayes (for Bayesian GLMs and posterior summaries)


### External Data
- **GBIF Data**: CSV files stored in `[https://figshare.com/articles/dataset/_b_Enhancing_Global_Avian_Influenza_Exposure_Risk_Assessment_with_Waterbird_Activity_Entropy_b_/28504868]`.
- **IUCN Boundaries**: Shapefiles in `[/root/autodl-tmp/IUCN_Data/](https://datazone.birdlife.org/species/requestdis]`.
- **BirdLife data zone**: BirdLife (https://datazone.birdlife.org/species/requestdis)
- **Environmental Layers**: 
  - Present Climate data (e.g.monthly climate variables from WorldClim).
  - Future Climate data (e.g.Future projections under combinations of RCPs and SSPs).
- **Trait data**:
  - Hand‑wing index (HWI)
  - Generation length
  - Diet and foraging traits
  - Functional group assignments
  (https://doi.org/10.6084/m9.figshare.27051040.v1)

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
#### (1) Arctic host distribution and richness
- **Richness calculation**: Present and future summer distributions (May–August) of confirmed AIV hosts
- **Hotspot identification**: Richness hotspots identified using quantile thresholds (Top 5%, Top 10%)
- **Future richness trends**: Linear regression of richness under three climate scenarios (RCP–SSP combinations)
- **Latitudinal summaries**: Richness change summarized across Tropics, Temperate, and Arctic zones
- **Visualization**: Polar richness maps, slope significance maps, and latitudinal percentage barplots
- *Output*: presentspNumTif.tif, slope rasters, percent-change tables, Fig.1
  
#### (2) CHEC quantification and bivariate mapping
- **CHEC definition**: CHEC = cross‑continental host encounter combinations, calculated as the diversity of species‑pair combinations between Eurasian and American hosts
- **Present CHEC**: CHEC computed from present distributions and visualized using polar projection
- **American vs Eurasian continental host comparison**: Richness normalized and combined into a 10×10 bivariate color scale
- **ΔCHEC calculation**:ΔCHEC = CHEC_future − CHEC_present for SSP1‑2.6, SSP3‑7.0, SSP5‑8.5, ΔCHEC categorized into 10 positive and 10 negative bins using custom diverging color scales
- **Regional statistic**:ΔCHEC summarized by: Arctic, Temperate, Tropical zones // American vs Eurasian continental host // SSP × Year combinations
- **Visualization**: Bivariate richness map and CHEC hotspot map, Polar ΔCHEC maps and stacked/mirrored barplots showing hemispheric contributions
- *Output*: presentCHEC_2newold.tif, bivariate map, ΔCHEC rasters,stats_deltaCHECArctic.csv, Fig.2

#### (3) Sector‑based comparison: Cross‑Pacific vs Cross‑Atlantic
- **Sector construction**: Two Arctic migratory corridors defined using custom polygon sectors:Cross‑Pacific sector (center 180°), Cross‑Atlantic sector (center −30°)
- **CHEC extraction**: Present and future CHEC values extracted within each sector
- **Relative change**: Comparison of ΔCHEC between sectors to identify emerging trans‑Arctic encounter pathways
- **Visualization**: Polar CHEC maps with sector boundaries and sector‑specific summaries
- *Output*: Sector shapefiles, sector CHEC tables, Fig.3

#### (4) Trait‑based modeling of CHEC (Bayesian GLM)
- **Goal**: Identify ecological and life‑history traits that predict cross‑continental encounter combinations (CHEC)
- **Traits included**:
   - Hand‑wing index (HWI)
   - Generation length
   - Diet categories
   - Foraging traits
   - Functional group (e.g., waterfowl, shorebirds, seabirds)
- **Model**: Bayesian generalized linear model fitted to CHEC values using posterior sampling
- **Visualization**: Trait effect forest plot showing posterior distributions and uncertainty
- *Output*: Trait posterior tables, effect summaries, Fig.4

  
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
   - Validation metrics (AUC scores).
   - Figures (PDF/PNG) for present, future change of richness and CHEC patterns.

---

## Notes
- **Data Paths**: Modify hardcoded paths (e.g., `/root/autodl-tmp/`) to match your local environment.
- **Parallelization**: Use `n_jobs=-1` in Python sections for multi-core processing.
- **Citation**: Ensure proper attribution for GBIF, IUCN, and WorldClim data sources.



