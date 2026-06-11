# WetLSP Explorer

Interactive visualization and exploration platform for the WetLSP (Wetland Land Surface Phenology) dataset.

## Overview

WetLSP provides global high-resolution wetland vegetation dynamics derived from PlanetScope imagery for methane science, ecological forecasting, and machine learning applications.

The dataset contains:  

  * 95 wetland flux tower sites worldwide
  * Daily gap-filled EVI trajectories
  * 3 m spatial resolution vegetation dynamics
  * 24 annual phenometric products following Moon et al. (2022)
  * Pixel-level and site-level products
  * Interactive visualization through a Shiny application

## Data Products

### Pixel-Level Time Series

Daily EVI trajectories for individual pixels.

Formats:

* Apache Parquet
* R objects (.rds)

### Annual Phenometric Products

Twenty-four annual vegetation phenometrics generated using the Moon et al. (2022) methodology.

Format:

* NetCDF

### Site Catalog

Metadata describing all WetLSP sites.

Formats:

* CSV
* JSON

### Data Access

Data are hosted through CyVerse under:

`/data-store/iplant/home/shared/esiil/ai_for_natural_methane_working_group/wetlsp-cyverse`

### Launching the Application

shiny::runApp()

See:

* docs/wetlsp_demo_guide.md
* docs/wetlsp_user_guide.md

## Citation

McNicol et al.

WetLSP: High-Resolution Global Wetland Land Surface Phenology Derived from PlanetScope Imagery.
