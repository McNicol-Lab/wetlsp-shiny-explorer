# WetLSP CyVerse Explorer

Starter Shiny app for exploring WetLSP PlanetScope EVI products stored in a portable folder layout:

```text
wetlsp-cyverse/
    wetlsp_cyverse_site_catalog_final.csv
    wetlsp_cyverse_site_catalog_final.json
    SITE1/
        SITE1-wetlsp-2021.nc
        SITE1-wetlsp-2022.nc
        SITE1-wetlsp-2023.nc
        SITE1-wetlsp-2024.nc
        pixels_geom_ds/
        pixels_meta_ds/
        pixels_timeseries_ds/
        README-parquet.json
        README-parquet.md
    SITE2/
    ...
```

The app is designed to work with either:

1. A local or mounted data directory, e.g. `/Volumes/CyVerse/wetlsp-cyverse`
2. A public HTTPS/WebDAV root, e.g. `https://data.cyverse.org/dav-anon/iplant/home/<user>/wetlsp-cyverse`
3. A future S3-compatible bucket, if paths are adapted.

## Core functions

- Site catalog browser
- Interactive site map
- EVI time-series plotting from Parquet
- Raw versus spline filtering
- Optional pixel subsampling
- Annual NetCDF phenometric metadata/layer inspection
- Direct data-path display for downloads/citation

## Install R packages

```r
install.packages(c(
  "shiny", "bslib", "DT", "dplyr", "readr", "ggplot2", "leaflet",
  "arrow", "duckdb", "DBI", "glue", "stringr", "lubridate",
  "tidyr", "jsonlite", "ncdf4", "terra"
))
```

`terra` and `ncdf4` are only needed for NetCDF layer inspection. The Parquet time-series explorer mainly uses `arrow`, `dplyr`, and `ggplot2`.

## Configure the data root

Copy `config.example.yml` to `config.yml` and edit:

```yaml
data_root: "/path/to/wetlsp-cyverse"
catalog_csv: "wetlsp_cyverse_site_catalog_final.csv"
mode: "local"
```

For public CyVerse WebDAV/HTTPS, use a URL root:

```yaml
data_root: "https://data.cyverse.org/dav/iplant/projects/esiil/ai_for_natural_methane_working_group/wetlsp-cyverse"
catalog_csv: "wetlsp_cyverse_site_catalog_final.csv"
mode: "http"
```

Important: full remote Parquet directory reads over WebDAV/HTTPS can be slow or unsupported depending on server listing behavior. For production, either:
- precompute site/year summary Parquet files, or
- serve from S3-compatible object storage, or
- run the app inside CyVerse/VICE near the data.

## Run

```r
shiny::runApp()
```

## Recommended production enhancement

Before deploying to public users, create summary tables such as:

```text
summaries/
    site_daily_evi_summary.parquet
    site_year_phenometric_summary.parquet
    site_pixel_sample_index.parquet
```

The app can load summaries instantly, then query pixel-level Parquet only when the user requests detailed subsets.
