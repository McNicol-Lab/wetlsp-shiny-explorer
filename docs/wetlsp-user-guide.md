# WetLSP User Guide

## Purpose

The WetLSP Explorer allows users to visualize vegetation dynamics and annual phenometrics for wetland flux tower sites.

## Site Catalog

The Overview tab contains:

* Site metadata
* Geographic coordinates
* Flux network information
* Annual product availability

## Phenometric Visualizer

The Phenometric Visualizer allows comparison of up to four annual NetCDF products simultaneously.

## Timing Metrics

Examples:

* OGI
* 50PCGI
* OGMx
* Peak
* OGD
* 50PCGD
* OGMn

These represent day-of-year timing metrics.

## Vegetation Metrics

Examples:

* EVImax
* EVIamp
* EVIarea

These describe vegetation magnitude and productivity.

## Pixel Time-Series Explorer

Users can:

* View randomly sampled pixels
* Select individual pixels
* Select groups of pixels using polygons
* Compare raw and spline-interpolated EVI trajectories

## Working with New Sites

To explore a new site efficiently:

cp -r "$CYVERSE_ROOT/SITE-ID" ~/wetlsp-cache/

Then reload the application.

## Data Products

NetCDF

Annual phenometrics.

Parquet

Pixel-level EVI trajectories.

JSON / CSV

Site catalog and metadata.

## Intended Uses

* Methane modeling
* Machine learning
* Phenological analysis
* Wetland ecology
* Flux tower interpretation