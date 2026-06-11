# WetLSP Demo Guide

This guide is intended for presentations, workshops, and demonstrations.

## Recommended Workflow

## Step 1: Cache Demonstration Sites

Before launching the application:

mkdir -p ~/wetlsp-cache
CYVERSE_ROOT="/data-store/iplant/home/shared/esiil/ai_for_natural_methane_working_group/wetlsp-cyverse"
cp "$CYVERSE_ROOT/wetlsp_cyverse_site_catalog_final.csv" ~/wetlsp-cache/
cp -r "$CYVERSE_ROOT/US-Atq" ~/wetlsp-cache/
cp -r "$CYVERSE_ROOT/PE-TNR" ~/wetlsp-cache/
cp -r "$CYVERSE_ROOT/CA-SCC" ~/wetlsp-cache/
cp -r "$CYVERSE_ROOT/BW-Nxr" ~/wetlsp-cache/
cp -r "$CYVERSE_ROOT/NZ-Kop" ~/wetlsp-cache/

## Step 2: Configure the App

data_root: "/home/rstudio/wetlsp-cache"
catalog_csv: "wetlsp_cyverse_site_catalog_final.csv"
mode: "local"

## Step 3: Launch

shiny::runApp()

## Demonstration Sequence

1. Overview tab
2. Site catalog
3. Phenometric visualizer
4. EVI time-series explorer
5. Pixel selection tools

## Recommended Settings

Maximum Pixels:

100–250

Year:

Single year at a time

## Avoid:

* Very large polygon selections
* Thousands of simultaneous pixels
* Uncached sites