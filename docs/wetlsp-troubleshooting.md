# WetLSP Troubleshooting

## Remote I/O Error

Example:

IOError: Could not open Parquet input source
Remote I/O error

Cause:

CyVerse Data Store latency during Parquet access.

Solution:

Cache the site locally:

cp -r "$CYVERSE_ROOT/SITE-ID" ~/wetlsp-cache/

⸻

## App Closes Unexpectedly

Possible causes:

* Excessive pixel selection
* Very large polygon selections
* Reading uncached Parquet datasets

Recommended actions:

* Reduce maximum pixel count
* Use cached sites
* Restrict date range

⸻

## Missing Basemap

Refresh the browser.

Alternatively switch basemap providers using the map controls.

⸻

## NetCDF Layers Not Loading

Verify:

ls SITE-ID/*.nc

Confirm that annual NetCDF products exist for the selected year.

⸻

## Site Not Available

Check:

* wetlsp_2021
* wetlsp_2022
* wetlsp_2023
* wetlsp_2024

within the site catalog.

A value of:

1

indicates availability.