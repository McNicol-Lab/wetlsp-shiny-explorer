load_site_catalog <- function(cfg) {
  path <- resolve_catalog_path(cfg)

  if (is_url(path)) {
    catalog <- readr::read_csv(path, show_col_types = FALSE)
  } else {
    if (!file.exists(path)) {
      stop(
        "Catalog not found: ", path, "\n",
        "Set data_root and catalog_csv in config.yml."
      )
    }
    catalog <- readr::read_csv(path, show_col_types = FALSE)
  }

  required <- c(
    "site_id", "site_name", "country", "lat", "lon",
    "parquet_timeseries_path", "parquet_geom_path",
    "netcdf_2021_file", "netcdf_2022_file", "netcdf_2023_file", "netcdf_2024_file"
  )

  missing <- setdiff(required, names(catalog))
  if (length(missing) > 0) {
    stop("Catalog is missing required columns: ", paste(missing, collapse = ", "))
  }

  catalog |>
    mutate(
      lat = as.numeric(.data$lat),
      lon = as.numeric(.data$lon),
      across(matches("^(base|wetlsp)_20[0-9]{2}$"), as.integer)
    ) |>
    arrange(.data$site_id)
}

validate_site_year_available <- function(site_row, year) {
  col <- paste0("wetlsp_", year)
  if (col %in% names(site_row)) {
    val <- site_row[[col]][[1]]
    if (!is.na(val) && as.integer(val) != 1L) {
      warning("Catalog indicates WetLSP product may not be available for ", site_row$site_id, " in ", year)
    }
  }
  invisible(TRUE)
}
