safe_read_parquet <- function(file, cols = NULL) {
  tryCatch({
    if (is.null(cols)) {
      arrow::read_parquet(file)
    } else {
      arrow::read_parquet(file, col_select = tidyselect::all_of(cols))
    }
  }, error = function(e) {
    warning(
      "Skipping unreadable parquet file: ",
      basename(file),
      " | ",
      conditionMessage(e)
    )
    NULL
  })
}

safe_list_parquet_files <- function(path) {
  files <- list.files(path, pattern = "\\.parquet$", full.names = TRUE)
  
  if (length(files) == 0) {
    warning("No parquet files found in: ", path)
  }
  
  files
}

read_evi_subset <- function(
    cfg,
    site_row,
    year,
    series = "spline",
    date_range = NULL,
    max_pixels = 250
) {
  requireNamespace("arrow", quietly = TRUE)
  requireNamespace("dplyr", quietly = TRUE)
  requireNamespace("tibble", quietly = TRUE)
  
  ts_path <- resolve_timeseries_path(cfg, site_row)
  geom_path <- resolve_geom_path(cfg, site_row)
  
  ts_files <- safe_list_parquet_files(ts_path)
  
  if (length(ts_files) == 0) {
    return(list(
      timeseries = tibble::tibble(),
      geometry = tibble::tibble(),
      paths = list(timeseries = ts_path, geometry = geom_path),
      skipped_files = character()
    ))
  }
  
  skipped_files <- character()
  
  # First pass: identify pixel IDs without requiring every parquet file to open.
  pixel_id_chunks <- lapply(ts_files, function(f) {
    dat <- tryCatch({
      arrow::read_parquet(
        f,
        col_select = tidyselect::all_of(c("pixel_id", "series", "year"))
      )
    }, error = function(e) {
      skipped_files <<- unique(c(skipped_files, f))
      warning(
        "Skipping unreadable parquet file during pixel-id scan: ",
        basename(f),
        " | ",
        conditionMessage(e)
      )
      return(NULL)
    })
    
    if (is.null(dat) || nrow(dat) == 0) return(NULL)
    
    dat |>
      dplyr::filter(
        .data$year == !!as.integer(year),
        .data$series %in% !!series
      ) |>
      dplyr::distinct(.data$pixel_id)
  })
  
  pixel_ids <- dplyr::bind_rows(pixel_id_chunks) |>
    dplyr::distinct(.data$pixel_id) |>
    head(max_pixels) |>
    dplyr::pull(.data$pixel_id)
  
  if (length(pixel_ids) == 0) {
    return(list(
      timeseries = tibble::tibble(),
      geometry = tibble::tibble(),
      paths = list(timeseries = ts_path, geometry = geom_path),
      skipped_files = basename(skipped_files)
    ))
  }
  
  # Second pass: read only rows needed for selected/sampled pixels.
  ts_chunks <- lapply(ts_files, function(f) {
    if (f %in% skipped_files) return(NULL)
    
    dat <- tryCatch({
      arrow::read_parquet(
        f,
        col_select = tidyselect::all_of(c("pixel_id", "series", "date", "year", "evi"))
      )
    }, error = function(e) {
      skipped_files <<- unique(c(skipped_files, f))
      warning(
        "Skipping unreadable parquet file during EVI read: ",
        basename(f),
        " | ",
        conditionMessage(e)
      )
      return(NULL)
    })
    
    if (is.null(dat) || nrow(dat) == 0) return(NULL)
    
    dat <- dat |>
      dplyr::filter(
        .data$year == !!as.integer(year),
        .data$series %in% !!series,
        .data$pixel_id %in% !!pixel_ids
      )
    
    if (!is.null(date_range) && length(date_range) == 2 && all(!is.na(date_range))) {
      dat <- dat |>
        dplyr::filter(
          .data$date >= as.Date(date_range[[1]]),
          .data$date <= as.Date(date_range[[2]])
        )
    }
    
    dat
  })
  
  ts <- dplyr::bind_rows(ts_chunks) |>
    dplyr::mutate(date = as.Date(.data$date))
  
  # Geometry is usually much smaller and less likely to fail, but still read safely.
  geom <- tryCatch({
    arrow::open_dataset(geom_path) |>
      dplyr::filter(.data$pixel_id %in% !!pixel_ids) |>
      dplyr::collect()
  }, error = function(e) {
    warning("Could not read geometry dataset: ", conditionMessage(e))
    tibble::tibble()
  })
  
  if (nrow(geom) > 0) {
    geom <- geom_to_lonlat_if_possible(geom, cfg, site_row)
  }
  
  list(
    timeseries = ts,
    geometry = geom,
    paths = list(timeseries = ts_path, geometry = geom_path),
    skipped_files = basename(skipped_files)
  )
}