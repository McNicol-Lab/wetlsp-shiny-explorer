read_evi_subset <- function(cfg, site_row, year, series = "spline", date_range = NULL, max_pixels = 250) {
  requireNamespace("arrow", quietly = TRUE)
  requireNamespace("dplyr", quietly = TRUE)

  ts_path <- resolve_timeseries_path(cfg, site_row)
  geom_path <- resolve_geom_path(cfg, site_row)

  # Open Arrow datasets. This works best locally, on mounted volumes, or with object storage.
  # Remote WebDAV directory listing may be unreliable; for production consider summary files or S3.
  ts_ds <- arrow::open_dataset(ts_path)
  geom_ds <- arrow::open_dataset(geom_path)

  pixel_ids <- ts_ds |>
    filter(.data$year == !!as.integer(year), .data$series %in% !!series) |>
    distinct(.data$pixel_id) |>
    head(max_pixels) |>
    collect() |>
    pull(.data$pixel_id)

  if (length(pixel_ids) == 0) {
    return(list(
      timeseries = tibble::tibble(),
      geometry = tibble::tibble(),
      paths = list(timeseries = ts_path, geometry = geom_path)
    ))
  }

  q <- ts_ds |>
    filter(
      .data$year == !!as.integer(year),
      .data$series %in% !!series,
      .data$pixel_id %in% !!pixel_ids
    )

  if (!is.null(date_range) && length(date_range) == 2 && all(!is.na(date_range))) {
    q <- q |>
      filter(
        .data$date >= as.Date(!!date_range[[1]]),
        .data$date <= as.Date(!!date_range[[2]])
      )
  }

  ts <- q |>
    select(.data$pixel_id, .data$series, .data$date, .data$year, .data$evi) |>
    collect() |>
    mutate(date = as.Date(.data$date))

  geom <- geom_ds |>
    filter(.data$pixel_id %in% !!pixel_ids) |>
    collect()

  geom <- geom_to_lonlat_if_possible(geom, cfg, site_row)

  list(
    timeseries = ts,
    geometry = geom,
    paths = list(timeseries = ts_path, geometry = geom_path)
  )
}

geom_to_lonlat_if_possible <- function(geom, cfg, site_row) {
  # If geometry already has lon/lat columns, use them.
  if (all(c("lon", "lat") %in% names(geom))) return(geom)

  # If only projected x/y are present, transform using CRS in meta.parquet if sf is available.
  if (!all(c("x", "y") %in% names(geom))) return(geom)

  if (!requireNamespace("sf", quietly = TRUE) || !requireNamespace("arrow", quietly = TRUE)) {
    # Fallback: use site location for all points so the map still renders.
    geom$lon <- site_row$lon[[1]]
    geom$lat <- site_row$lat[[1]]
    return(geom)
  }

  meta_path <- resolve_meta_path(cfg, site_row)

  crs_wkt <- tryCatch({
    meta <- arrow::read_parquet(meta_path)
    meta$value[meta$key == "crs_wkt"][[1]]
  }, error = function(e) NA_character_)

  if (is.na(crs_wkt) || !nzchar(crs_wkt)) {
    geom$lon <- site_row$lon[[1]]
    geom$lat <- site_row$lat[[1]]
    return(geom)
  }

  pts <- sf::st_as_sf(geom, coords = c("x", "y"), crs = crs_wkt, remove = FALSE)
  pts_ll <- sf::st_transform(pts, 4326)
  coords <- sf::st_coordinates(pts_ll)

  geom$lon <- coords[, 1]
  geom$lat <- coords[, 2]
  geom
}
