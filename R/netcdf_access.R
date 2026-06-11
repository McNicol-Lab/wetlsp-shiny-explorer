safe_netcdf_layers <- function(path) {
  if (!requireNamespace("terra", quietly = TRUE)) {
    return(tibble::tibble(
      layer = NA_character_,
      note = "Install package 'terra' to inspect NetCDF layers."
    ))
  }

  if (is_url(path)) {
    # terra can sometimes read GDAL-supported URLs, but HTTPS/WebDAV NetCDF access is variable.
    # For production, cache selected NetCDFs locally or precompute layer summaries.
    return(tibble::tibble(
      layer = NA_character_,
      path = path,
      note = "Remote NetCDF layer inspection is disabled in this starter app. Cache locally or precompute summaries."
    ))
  }

  if (!file.exists(path)) {
    return(tibble::tibble(
      layer = NA_character_,
      path = path,
      note = "NetCDF file not found."
    ))
  }

  tryCatch({
    r <- terra::rast(path)
    tibble::tibble(
      layer_index = seq_len(terra::nlyr(r)),
      layer = names(r),
      nrow = terra::nrow(r),
      ncol = terra::ncol(r),
      xmin = terra::ext(r)$xmin,
      xmax = terra::ext(r)$xmax,
      ymin = terra::ext(r)$ymin,
      ymax = terra::ext(r)$ymax
    )
  }, error = function(e) {
    tibble::tibble(
      layer = NA_character_,
      path = path,
      note = paste("Could not read NetCDF:", conditionMessage(e))
    )
  })
}
