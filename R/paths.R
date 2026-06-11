is_url <- function(x) {
  grepl("^https?://", x)
}

path_join_portable <- function(...) {
  parts <- list(...)
  parts <- unlist(parts, use.names = FALSE)
  parts <- parts[!is.na(parts) & nzchar(parts)]

  if (length(parts) == 0) return("")

  first <- parts[[1]]
  if (is_url(first)) {
    first <- sub("/+$", "", first)
    rest <- parts[-1]
    rest <- gsub("^/+|/+$", "", rest)
    paste(c(first, rest), collapse = "/")
  } else {
    do.call(file.path, as.list(parts))
  }
}

resolve_catalog_path <- function(cfg) {
  path_join_portable(cfg$data_root, cfg$catalog_csv)
}

resolve_site_path <- function(cfg, site_row, rel_col) {
  rel <- site_row[[rel_col]][[1]]
  path_join_portable(cfg$data_root, rel)
}

resolve_timeseries_path <- function(cfg, site_row) {
  resolve_site_path(cfg, site_row, "parquet_timeseries_path")
}

resolve_geom_path <- function(cfg, site_row) {
  resolve_site_path(cfg, site_row, "parquet_geom_path")
}

resolve_meta_path <- function(cfg, site_row) {
  resolve_site_path(cfg, site_row, "parquet_meta_path")
}

resolve_netcdf_path <- function(cfg, site_row, year) {
  col <- paste0("netcdf_", year, "_file")
  if (!col %in% names(site_row)) {
    stop("Catalog does not contain column: ", col)
  }
  path_join_portable(cfg$data_root, site_row[[col]][[1]])
}
