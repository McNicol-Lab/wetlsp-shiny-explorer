`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || is.na(x)) y else x
}

read_app_config <- function(path = "config.yml") {
  defaults <- list(
    data_root = getwd(),
    catalog_csv = "wetlsp_cyverse_site_catalog_final.csv",
    mode = "local",
    default_site_id = NULL,
    default_year = 2021,
    max_pixels_default = 250,
    max_pixels_limit = 5000
  )

  if (file.exists(path)) {
    if (!requireNamespace("yaml", quietly = TRUE)) {
      stop("Package 'yaml' is required to read config.yml. Install it or remove config.yml.")
    }
    user_cfg <- yaml::read_yaml(path)
    defaults[names(user_cfg)] <- user_cfg
  }

  defaults
}
