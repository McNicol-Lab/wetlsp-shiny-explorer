library(shiny)
library(bslib)
library(dplyr)
library(ggplot2)
library(leaflet)
library(leaflet.extras)
library(htmlwidgets)
library(sf)
library(terra)
library(DT)

source("R/config.R")
source("R/paths.R")
source("R/catalog.R")
source("R/parquet_access.R")
source("R/netcdf_access.R")

cfg <- read_app_config()
catalog <- load_site_catalog(cfg)

ui <- page_sidebar(
  title = "WetLSP PlanetScope EVI Explorer",
  theme = bs_theme(version = 5, bootswatch = "flatly"),
  sidebar = sidebar(
    width = 360,
    selectizeInput(
      "site_id",
      "Site",
      choices = setNames(catalog$site_id, paste0(catalog$site_id, " — ", catalog$site_name)),
      selected = cfg$default_site_id %||% catalog$site_id[[1]]
    ),
    selectInput(
      "year",
      "Year",
      choices = c("2021", "2022", "2023", "2024"),
      selected = as.character(cfg$default_year %||% 2021)
    ),
    checkboxGroupInput(
      "series",
      "Series",
      choices = c("spline", "raw"),
      selected = "spline",
      inline = TRUE
    ),
    sliderInput(
      "max_pixels",
      "Maximum pixels to sample",
      min = 0,
      max = cfg$max_pixels_limit %||% 5000,
      value = cfg$max_pixels_default %||% 250,
      step = 25
    ),
    checkboxInput(
      "use_pixel_selection",
      "Use selected pixels for EVI plot",
      value = FALSE
    ),
    actionButton(
      "clear_pixel_selection",
      "Clear selected pixels",
      class = "btn-outline-secondary btn-sm"
    ),
    helpText(
      "Pixel selection: on the Pixel map, use the rectangle/polygon tool to select a group, or right-click / secondary-click individual pixels to toggle them."
    ),
    dateRangeInput(
      "date_range",
      "Date range",
      start = "2021-01-01",
      end = "2021-12-31"
    ),
    actionButton("load_data", "Load / refresh time series", class = "btn-primary"),
    hr(),
    helpText("For large sites, start with 100–500 pixels. Use summaries for public deployment.")
  ),
  navset_card_tab(
    nav_panel(
      "Overview",
      layout_columns(
        card(
          card_header("Site map"),
          leafletOutput("site_map", height = 420)
        ),
        card(
          card_header("Selected site"),
          DTOutput("site_info")
        ),
        col_widths = c(7, 5)
      )
    ),
    nav_panel(
      "EVI time series",
      card(
        card_header("Pixel-level EVI trajectories and daily site summary"),
        plotOutput("evi_plot", height = 520),
        verbatimTextOutput("query_status")
      )
    ),
    nav_panel(
      "Pixel map",
      card(
        card_header("Pixel selection map"),
        leafletOutput("pixel_map", height = 560)
      )
    ),
    nav_panel(
      "NetCDF phenometrics",
      layout_columns(
        card(
          card_header("NetCDF layers"),
          DTOutput("nc_layers")
        ),
        card(
          card_header("NetCDF path"),
          verbatimTextOutput("nc_path")
        ),
        col_widths = c(7, 5)
      )
    ),
    nav_panel(
      "Phenometric visualizer",
      layout_columns(
        card(
          card_header("Controls"),
          selectizeInput(
            "pheno_nc_files",
            "Select up to four annual NetCDF files",
            choices = NULL,
            multiple = TRUE,
            options = list(maxItems = 4, placeholder = "Choose one to four site-year NetCDF files")
          ),
          uiOutput("pheno_layer_ui"),
          numericInput(
            "pheno_max_cells",
            "Maximum raster cells per panel",
            value = 50000,
            min = 5000,
            max = 250000,
            step = 5000
          ),
          actionButton("load_pheno", "Load phenometric rasters", class = "btn-primary"),
          verbatimTextOutput("pheno_status")
        ),
        card(
          card_header("Phenometric rasters"),
          plotOutput("pheno_plot_grid", height = 720)
        ),
        col_widths = c(4, 8)
      )
    ),
    nav_panel(
      "Catalog",
      card(
        card_header("Full site catalog"),
        DTOutput("catalog_table")
      )
    )
  )
)


# Read the full pixel geometry table for the selected site.
# This is separate from read_evi_subset(), which intentionally samples pixels for speed.
read_all_pixel_geometry <- function(cfg, site_row) {
  requireNamespace("arrow", quietly = TRUE)

  geom_path <- resolve_geom_path(cfg, site_row)
  geom <- arrow::open_dataset(geom_path) |>
    dplyr::select(dplyr::any_of(c("pixel_id", "cell", "x", "y", "lon", "lat"))) |>
    dplyr::collect()

  geom_to_lonlat_if_possible(geom, cfg, site_row)
}


# App-level EVI reader that supports either random pixel sampling or explicit
# user-selected pixels from the Leaflet map.
read_evi_subset_for_pixels <- function(
    cfg,
    site_row,
    year,
    series = "spline",
    date_range = NULL,
    max_pixels = 250,
    selected_pixel_ids = NULL
) {
  requireNamespace("arrow", quietly = TRUE)

  ts_path <- resolve_timeseries_path(cfg, site_row)
  geom_path <- resolve_geom_path(cfg, site_row)

  ts_ds <- arrow::open_dataset(ts_path)
  geom_ds <- arrow::open_dataset(geom_path)

  if (!is.null(selected_pixel_ids) && length(selected_pixel_ids) > 0) {
    pixel_ids <- unique(as.integer(selected_pixel_ids))
    pixel_ids <- pixel_ids[!is.na(pixel_ids)]
  } else {
    pixel_ids <- ts_ds |>
      dplyr::filter(.data$year == !!as.integer(year), .data$series %in% !!series) |>
      dplyr::distinct(.data$pixel_id) |>
      head(max_pixels) |>
      dplyr::collect() |>
      dplyr::pull(.data$pixel_id)
  }

  if (length(pixel_ids) == 0) {
    return(list(
      timeseries = tibble::tibble(),
      geometry = tibble::tibble(),
      paths = list(timeseries = ts_path, geometry = geom_path)
    ))
  }

  q <- ts_ds |>
    dplyr::filter(
      .data$year == !!as.integer(year),
      .data$series %in% !!series,
      .data$pixel_id %in% !!pixel_ids
    )

  if (!is.null(date_range) && length(date_range) == 2 && all(!is.na(date_range))) {
    q <- q |>
      dplyr::filter(
        .data$date >= as.Date(!!date_range[[1]]),
        .data$date <= as.Date(!!date_range[[2]])
      )
  }

  ts <- q |>
    dplyr::select(.data$pixel_id, .data$series, .data$date, .data$year, .data$evi) |>
    dplyr::collect() |>
    dplyr::mutate(date = as.Date(.data$date))

  geom <- geom_ds |>
    dplyr::filter(.data$pixel_id %in% !!pixel_ids) |>
    dplyr::collect()

  geom <- geom_to_lonlat_if_possible(geom, cfg, site_row)

  list(
    timeseries = ts,
    geometry = geom,
    paths = list(timeseries = ts_path, geometry = geom_path)
  )
}

# A more reliable leaflet basemap helper: OSM is the default because it is often
# the least fragile; Carto and Esri are provided as alternate layers.
add_reliable_base_tiles <- function(map) {
  map |>
    addTiles(
      group = "OpenStreetMap",
      urlTemplate = "https://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png",
      attribution = '&copy; <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a> contributors',
      options = tileOptions(noWrap = TRUE)
    ) |>
    addProviderTiles(providers$CartoDB.Positron, group = "CartoDB Positron") |>
    addProviderTiles(providers$Esri.WorldImagery, group = "Esri World Imagery")
}



cap_selected_pixels <- function(pixel_ids, limit = 500) {
  pixel_ids <- unique(as.character(pixel_ids))
  limit <- suppressWarnings(as.integer(limit))
  if (is.na(limit) || limit < 1) limit <- 500L

  if (length(pixel_ids) > limit) {
    list(
      pixel_ids = pixel_ids[seq_len(limit)],
      n_original = length(pixel_ids),
      n_kept = limit,
      was_limited = TRUE
    )
  } else {
    list(
      pixel_ids = pixel_ids,
      n_original = length(pixel_ids),
      n_kept = length(pixel_ids),
      was_limited = FALSE
    )
  }
}



# Resolve all annual NetCDF files available for a selected site row.
available_netcdf_choices <- function(cfg, site_row) {
  years <- c(2021, 2022, 2023, 2024)
  out <- lapply(years, function(yr) {
    col <- paste0("netcdf_", yr, "_file")
    wet_col <- paste0("wetlsp_", yr)

    if (!col %in% names(site_row)) return(NULL)

    rel <- site_row[[col]][[1]]
    available_flag <- if (wet_col %in% names(site_row)) site_row[[wet_col]][[1]] else 1

    if (is.na(rel) || !nzchar(rel) || is.na(available_flag) || as.integer(available_flag) != 1L) {
      return(NULL)
    }

    path <- path_join_portable(cfg$data_root, rel)
    label <- paste0(site_row$site_id[[1]], " — ", yr)
    data.frame(label = label, year = yr, path = path, stringsAsFactors = FALSE)
  })

  dplyr::bind_rows(out)
}

safe_read_netcdf_rast <- function(path) {
  tryCatch({
    terra::rast(path)
  }, error = function(e) {
    warning("Could not read NetCDF: ", path, " | ", conditionMessage(e))
    NULL
  })
}

netcdf_layer_names <- function(paths) {
  paths <- paths[!is.na(paths) & nzchar(paths)]
  if (length(paths) == 0) return(character())

  layer_sets <- lapply(paths, function(p) {
    r <- safe_read_netcdf_rast(p)
    if (is.null(r)) character() else names(r)
  })

  layer_sets <- layer_sets[lengths(layer_sets) > 0]
  if (length(layer_sets) == 0) return(character())

  Reduce(intersect, layer_sets)
}

read_one_netcdf_layer_df <- function(path, layer_name, max_cells = 50000) {
  r <- safe_read_netcdf_rast(path)
  if (is.null(r)) return(NULL)

  if (!layer_name %in% names(r)) {
    warning("Layer not found in NetCDF: ", layer_name)
    return(NULL)
  }

  lyr <- r[[layer_name]]

  # Downsample for Shiny plotting responsiveness.
  n <- terra::ncell(lyr)
  if (!is.na(n) && n > max_cells) {
    fact <- ceiling(sqrt(n / max_cells))
    lyr <- terra::aggregate(lyr, fact = fact, fun = mean, na.rm = TRUE)
  }

  df <- terra::as.data.frame(lyr, xy = TRUE, na.rm = FALSE)
  names(df) <- c("x", "y", "value")
  df
}

phenometric_scale_type <- function(layer_name) {
  nm <- tolower(layer_name %||% "")

  # Date/timing phenometrics: day-of-year values, season timing, peak timing.
  date_patterns <- c(
    "sos", "eos", "pos", "pop", "greenup", "green_up", "green_down",
    "senescence", "dormancy", "maturity", "midgreen", "mid_green",
    "date", "doy", "day", "onset", "offset", "timing", "peak_date",
    "max_date", "min_date"
  )

  # Greenness / magnitude / integral phenometrics.
  green_patterns <- c(
    "evi", "max", "min", "amp", "amplitude", "base", "peak",
    "auc", "integral", "small_integral", "large_integral",
    "area", "value", "mean", "median"
  )

  if (any(vapply(date_patterns, grepl, logical(1), x = nm, fixed = TRUE))) {
    return("date")
  }

  if (any(vapply(green_patterns, grepl, logical(1), x = nm, fixed = TRUE))) {
    return("greenness")
  }

  "generic"
}

plot_netcdf_layer <- function(df, title = NULL, layer_name = NULL) {
  if (is.null(df) || nrow(df) == 0) {
    return(
      ggplot2::ggplot() +
        ggplot2::theme_void() +
        ggplot2::labs(title = title %||% "No data")
    )
  }

  scale_type <- phenometric_scale_type(layer_name %||% title)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$x, y = .data$y, fill = .data$value)) +
    ggplot2::geom_raster() +
    ggplot2::coord_equal(expand = FALSE) +
    ggplot2::labs(title = title, x = NULL, y = NULL, fill = NULL) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.key.width = grid::unit(1.3, "cm"),
      plot.title = ggplot2::element_text(size = 10)
    )

  if (identical(scale_type, "date")) {
    # Moon et al.-style high-contrast timing palette for DOY/date phenometrics.
    # Approximate visual: blue/cyan/green/yellow/orange/red/magenta contrast.
    p <- p +
      ggplot2::scale_fill_gradientn(
        colours = c(
          "#2c00ff", "#004cff", "#00b7ff", "#00e5a8",
          "#7dff00", "#ffff00", "#ff9e00", "#ff0000", "#9e0000"
        ),
        na.value = "transparent",
        guide = ggplot2::guide_colourbar(title.position = "top")
      )
  } else if (identical(scale_type, "greenness")) {
    # Moon et al.-style broad green scale for EVI magnitude/amplitude/integral metrics.
    p <- p +
      ggplot2::scale_fill_gradientn(
        colours = c(
          "#f7fcf5", "#e5f5e0", "#c7e9c0", "#a1d99b",
          "#74c476", "#41ab5d", "#238b45", "#006d2c", "#00441b"
        ),
        na.value = "transparent",
        guide = ggplot2::guide_colourbar(title.position = "top")
      )
  } else {
    p <- p +
      ggplot2::scale_fill_viridis_c(
        na.value = "transparent",
        guide = ggplot2::guide_colourbar(title.position = "top")
      )
  }

  p
}

server <- function(input, output, session) {

  selected_site <- reactive({
    catalog |> filter(.data$site_id == input$site_id) |> slice(1)
  })

  selected_pixel_ids <- reactiveVal(character())

  observeEvent(input$site_id, {
    selected_pixel_ids(character())
  })

  observeEvent(input$clear_pixel_selection, {
    selected_pixel_ids(character())
  })

  observeEvent(input$pixel_map_marker_contextmenu, {
    px <- input$pixel_map_marker_contextmenu$id
    if (!is.null(px) && nzchar(px)) {
      current <- as.character(selected_pixel_ids())
      px <- as.character(px)

      if (px %in% current) {
        selected_pixel_ids(setdiff(current, px))
      } else {
        capped <- cap_selected_pixels(unique(c(current, px)), input$selection_pixel_limit)
        if (isTRUE(capped$was_limited)) {
          showNotification(
            paste0("Selection limit reached: ", capped$n_kept, " pixels."),
            type = "warning",
            duration = 4
          )
        }
        selected_pixel_ids(capped$pixel_ids)
      }
    }
  })

  observeEvent(input$pixel_map_draw_new_feature, {
    feature <- input$pixel_map_draw_new_feature
    if (is.null(feature$geometry)) return()

    geom_all <- all_pixel_geom()
    if (is.null(geom_all) || nrow(geom_all) == 0) return()

    # Robust GeoJSON polygon/rectangle selection.
    # Leaflet draw sends coordinates as lon/lat GeoJSON. Using sf avoids the
    # earlier failure mode where only the horizontal/longitude domain was honored.
    inside <- tryCatch({
      coords <- feature$geometry$coordinates

      # Polygon and rectangle GeoJSON are usually: list(ring), where ring is
      # list(c(lon, lat), c(lon, lat), ...). MultiPolygon adds one more level.
      if (identical(feature$geometry$type, "MultiPolygon")) {
        ring <- coords[[1]][[1]]
      } else {
        ring <- coords[[1]]
      }

      xy <- do.call(rbind, lapply(ring, function(z) {
        c(as.numeric(z[[1]]), as.numeric(z[[2]]))
      }))

      xy <- xy[stats::complete.cases(xy), , drop = FALSE]
      if (nrow(xy) < 3) stop("Drawn geometry has fewer than three valid vertices.")

      # Ensure closed polygon ring.
      if (!all(xy[1, ] == xy[nrow(xy), ])) {
        xy <- rbind(xy, xy[1, ])
      }

      poly <- sf::st_sfc(sf::st_polygon(list(xy)), crs = 4326)

      pts <- geom_all |>
        dplyr::filter(!is.na(.data$lon), !is.na(.data$lat)) |>
        sf::st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE)

      hit <- as.logical(sf::st_intersects(pts, poly, sparse = FALSE)[, 1])

      pts |>
        dplyr::filter(hit) |>
        sf::st_drop_geometry() |>
        dplyr::pull(.data$pixel_id) |>
        as.character()
    }, error = function(e) {
      warning("Pixel draw selection failed: ", conditionMessage(e))
      character()
    })

    new_selection <- unique(c(as.character(selected_pixel_ids()), inside))
    capped <- cap_selected_pixels(new_selection, input$selection_pixel_limit)

    if (isTRUE(capped$was_limited)) {
      showNotification(
        paste0(
          "Polygon selected ", capped$n_original, " pixels; keeping first ",
          capped$n_kept, " to avoid loading an excessive time series."
        ),
        type = "warning",
        duration = 8
      )
    } else {
      showNotification(
        paste0("Selected ", capped$n_kept, " pixels."),
        type = "message",
        duration = 4
      )
    }

    selected_pixel_ids(capped$pixel_ids)
  })


  observeEvent(input$year, {
    updateDateRangeInput(
      session,
      "date_range",
      start = paste0(input$year, "-01-01"),
      end = paste0(input$year, "-12-31")
    )
  }, ignoreInit = FALSE)

  output$site_map <- renderLeaflet({
    leaflet(catalog, options = leafletOptions(preferCanvas = TRUE)) |>
      add_reliable_base_tiles() |>
      addCircleMarkers(
        lng = ~lon, lat = ~lat,
        radius = 5,
        stroke = TRUE,
        weight = 1,
        color = "#1f78b4",
        fillColor = "#1f78b4",
        fillOpacity = 0.8,
        label = ~paste0(site_id, " — ", site_name),
        layerId = ~site_id,
        group = "WetLSP sites"
      ) |>
      addLayersControl(
        baseGroups = c("OpenStreetMap", "CartoDB Positron", "Esri World Imagery"),
        overlayGroups = c("WetLSP sites"),
        options = layersControlOptions(collapsed = TRUE)
      ) |>
      fitBounds(
        lng1 = min(catalog$lon, na.rm = TRUE),
        lat1 = min(catalog$lat, na.rm = TRUE),
        lng2 = max(catalog$lon, na.rm = TRUE),
        lat2 = max(catalog$lat, na.rm = TRUE)
      )
  })

  observeEvent(input$site_map_marker_click, {
    click <- input$site_map_marker_click
    if (!is.null(click$id)) updateSelectizeInput(session, "site_id", selected = click$id)
  })

  output$site_info <- renderDT({
    selected_site() |>
      select(any_of(c(
        "site_id", "site_name", "country", "lat", "lon", "base_network",
        "final_fetch_r", "tower_height_m", "canopy_height_m",
        "base_2021", "base_2022", "base_2023", "base_2024",
        "wetlsp_2021", "wetlsp_2022", "wetlsp_2023", "wetlsp_2024"
      ))) |>
      tidyr::pivot_longer(everything(), names_to = "field", values_to = "value") |>
      datatable(options = list(dom = "t", pageLength = 30), rownames = FALSE)
  })

  all_pixel_geom <- reactive({
    req(input$site_id)
    site <- selected_site()

    tryCatch(
      read_all_pixel_geometry(cfg, site),
      error = function(e) tibble::tibble()
    )
  })

  evi_data <- eventReactive(input$load_data, {
    req(input$site_id, input$year, input$series)
    site <- selected_site()
    validate_site_year_available(site, input$year)

    selected_ids <- if (isTRUE(input$use_pixel_selection) && length(selected_pixel_ids()) > 0) {
      capped <- cap_selected_pixels(selected_pixel_ids(), input$selection_pixel_limit)
      if (isTRUE(capped$was_limited)) {
        showNotification(
          paste0("Loading only ", capped$n_kept, " selected pixels out of ", capped$n_original, "."),
          type = "warning",
          duration = 6
        )
      }
      capped$pixel_ids
    } else {
      NULL
    }

    read_evi_subset_for_pixels(
      cfg = cfg,
      site_row = site,
      year = as.integer(input$year),
      series = input$series,
      date_range = input$date_range,
      max_pixels = input$max_pixels,
      selected_pixel_ids = selected_ids
    )
  }, ignoreInit = TRUE)

  output$query_status <- renderPrint({
    dat <- evi_data()
    cat("Rows loaded:", nrow(dat$timeseries), "\n")
    pixels_loaded <- if ("pixel_id" %in% names(dat$timeseries)) dplyr::n_distinct(dat$timeseries$pixel_id) else 0
    cat("Pixels loaded:", pixels_loaded, "\n")
    cat("Pixels currently selected:", length(selected_pixel_ids()), "\n")
    cat("Selection pixel load limit:", input$selection_pixel_limit, "\n")
    cat("Using selected pixels:", isTRUE(input$use_pixel_selection) && length(selected_pixel_ids()) > 0, "\n")
    cat("Time-series path:", dat$paths$timeseries, "\n")
  })

  output$evi_plot <- renderPlot({
    dat <- evi_data()
    ts <- dat$timeseries

    validate(need(nrow(ts) > 0, "No rows returned for this filter."))

    daily <- ts |>
      group_by(date, series) |>
      summarize(
        evi_mean = mean(evi, na.rm = TRUE),
        evi_q25 = quantile(evi, 0.25, na.rm = TRUE),
        evi_q75 = quantile(evi, 0.75, na.rm = TRUE),
        .groups = "drop"
      )

    ggplot() +
      geom_line(
        data = ts,
        aes(x = date, y = evi, group = interaction(pixel_id, series)),
        alpha = 0.08
      ) +
      geom_ribbon(
        data = daily,
        aes(x = date, ymin = evi_q25, ymax = evi_q75, fill = series),
        alpha = 0.18
      ) +
      geom_line(
        data = daily,
        aes(x = date, y = evi_mean, color = series),
        linewidth = 1.0
      ) +
      labs(
        x = NULL,
        y = "EVI",
        title = paste(input$site_id, input$year, "EVI time series"),
        subtitle = if (isTRUE(input$use_pixel_selection) && length(selected_pixel_ids()) > 0) "Thin lines are selected pixels; bold line is daily mean; ribbon is interquartile range." else "Thin lines are sampled pixels; bold line is daily mean; ribbon is interquartile range."
      ) +
      theme_minimal(base_size = 13)
  })

  output$pixel_map <- renderLeaflet({
    all_geom <- all_pixel_geom()

    validate(need(nrow(all_geom) > 0, "No full-site geometry rows returned."))

    leaflet(options = leafletOptions(preferCanvas = TRUE)) |>
      add_reliable_base_tiles() |>
      addCircleMarkers(
        data = all_geom,
        lng = ~lon, lat = ~lat,
        radius = 2,
        stroke = FALSE,
        fillColor = "#bdbdbd",
        fillOpacity = 0.30,
        label = ~paste0("pixel_id: ", pixel_id),
        layerId = ~as.character(pixel_id),
        group = "All pixels"
      ) |>
      addDrawToolbar(
        targetGroup = "Selection boxes",
        polylineOptions = FALSE,
        markerOptions = FALSE,
        circleMarkerOptions = FALSE,
        circleOptions = FALSE,
        polygonOptions = drawPolygonOptions(showArea = TRUE),
        rectangleOptions = drawRectangleOptions(showArea = TRUE),
        editOptions = editToolbarOptions(edit = FALSE, remove = TRUE)
      ) |>
      addLayersControl(
        baseGroups = c("OpenStreetMap", "CartoDB Positron", "Esri World Imagery"),
        overlayGroups = c("All pixels", "Current EVI pixels", "Selected pixels", "Selection boxes"),
        options = layersControlOptions(collapsed = TRUE)
      ) |>
      fitBounds(
        lng1 = min(all_geom$lon, na.rm = TRUE),
        lat1 = min(all_geom$lat, na.rm = TRUE),
        lng2 = max(all_geom$lon, na.rm = TRUE),
        lat2 = max(all_geom$lat, na.rm = TRUE)
      ) |>
      htmlwidgets::onRender("
        function(el, x) {
          var map = this;
          map.eachLayer(function(layer) {
            if (layer.options && layer.options.layerId) {
              layer.off('contextmenu');
              layer.on('contextmenu', function(e) {
                if (e.originalEvent) {
                  e.originalEvent.preventDefault();
                }
                Shiny.setInputValue(
                  el.id + '_marker_contextmenu',
                  {id: layer.options.layerId, nonce: Math.random()},
                  {priority: 'event'}
                );
              });
            }
          });
        }
      ")
  })

  # Update selected-pixel overlay without rebuilding the whole Leaflet map.
  # This is much more stable for large pixel clouds than making selected_pixel_ids()
  # a dependency of renderLeaflet().
  observeEvent(selected_pixel_ids(), {
    all_geom <- all_pixel_geom()
    ids <- as.character(selected_pixel_ids())

    selected_geom <- all_geom |>
      dplyr::filter(as.character(.data$pixel_id) %in% ids)

    proxy <- leafletProxy("pixel_map") |>
      clearGroup("Selected pixels")

    if (nrow(selected_geom) > 0) {
      proxy |>
        addCircleMarkers(
          data = selected_geom,
          lng = ~lon, lat = ~lat,
          radius = 5,
          stroke = TRUE,
          weight = 2,
          color = "#08306b",
          fillColor = "#2171b5",
          fillOpacity = 1,
          label = ~paste0("selected pixel_id: ", pixel_id),
          layerId = ~as.character(pixel_id),
          group = "Selected pixels"
        )
    }
  }, ignoreInit = TRUE)

  # Update current EVI/query-pixel overlay only after the user explicitly loads data.
  observeEvent(evi_data(), {
    dat <- evi_data()
    sampled_geom <- dat$geometry

    proxy <- leafletProxy("pixel_map") |>
      clearGroup("Current EVI pixels")

    if (!is.null(sampled_geom) && nrow(sampled_geom) > 0 && "pixel_id" %in% names(sampled_geom)) {
      proxy |>
        addCircleMarkers(
          data = sampled_geom,
          lng = ~lon, lat = ~lat,
          radius = 4,
          stroke = TRUE,
          weight = 1,
          color = "#1f78b4",
          fillColor = "#1f78b4",
          fillOpacity = 0.85,
          label = ~paste0("current EVI pixel_id: ", pixel_id),
          layerId = ~as.character(pixel_id),
          group = "Current EVI pixels"
        )
    }
  }, ignoreInit = TRUE)

  output$nc_path <- renderPrint({
    site <- selected_site()
    cat(resolve_netcdf_path(cfg, site, input$year), "\n")
  })

  output$nc_layers <- renderDT({
    site <- selected_site()
    path <- resolve_netcdf_path(cfg, site, input$year)
    layers <- safe_netcdf_layers(path)

    datatable(layers, options = list(pageLength = 24), rownames = FALSE)
  })


  observeEvent(selected_site(), {
    choices_df <- available_netcdf_choices(cfg, selected_site())

    choices <- if (nrow(choices_df) > 0) {
      stats::setNames(choices_df$path, choices_df$label)
    } else {
      character()
    }

    updateSelectizeInput(
      session,
      "pheno_nc_files",
      choices = choices,
      selected = head(unname(choices), 1),
      server = TRUE
    )
  }, ignoreInit = FALSE)

  output$pheno_layer_ui <- renderUI({
    req(input$pheno_nc_files)
    layers <- netcdf_layer_names(input$pheno_nc_files)

    if (length(layers) == 0) {
      return(helpText("No common NetCDF layers found across the selected file(s), or files could not be read."))
    }

    selectizeInput(
      "pheno_layer",
      "Phenometric layer",
      choices = layers,
      selected = layers[[1]],
      multiple = FALSE
    )
  })

  pheno_rasters <- eventReactive(input$load_pheno, {
    req(input$pheno_nc_files, input$pheno_layer)

    paths <- input$pheno_nc_files
    if (length(paths) > 4) paths <- paths[seq_len(4)]

    max_cells <- suppressWarnings(as.integer(input$pheno_max_cells))
    if (is.na(max_cells) || max_cells < 1000) max_cells <- 50000L

    lapply(paths, function(p) {
      df <- read_one_netcdf_layer_df(p, input$pheno_layer, max_cells = max_cells)
      list(path = p, data = df)
    })
  }, ignoreInit = TRUE)

  output$pheno_status <- renderPrint({
    req(input$pheno_nc_files)
    cat("Selected NetCDF files:", length(input$pheno_nc_files), "\n")
    cat("Layer:", input$pheno_layer %||% NA_character_, "\n")
    cat("Max cells per panel:", input$pheno_max_cells, "\n")
    cat("Paths:\n")
    cat(paste0(" - ", input$pheno_nc_files, collapse = "\n"), "\n")
  })

  output$pheno_plot_grid <- renderPlot({
    rasters <- pheno_rasters()
    validate(need(length(rasters) > 0, "Select one to four NetCDF files and click Load phenometric rasters."))

    plots <- lapply(rasters, function(x) {
      nm <- basename(x$path)
      plot_netcdf_layer(x$data, title = paste0(nm, " | ", input$pheno_layer), layer_name = input$pheno_layer)
    })

    # Use patchwork if available; otherwise fall back to gridExtra.
    if (requireNamespace("patchwork", quietly = TRUE)) {
      Reduce(`+`, plots) + patchwork::plot_layout(ncol = min(2, length(plots)))
    } else if (requireNamespace("gridExtra", quietly = TRUE)) {
      do.call(gridExtra::grid.arrange, c(plots, ncol = min(2, length(plots))))
    } else {
      showNotification("Install 'patchwork' or 'gridExtra' for multi-panel plotting.", type = "warning")
      print(plots[[1]])
    }
  })


  output$catalog_table <- renderDT({
    datatable(
      catalog,
      filter = "top",
      options = list(pageLength = 20, scrollX = TRUE),
      rownames = FALSE
    )
  })
}

shinyApp(ui, server)
