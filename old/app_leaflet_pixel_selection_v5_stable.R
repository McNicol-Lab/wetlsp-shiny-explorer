library(shiny)
library(bslib)
library(dplyr)
library(ggplot2)
library(leaflet)
library(leaflet.extras)
library(htmlwidgets)
library(sf)
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
    pixel_ids <- as.integer(selected_pixel_ids)
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
        selected_pixel_ids(unique(c(current, px)))
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

    selected_pixel_ids(unique(c(as.character(selected_pixel_ids()), inside)))
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
      selected_pixel_ids()
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
