library(shiny)
library(bslib)
library(dplyr)
library(ggplot2)
library(leaflet)
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
      min = 25,
      max = cfg$max_pixels_limit %||% 5000,
      value = cfg$max_pixels_default %||% 250,
      step = 25
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
        card_header("Sampled pixel centroids"),
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

server <- function(input, output, session) {

  selected_site <- reactive({
    catalog |> filter(.data$site_id == input$site_id) |> slice(1)
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
    leaflet(catalog) |>
      addProviderTiles(providers$CartoDB.Positron) |>
      addCircleMarkers(
        lng = ~lon, lat = ~lat,
        radius = 5,
        stroke = FALSE,
        fillOpacity = 0.75,
        label = ~paste0(site_id, " — ", site_name),
        layerId = ~site_id
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

  evi_data <- eventReactive(input$load_data, {
    req(input$site_id, input$year, input$series)
    site <- selected_site()
    validate_site_year_available(site, input$year)

    read_evi_subset(
      cfg = cfg,
      site_row = site,
      year = as.integer(input$year),
      series = input$series,
      date_range = input$date_range,
      max_pixels = input$max_pixels
    )
  }, ignoreInit = TRUE)

  output$query_status <- renderPrint({
    dat <- evi_data()
    cat("Rows loaded:", nrow(dat$timeseries), "\n")
    cat("Pixels loaded:", dplyr::n_distinct(dat$timeseries$pixel_id), "\n")
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
        subtitle = "Thin lines are sampled pixels; bold line is daily mean; ribbon is interquartile range."
      ) +
      theme_minimal(base_size = 13)
  })

  output$pixel_map <- renderLeaflet({
    dat <- evi_data()
    geom <- dat$geometry

    validate(need(nrow(geom) > 0, "No geometry rows returned."))

    leaflet(geom) |>
      addProviderTiles(providers$CartoDB.Positron) |>
      addCircleMarkers(
        lng = ~lon, lat = ~lat,
        radius = 3,
        stroke = FALSE,
        fillOpacity = 0.6,
        label = ~paste0("pixel_id: ", pixel_id)
      )
  })

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
