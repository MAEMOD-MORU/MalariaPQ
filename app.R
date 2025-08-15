library(shiny)
library(shinydashboard)
library(deSolve)

# Read model
source("model_D0_33_D1_67_baseline_change_birth_death_2betaSets.r")
source("init_parameter_equilibrium_birth_2_rate_2beta_res_fitted.R")
uganda_Incidence <- read.csv("data/uganda_Incidence.csv")

events <- list(
  func = function(time, state, parameters) {
    if (time == 12 * 14) {
      state["Ir0"] <- state["Ir0"] + state["Is0"] * parameters$prob_res
      state["Ir1"] <- state["Ir1"] + state["Is1"] * parameters$prob_res
      state["Ir2"] <- state["Ir2"] + state["Is2"] * parameters$prob_res
      state["Ar"]  <- state["Ar"]  + state["As"]  * parameters$prob_res
      
      state["Is0"] <- state["Is0"] *(1-parameters$prob_res)
      state["Is1"] <- state["Is1"] *(1-parameters$prob_res)
      state["Is2"] <- state["Is2"] *(1-parameters$prob_res)
      state["As"]  <- state["As"]  *(1-parameters$prob_res)
    }
    return(state)
  },
  time = 12 * 14
)

init_state <- readRDS("init_state.rds")

summarise_by_year <- function(out, vars_all, n_years, start_year = 2000) {
  out_in_year <- data.frame(matrix(NA, nrow = n_years, ncol = length(vars_all) + 1))
  colnames(out_in_year) <- c("year", vars_all)
  for (i in seq_len(n_years)) {
    start_idx <- 12 * (i - 1) + 1
    end_idx <- 12 * i
    out_in_year[i, "year"] <- start_year + (i - 1)
    for (v in vars_all) {
      if (v %in% colnames(out)) {
        if (v == "N") {
          out_in_year[i, v] <- out[end_idx, v]
        } else {
          out_in_year[i, v] <- sum(out[start_idx:end_idx, v], na.rm = TRUE)
        }
      } else {
        out_in_year[i, v] <- NA
      }
    }
  }
  return(out_in_year)
}

ui <- dashboardPage(
  dashboardHeader(title = "Malaria"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Incidence", tabName = "incidence", icon = icon("chart-line")),
      menuItem("Gametocyte", tabName = "gametocyte", icon = icon("bug")),
      br(),
      sliderInput("start_res","Number of res (%):",min = 0,max = 50,step = 1,value = parameters$prob_res*100),
      sliderInput("prop_sym","proportion of symptomatic :",min = 0.01,max = 1,step = 0.01,value = (25/75)),
      h4("Proportion of D0,D1,D2"),
      sliderInput("D0", "proportion of D0", min = 0, max = 1, value = 0.34, step = 0.01),
      sliderInput("D1", "proportion of D1", min = 0, max = 1, value = 0.33, step = 0.01),
      sliderInput("D2", "proportion of D2", min = 0, max = 1, value = 0.33, step = 0.01),
      sliderInput("D2_start","start of D2 (year) :",min = 2023,max = 2033,step = 1,value = 2025)
    )
  ),
  dashboardBody(
    tabItems(
      tabItem(tabName = "incidence",
              fluidRow(box(width = 12, plotOutput("distPlot"))),
              fluidRow(
                box(width = 6, plotOutput("Plot_sym_asym")),
                box(width = 6, plotOutput("Plot_s_r"))
              )
      ),
      tabItem(tabName = "gametocyte",
              fluidRow(box(width = 12, plotOutput("Plot_Gasym_sym"))),
              fluidRow(box(width = 12, plotOutput("Plot_GS_GR")))
      )
    )
  )
)

server <- function(input, output, session) {
  
  values <- reactiveValues(
    out_in_year = NULL,
    init_state = init_state,
    out = NULL,
    parameters = parameters
  )
  
  last_mod <- reactiveVal(list(id = "D0", value = 0))
  observeEvent(input$D0, { last_mod(list(id = "D0", value = input$D0)) })
  observeEvent(input$D1, { last_mod(list(id = "D1", value = input$D1)) })
  observeEvent(input$D2, { last_mod(list(id = "D2", value = input$D2)) })
  
  observeEvent(c(input$D0, input$D1, input$D2), {
    d0 <- isolate(input$D0); d1 <- isolate(input$D1); d2 <- isolate(input$D2)
    total <- d0 + d1 + d2
    last_mod_id <- last_mod()$id
    
    if (total > 1) {
      overvalue <- total - 1
      if (last_mod_id == "D0") {
        if (d1 >= overvalue) {
          updateSliderInput(session, "D1", value = d1 - overvalue)
        } else {
          updateSliderInput(session, "D1", value = 0)
          updateSliderInput(session, "D2", value = d2 - (overvalue - d1))
        }
      } else if (last_mod_id == "D1") {
        if (d0 >= overvalue) {
          updateSliderInput(session, "D0", value = d0 - overvalue)
        } else {
          updateSliderInput(session, "D0", value = 0)
          updateSliderInput(session, "D2", value = d2 - (overvalue - d0))
        }
      } else if (last_mod_id == "D2") {
        if (d0 >= overvalue) {
          updateSliderInput(session, "D0", value = d0 - overvalue)
        } else {
          updateSliderInput(session, "D0", value = 0)
          updateSliderInput(session, "D1", value = d1 - (overvalue - d0))
        }
      }
    } else if (total < 1) {
      remainder <- 1 - total
      if (last_mod_id == "D0") {
        if (d1 + remainder <= 1) {
          updateSliderInput(session, "D1", value = d1 + remainder)
        } else {
          remainder_for_d2 <- remainder - (1 - d1)
          updateSliderInput(session, "D1", value = 1)
          updateSliderInput(session, "D2", value = d2 + remainder_for_d2)
        }
      } else if (last_mod_id == "D1") {
        if (d0 + remainder <= 1) {
          updateSliderInput(session, "D0", value = d0 + remainder)
        } else {
          remainder_for_d2 <- remainder - (1 - d0)
          updateSliderInput(session, "D0", value = 1)
          updateSliderInput(session, "D2", value = d2 + remainder_for_d2)
        }
      } else if (last_mod_id == "D2") {
        if (d0 + remainder <= 1) {
          updateSliderInput(session, "D0", value = d0 + remainder)
        } else {
          remainder_for_d1 <- remainder - (1 - d0)
          updateSliderInput(session, "D0", value = 1)
          updateSliderInput(session, "D1", value = d1 + remainder_for_d1)
        }
      }
    }
    
    values$parameters$propIs <- c(input$D0, input$D1, input$D2)
    values$parameters$propIr <- c(input$D0, input$D1, input$D2)
  })
  
  observeEvent(input$D2_start, {
    year_d2 <- input$D2_start - 1999
    values$parameters$start_d <- 12 * year_d2
  })
  
  observeEvent(input$start_res, {
    values$parameters$prob_res <- input$start_res / 100
    events <- list(
      func = function(time, state, parameters) {
        if (time == 12 * 14) {
          state["Ir0"] <- state["Ir0"] + state["Is0"] * input$start_res
          state["Ir1"] <- state["Ir1"] + state["Is1"] * input$start_res
          state["Ir2"] <- state["Ir2"] + state["Is2"] * input$start_res
          state["Ar"]  <- state["Ar"]  + state["As"]  * input$start_res
          
          state["Is0"] <- state["Is0"] *(1-input$start_res)
          state["Is1"] <- state["Is1"] *(1-input$start_res)
          state["Is2"] <- state["Is2"] *(1-input$start_res)
          state["As"]  <- state["As"]  *(1-input$start_res)
        }
        return(state)
      },
      time = 12 * 14
    )
  })
  
  observeEvent(input$prop_sym, {
    values$parameters$prob_sym_s <- input$prop_sym
    values$parameters$prob_sym_r <- input$prop_sym
  })
  
  ode_results <- reactive({
    out <- ode(
      y = values$init_state,
      times = times,
      func = Malaria_model_with_Array,
      parms = values$parameters,
      events = events
    )
    vars_all <- colnames(out)
    summarise_by_year(out, vars_all, 36)
  })
  
  ode_baseline <- reactive({
    parameters_baseline <- values$parameters
    parameters_baseline$start_d <- 6000
    out <- ode(
      y = values$init_state,
      times = times,
      func = Malaria_model_with_Array,
      parms = parameters_baseline,
      events = events
    )
    vars_all <- colnames(out)
    summarise_by_year(out, vars_all, 36)
  })
  
  output$distPlot <- renderPlot({
    out_in_year <- ode_results()
    out_in_year_baseline <- ode_baseline()
    plot(out_in_year[, 1], out_in_year[, "inc"], type = "l",
         xlab = "Years", ylab = "Incidence",
         main = "Total Incidence by Year",
         xlim = c(2023, 2035), lwd = 2, lty = 2,
         col = 1, ylim = c(0, 1.2 * max(out_in_year[1:35, "inc"])))
    points(2000:2022, uganda_Incidence[1:23, 2], col = "red", pch = 19)
    lines(out_in_year_baseline[, 1], out_in_year_baseline[, "inc"], lwd = 2)
    abline(v = input$D2_start, lty = 2)
    text(input$D2_start - 1, 2e07, "Start D2")
  })
  
  output$Plot_sym_asym <- renderPlot({
    out_in_year <- ode_results()
    out_in_year_baseline <- ode_baseline()
    plot(out_in_year[, 1], out_in_year[, "inc_sym"], type = "l",
         xlab = "Years", ylab = "Incidence",
         main = "Symptomatic/Asymptomatic Incidence by Year",
         xlim = c(2023, 2035), lwd = 2, lty = 2,
         col = "blue2", ylim = c(0, 1.2 * max(c(out_in_year[1:35, "inc_sym"], out_in_year[1:35, "inc_asym"]))))
    lines(out_in_year[, 1], out_in_year[, "inc_asym"], col = "green4", lty = 2, lwd = 2)
    lines(out_in_year_baseline[, 1], out_in_year_baseline[, "inc_sym"], col = "blue2", lwd = 2)
    lines(out_in_year_baseline[, 1], out_in_year_baseline[, "inc_asym"], col = "green4", lwd = 2)
    legend("bottomright", c("Symptomatic", "Asymptomatic"), col = c("blue2", "green4"), lty = c(1, 1), lwd = 2)
    abline(v = input$D2_start, lty = 2)
    text(input$D2_start - 1, 2e07, "Start D2")
  })
  
  output$Plot_s_r <- renderPlot({
    out_in_year <- ode_results()
    out_in_year_baseline <- ode_baseline()
    plot(out_in_year[, 1], out_in_year[, "inc_s"], type = "l",
         xlab = "Years", ylab = "Incidence",
         main = "Sensitive/Resistant Incidence by Year",
         xlim = c(2023, 2035), lwd = 2, lty = 2,
         col = "black", ylim = c(0, 1.2 * max(c(out_in_year[1:35, "inc_s"], out_in_year[1:35, "inc_r"]))))
    lines(out_in_year[, 1], out_in_year[, "inc_r"], col = "red", lty = 2, lwd = 2)
    lines(out_in_year_baseline[, 1], out_in_year_baseline[, "inc_r"], col = "red", lwd = 2)
    lines(out_in_year_baseline[, 1], out_in_year_baseline[, "inc_s"], lwd = 2)
    legend("bottomright", c("Sensitive", "Resistant"), col = c("black", "red"), lty = c(1, 1), lwd = 2)
    abline(v = input$D2_start, lty = 2)
    text(input$D2_start - 1, 2e07, "Start D2")
  })
  
  output$Plot_Gasym_sym <- renderPlot({
    out_in_year <- ode_results()
    out_in_year_baseline <- ode_baseline()
    plot(out_in_year[, 1], out_in_year[, "Gasym_inc"], type = "l",
         xlab = "Years", ylab = "Incidence",lty = 2,
         main = "Gametocyte Infections Asymptomatic/Symptomatic by Year",
         xlim = c(2023, 2035), lwd = 2, col = 3,
         ylim = c(0, 1.2 * max(out_in_year[1:35, "Gasym_inc"])))
    lines(out_in_year_baseline[, 1], out_in_year_baseline[, "Gasym_inc"], lty = 1, lwd = 2, col = 3)
    lines(out_in_year[, 1], out_in_year[, "Gsym_inc"], lwd = 2,lty = 2, col = "blue")
    lines(out_in_year_baseline[, 1], out_in_year_baseline[, "Gsym_inc"], lwd = 2, col = "blue")
    legend("topleft", c("Asymptomatic", "Symptomatic"), col = c(3, 4), lty = c(1, 1))
    abline(v = input$D2_start, lty = 2)
    text(input$D2_start - 1, 4e05, "Start D2")
  })
  
  output$Plot_GS_GR <- renderPlot({
    out_in_year <- ode_results()
    out_in_year_baseline <- ode_baseline()
    plot(out_in_year[, 1], out_in_year[, "GS_inc"], type = "l",
         xlab = "Years", ylab = "Incidence",
         main = "Gametocyte Infections Sensitive/Resistant by Year",
         xlim = c(2023, 2035), lwd = 2, col = 1,lty = 2,
         ylim = c(0, 1.2 * max(out_in_year[1:35, "GR_inc"])))
    lines(out_in_year_baseline[, 1], out_in_year_baseline[, "GS_inc"], lty = 1, lwd = 2)
    lines(out_in_year[, 1], out_in_year[, "GR_inc"], lwd = 2,lty = 2, col = "red")
    lines(out_in_year_baseline[, 1], out_in_year_baseline[, "GR_inc"], lwd = 2, col = "red")
    legend("topleft", c("Sensitive", "Resistant"), col = c(1, 2), lty = c(1, 1))
    abline(v = input$D2_start, lty = 2)
    text(input$D2_start - 1, 4e05, "Start D2")
  })
}

shinyApp(ui = ui, server = server)
