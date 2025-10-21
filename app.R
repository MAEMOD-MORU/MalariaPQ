library(shiny)
library(shinydashboard)
library(htmltools)
library(shinyBS)
library(deSolve)
library(shinyjs)
library(bsplus)
library(shinycssloaders)

# Read model
source("model_D0_33_D1_67_baseline_change_birth_death_2betaSets.R")
source("init_parameter_equilibrium_birth_2_rate_2beta_res_fitted.R")
uganda_Incidence <- read.csv("data/uganda_Incidence.csv")
init_state <- readRDS("init_state.rds")

source("init_parameter_equilibrium_birth_2_rate_2beta_res_fitted_Tanzania.R")
tanzania_Incidence <- read.csv("data/Tanzania_incidence.csv")
init_state_tanzani <- readRDS("init_state_Tanzania.rds")

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
        }
        else if(v == "p_resistant_inc" | v == "F_T" | v == "CT_total" | v == "CU_total"){
          out_in_year[i, v] <- mean(out[start_idx:end_idx, v], na.rm = TRUE)
        }
        else {
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

  dashboardHeader(title = "MalariaPQ",
                  tags$li(actionLink("goto_intro", HTML("<b>Introduction</b>"), class = "btn btn-default"), class = "dropdown"),
                  tags$li(actionLink("goto_simulation", HTML("<b>Simulation</b>"), class = "btn btn-default"), class = "dropdown"),
                  tags$li(actionLink("goto_about", HTML("<b>About us</b>"), class = "btn btn-default"), class = "dropdown"),
                  tags$li(actionLink("goto_data", HTML("<b>Source Data</b>"), class = "btn btn-default"), class = "dropdown")
                  ),
  dashboardSidebar(
    sidebarMenu(

      # useShinyjs(),
      actionButton("reset_all", "Reset", icon = icon("undo")),
      radioButtons("country" , "Data Country",choices =c("Uganda","Tanzania") ),
      
      # Transmissibility group
      menuItem("Transmissibility", tabName = "transmissibility", icon = icon("bug"),
               bsPopover(
                 id = "Sensitive_2000_2014",
                 title = NULL,
                 content = "β is defined separately for sensitive/resistant strains. Sensitive infections are split into pre- and post-.",
                 placement = "right",
                 trigger = "hover"
               ),
               menuItem(tagList(tags$span("Sensitive", tags$i(class = "fa fa-circle-info"), id = "Sensitive_2000_2014")), 
                        tabName = "transmissibility_sensitive",
                        bsPopover(
                          id = "beta1",
                          title = NULL,
                          content = "Uganda: β1 (2000–2014); Tanzania: β1 (2015–2020)",
                          placement = "right",
                          trigger = "hover"
                        ),
                        menuItem(tagList(tags$span("Beta1", tags$i(class = "fa fa-circle-info"), id = "beta1")), 
                                 tabName = "Beta_1",
                                 sliderInput("Bss", "β (Symptomatic)", 
                                             min = 0, max = 10, step = 0.01, value = parameters$beta_s[1]),
                                 sliderInput("Bas", "β (Asymptomatic)", 
                                             min = 0, max = 10, step = 0.01, value = parameters$beta_as)
                        ),
                        bsPopover(
                          id = "beta2",
                          title = NULL,
                          content = "Uganda: β2 (2014–2035); Tanzania: β2 (2020–2035)",
                          placement = "right",
                          trigger = "hover"
                        ),
                        menuItem(tagList(tags$span("Beta2", tags$i(class = "fa fa-circle-info"), id = "beta2")), 
                                 tabName = "2014–2035",
                                 sliderInput("Bss2", "β (Symptomatic)", 
                                             min = 0, max = 10, step = 0.01, value = parameters$beta_s_2[1]),
                                 sliderInput("Bas2", "β (Asymptomatic)", 
                                             min = 0, max = 10, step = 0.01, value = parameters$beta_as_2)
                        )
               ),
               menuItem("Resistant", tabName = "transmissibility_resistant",
                           sliderInput("Bsr", "β (Symptomatic)", 
                                       min = 0, max = 10, step = 0.01, value = parameters$beta_r[1]),
                           sliderInput("Bar", "β (Asymptomatic)", 
                                       min = 0, max = 10, step = 0.01, value = parameters$beta_ar)
               ),
               menuItem("Treatment-associated",tabName = "Treatment_associated",
                        sliderInput("s_hat","s_hat", 
                                    min = 0.001, max = 4, step = 0.001, value = 0.35)
               )
      ),
      
      # Symptomatic group
      menuItem("Symptomatic", tabName = "symptomatic", icon = icon("user-md"),
               sliderInput("prop_sym","Proportion Symptomatic", 
                           min = 0.01, max = 1, step = 0.01, value = 0.25)
      ),
      
      # Resistance group
      menuItem("Resistance", tabName = "resistance", icon = icon("dna"),
               sliderInput("start_res_year","Start Year of Resistance", 
                           min = 2010, max = 2020, step = 1, value = 2014,sep = ""),
               sliderInput("start_res","Initial Resistance (%)", 
                           min = 0, max = 50, step = 1, value = parameters$prob_res*100)
      ),
      
      # Gametocyte group
      menuItem("Gametocyte", tabName = "gametocyte", icon = icon("tint"),
               menuItem("Sensitive", tabName = "gametocyte_sensitive",
               sliderInput("g_is", "Gametocyte (Symptomatic)", 
                           min = 0, max = 10, step = 0.01, value = parameters$g_is[1]),
               sliderInput("g_as", "Gametocyte (Asymptomatic)", 
                           min = 0, max = 10, step = 0.01, value = parameters$g_as[1]),
               sliderInput("g_infs", "Gametocyte (Treatm0ent, Not Failed)", 
                           min = 0, max = 10, step = 0.01, value = parameters$g_infs[1]),
               sliderInput("g_ifs", "Gametocyte (Treatment, Failed)", 
                           min = 0, max = 10, step = 0.01, value = parameters$g_ifs[1])
               ),
               menuItem("Resistant", tabName = "gametocyte_resistant",
               sliderInput("g_ir", "Gametocyte (Symptomatic)", 
                           min = 0, max = 10, step = 0.01, value = parameters$g_ir[1]),
               sliderInput("g_ar", "Gametocyte (Asymptomatic)", 
                           min = 0, max = 10, step = 0.01, value = parameters$g_ar[1]),
               sliderInput("g_infr", "Gametocyte (Treatment, Not Failed)", 
                           min = 0, max = 10, step = 0.01, value = parameters$g_infr[1]),
               sliderInput("g_ifr", "Gametocyte (Treatment, Failed)", 
                           min = 0, max = 10, step = 0.01, value = parameters$g_ifr[1])
               )
      ),
      
      # Treatment (Drugs) group
      menuItem("Treatment", tabName = "treatment", icon = icon("pills"),
               menuItem("Start Time", tabName = "treatment_time",
               sliderInput("D2_start","Start Year of D2", 
                           min = 2020, max = 2033, step = 1, value = 2025,sep = "")
               ),
               menuItem("Clearance Gametocyte Time", tabName = "treatment_clearance",
               sliderInput("D0_treatment","Clearance Time (D0) [days]",min=10,max=120, value = 90, step = 0.01),
               sliderInput("D1_treatment","Clearance Time (D1) [days]",min=1,max=90, value = 44.1, step = 0.01),
               sliderInput("D2_treatment","Clearance Time (D2) [days]",min=0.1,max=60, value = 4.88, step = 0.01)
               ),
               menuItem("Percentage", tabName = "treatment_percentage",
                        sliderInput("D0", "% of No Treatment (D0)", 
                                    min = 0, max = 100, value = 33.34, step = 0.01),
                        sliderInput("D1", "% of ACT (D1)", 
                                    min = 0, max = 100, value = 33.33, step = 0.01),
                        sliderInput("D2", "% of ACT + PQ (D2)", 
                                    min = 0, max = 100, value = 33.33, step = 0.01)
               ),
               menuItem("Treatment Failure", tabName = "treatment_fail",
               sliderInput("fail_rate_s","Sensitive ACT Treatment Failure (%)", 
                                    min = 0, max = 50, step = 0.1, value = parameters$Fail_rate_s*100),
               sliderInput("fail_rate_r","Resistant ACT Treatment Failure (%)", 
                           min = 0, max = 50, step = 0.1, value = parameters$Fail_rate_r*100)
               )

      ),
      
      # Visualization group
      menuItem("Visualization", tabName = "visualization", icon = icon("chart-line"),
               sliderInput("time_x","X-axis Time Range (years)", 
                           min = 2000, max = 2035, step = 1, value = 2014,sep = "")
      )
    )
  ),
  
  dashboardBody(
    useShinyjs(),
    tags$head(
      tags$link(
        rel = "stylesheet", 
        type = "text/css", 
        href = "css/style.css")
    ),
    div(id="intro",
        h2("Model Overview"),
        HTML("<p>We developed a deterministic compartmental model to capture the transmission dynamics of <em>Plasmodium falciparum</em> in a human population, with the goal of evaluating the potential impact of artemisinin-based combination therapy (ACT) and ACT combined with a single low-dose of primaquine (PQ) on gametocyte prevalence and resistance propagation. The model stratifies infections by clinical presentation (asymptomatic or symptomatic), treatment status (D0 = no treatment, D1 = ACT, D2 = ACT+PQ), and parasite phenotype (artemisinin-sensitive vs. artemisinin-resistant).</p>
             
             <p>The model is structured using ordinary differential equations (ODEs) and incorporates gametocyte-producing compartments to explicitly represent the infectious reservoir that contributes to onward transmission. Seasonal variation in transmission is included to approximate vector abundance changes, though vector dynamics are not modelled explicitly.</p>
             
             "),
        
        
          h3("Figure Model"),
          
          # Figure 1
          tags$figure(
            img(src = "img/figure1.png", style = "max-width:100%;"),
            tags$figcaption(
              "Figure 1. Compartmental model of ",
              tags$em("Plasmodium falciparum"),
              " transmission incorporating treatment with ACT and low-dose primaquine. 
      The diagram illustrates the model structure for artemisinin-sensitive strains. 
      The human population is divided into susceptible (S), asymptomatic (A), and symptomatic (I) compartments. 
      Symptomatic infections are stratified by treatment status: untreated (D0), treated with ACT alone (D1), 
      or treated with ACT plus low-dose primaquine (D2). Each infected class may produce gametocytes (G), which 
      contribute to onward transmission. Treated individuals may either recover (ST: success) or experience 
      treatment failure (F), with failed treatments associated with increased gametocyte carriage. Recovered 
      individuals (R) can lose immunity and return to susceptibility. Seasonal variation influences the force 
      of infection."
            )
          ),
          
          br(), br(),
          
          # Figure 2
          tags$figure(
            img(src = "img/figure2.jpg", style = "max-width:100%;"),
            tags$figcaption(
              "Figure 2. Full compartmental model of ",
              tags$em("Plasmodium falciparum"),
              " transmission incorporating treatment with ACT and low-dose primaquine. 
      This complete model includes both artemisinin-sensitive and artemisinin-resistant parasite strains. 
      The human population is divided into susceptible (S), asymptomatic (A), and symptomatic (I) compartments. 
      Symptomatic infections are stratified by treatment status: untreated (D0), treated with ACT alone (D1), 
      or treated with ACT plus low-dose primaquine (D2). Each infected class may produce gametocytes (G), which 
      contribute to onward transmission. Treated individuals may either recover (ST: success) or experience 
      treatment failure (F), with failures associated with higher gametocyte carriage. Recovered individuals (R) 
      can lose immunity and return to susceptibility. Both strains share the same susceptible pool but differ in 
      treatment failure rates and gametocyte dynamics. Seasonal variation modulates the force of infection for both strains."
            )
          ),
        h3("About the Simulation Scenarios"),
        "This app explores how different malaria treatment strategies might affect the spread of artemisinin resistance in Uganda. The model assumes a population of ~24 million (as in 2000) with stable, seasonal malaria transmission. Users can compare three treatment pathways for symptomatic patients:",
        tags$ul(
          tags$li(tags$b("D0:")," No treatment"),
          tags$li(tags$b("D1:")," Artemisinin-based combination therapy (ACT) only"),
          tags$li(tags$b("D2:")," ACT plus low-dose primaquine (0.25 mg/kg)")
        ),
        "The simulations track key outcomes over time, including:",
        tags$ul(
          tags$li("Incidence of sensitive and resistant infections"),
          tags$li("Numbers of asymptomatic and symptomatic cases"),
          tags$li("Gametocyte carriage by strain type")
        ),
        "The model uses demographic data for Uganda (2000–2050) and incorporates changes in transmission patterns before and after 2015, reflecting the scale-up of insecticide-treated nets and indoor residual spraying.

By comparing scenarios with different proportions of D0, D1, and D2, the app highlights the potential role of low-dose primaquine in reducing onward transmission of resistant malaria parasites."
        
    ),
    div(id="simulation",
    tabsetPanel(id="main",
                
                tabPanel(title="Total Incidence",
                         withSpinner(plotOutput("distPlot"))
                         ),
                tabPanel(title="Fraction of transmission & Resistant fraction",
                         withSpinner(plotOutput("Plot_ft")),
                         br(),
                         withSpinner(verbatimTextOutput("kT_star_text")),
                         hr(),
                         br(),
                         withSpinner(plotOutput("Plot_p_res")),
                ),
                tabPanel(title="Symptomatic/Asymptomatic Incidence",
                         withSpinner(plotOutput("Plot_sym_asym")),
                         br(),
                         hr(),
                         br(),
                         withSpinner(plotOutput("Plot_sym_asym_ratio")),
                         ),
                tabPanel(title="Sensitive/Resistant Incidence",
                         withSpinner(plotOutput("Plot_s_r")),
                         br(),
                         hr(),
                         br(),
                         withSpinner(plotOutput("Plot_s_r_ratio"))
                         ),
                tabPanel(title="Gametocyte Infections (Symptomatic vs Asymptomatic)",
                         withSpinner(plotOutput("Plot_Gasym_sym")),
                         br(),
                         hr(),
                         br(),
                         withSpinner(plotOutput("Plot_Gasym_sym_ratio"))
                         ),
                tabPanel(title="Gametocyte Infections (Sensitive vs Resistant)",
                         withSpinner( plotOutput("Plot_GS_GR")),
                         br(),
                         hr(),
                         br(),
                         withSpinner(plotOutput("Plot_GS_GR_ratio"))
                         ),
                tabPanel(title="Treated/untreated channel contribution",
                         withSpinner( plotOutput("Plot_CT_total")),
                         br(),
                         hr(),
                         br(),
                         withSpinner(plotOutput("Plot_CU_total"))
                )
                ),
    ),
                div(id="about_us",
                    tags$b("Model Development Team:"),
                    br(),
                    "Assoc. Prof. Wirichada Pan-ngum (email: ", tags$a(href="mailto:pan@tropmedres.ac", "pan@tropmedres.ac"),")",
                    br(),
                    "Dr. Sompob Saralamba (email: ", tags$a(href="mailto:Sompob@tropmedres.ac", "Sompob@tropmedres.ac"),")",
                    br(),
                    "Mr. Tanaphum Wichaita (email: ", tags$a(href="mailto:Tanaphum@tropmedres.ac", "Tanaphum@tropmedres.ac"),")",
                    br(),
                    br(),
                    tags$b("Scientific and Clinical Advisors:"),
                    br(),
                    "Prof. Arjen Dondorp (email: ", tags$a(href="mailto:arjen@tropmedres.ac", "arjen@tropmedres.ac"),")",
                    br(),
                    "Dr. Thomas J Peto (email: ", tags$a(href="mailto:Tom@tropmedres.ac", "Tom@tropmedres.ac"),")",
                    ),
    div(id="Data",
        h3("Source Data"),
        p("Below is the list of source data and their citations:"),
        tableOutput("source_table")
    ),
    
  )
)

server <- function(session,input, output ) {
  useShinyjs()
  values <- reactiveValues(
    out_in_year = NULL,
    out_in_year_baseline = NULL,
    init_state = init_state,
    parameters = parameters,
    times =times,
    events = events,
    points_time = 2000:2022,
    obs = uganda_Incidence[1:23, 2]
  )
  
  
  observeEvent(input$country,{
    

    if(input$country == "Uganda"){
      
      values$init_state <- init_state
      values$parameters <- parameters
      updateSliderInput(session, "start_res_year", value = 2014)
      
      year_r <- input$start_res_year - 2000
      values$points_time <- 2000:2022
      values$obs <- uganda_Incidence[1:23, 2]
      year_d2 <- input$D2_start - 1999
      values$parameters$start_d <- 12 * year_d2

    }else if(input$country == "Tanzania"){

      values$init_state <- init_state_tanzani
      values$parameters <- parameters_Tanzania
      updateSliderInput(session, "start_res_year", value = 2016)
      
      year_r <- input$start_res_year - 2015
      values$points_time <- 2015:2023
      values$obs <- tanzania_Incidence[6:14,4]
      year_d2 <- input$D2_start - 2014
      values$parameters$start_d <- 12 * year_d2

    }

    values$parameters$start_res_year <- 12 * year_r

    values$events <- list(
      func = function(time, state, parameters) {
        if (time == values$parameters$start_res_year) {
          state["Ir0"] <- state["Ir0"] + state["Is0"] * input$start_res/100
          state["Ir1"] <- state["Ir1"] + state["Is1"] * input$start_res/100
          state["Ir2"] <- state["Ir2"] + state["Is2"] * input$start_res/100
          state["Ar"]  <- state["Ar"]  + state["As"]  * input$start_res/100
          
          state["Is0"] <- state["Is0"] *(1-input$start_res/100)
          state["Is1"] <- state["Is1"] *(1-input$start_res/100)
          state["Is2"] <- state["Is2"] *(1-input$start_res/100)
          state["As"]  <- state["As"]  *(1-input$start_res/100)
        }
        return(state)
      },
      time = values$parameters$start_res_year
    )
    # Transmissibility
    updateSliderInput(session, "Bss", value = values$parameters$beta_s[1])
    updateSliderInput(session, "Bas", value = values$parameters$beta_as)
    updateSliderInput(session, "Bss2", value = values$parameters$beta_s_2[1])
    updateSliderInput(session, "Bas2", value = values$parameters$beta_as_2)
    updateSliderInput(session, "Bsr", value = values$parameters$beta_r[1])
    updateSliderInput(session, "Bar", value = values$parameters$beta_ar)
    
    # Symptomatic
    updateSliderInput(session, "prop_sym", value = 0.25)
    
    # Resistance
    updateSliderInput(session, "start_res", value = values$parameters$prob_res*100)
    
    # Gametocyte
    updateSliderInput(session, "g_is", value = values$parameters$g_is[1])
    updateSliderInput(session, "g_ir", value = values$parameters$g_ir[1])
    updateSliderInput(session, "g_as", value = values$parameters$g_as[1])
    updateSliderInput(session, "g_ar", value = values$parameters$g_ar[1])
    updateSliderInput(session, "g_infs", value = values$parameters$g_infs[1])
    updateSliderInput(session, "g_ifs", value = values$parameters$g_ifs[1])
    updateSliderInput(session, "g_infr", value = values$parameters$g_infr[1])
    updateSliderInput(session, "g_ifr", value = values$parameters$g_ifr[1])
    
    # Treatment
    updateSliderInput(session, "D2_start", value = 2025)
    updateSliderInput(session, "D0_treatment", value = 90)
    updateSliderInput(session, "D1_treatment", value = 44.1)
    updateSliderInput(session, "D2_treatment", value = 4.88)
    updateSliderInput(session, "D0", value = 33.34)
    updateSliderInput(session, "D1", value = 33.33)
    updateSliderInput(session, "D2", value = 33.33)
    updateSliderInput(session, "fail_rate_s", value = values$parameters$Fail_rate_s*100)
    updateSliderInput(session, "fail_rate_r", value = values$parameters$Fail_rate_r*100)
    
    # Visualization
    updateSliderInput(session, "time_x", value = 2014)
  })

  
  observeEvent("", {
    show("intro")
    hide("simulation")
    hide("about_us")
    hide("Data")
  }, once = TRUE)
  
  observeEvent(input$goto_intro,{
    show("intro")
    hide("simulation")
    hide("about_us")
    hide("Data")
  })
  
  observeEvent(input$goto_simulation, {
    hide("intro")
    show("simulation")
    hide("about_us")
    hide("Data")
  })
  
  observeEvent(input$goto_about, {
    hide("intro")
    hide("simulation")
    show("about_us")
    hide("Data")
  })
  observeEvent(input$goto_data, {
    hide("intro")
    hide("simulation")
    hide("about_us")
    show("Data")
  })
  
  observeEvent(input$reset_all, {
    if(input$country == "Uganda"){
    # Transmissibility
    updateSliderInput(session, "Bss", value = parameters$beta_s[1])
    updateSliderInput(session, "Bas", value = parameters$beta_as)
    updateSliderInput(session, "Bss2", value = parameters$beta_s_2[1])
    updateSliderInput(session, "Bas2", value = parameters$beta_as_2)
    updateSliderInput(session, "Bsr", value = parameters$beta_r[1])
    updateSliderInput(session, "Bar", value = parameters$beta_ar)
    
    # Symptomatic
    updateSliderInput(session, "prop_sym", value = 0.25)
    
    # Resistance
    updateSliderInput(session, "start_res_year", value = 2014)
    updateSliderInput(session, "start_res", value = parameters$prob_res*100)
    
    # Gametocyte
    updateSliderInput(session, "g_is", value = parameters$g_is[1])
    updateSliderInput(session, "g_ir", value = parameters$g_ir[1])
    updateSliderInput(session, "g_as", value = parameters$g_as[1])
    updateSliderInput(session, "g_ar", value = parameters$g_ar[1])
    updateSliderInput(session, "g_infs", value = parameters$g_infs[1])
    updateSliderInput(session, "g_ifs", value = parameters$g_ifs[1])
    updateSliderInput(session, "g_infr", value = parameters$g_infr[1])
    updateSliderInput(session, "g_ifr", value = parameters$g_ifr[1])
    
    # Treatment
    updateSliderInput(session, "D2_start", value = 2025)
    updateSliderInput(session, "D0_treatment", value = 90)
    updateSliderInput(session, "D1_treatment", value = 44.1)
    updateSliderInput(session, "D2_treatment", value = 4.88)
    updateSliderInput(session, "D0", value = 33.34)
    updateSliderInput(session, "D1", value = 33.33)
    updateSliderInput(session, "D2", value = 33.33)
    updateSliderInput(session, "fail_rate_s", value = parameters$Fail_rate_s*100)
    updateSliderInput(session, "fail_rate_r", value = parameters$Fail_rate_r*100)
    
    # Visualization
    updateSliderInput(session, "time_x", value = 2014)
    }else{
      # Transmissibility
      updateSliderInput(session, "Bss", value = parameters_Tanzania$beta_s[1])
      updateSliderInput(session, "Bas", value = parameters_Tanzania$beta_as)
      updateSliderInput(session, "Bss2", value = parameters_Tanzania$beta_s_2[1])
      updateSliderInput(session, "Bas2", value = parameters_Tanzania$beta_as_2)
      updateSliderInput(session, "Bsr", value = parameters_Tanzania$beta_r[1])
      updateSliderInput(session, "Bar", value = parameters_Tanzania$beta_ar)
      
      # Symptomatic
      updateSliderInput(session, "prop_sym", value = 0.25)
      
      # Resistance
      updateSliderInput(session, "start_res_year", value = 2016)
      updateSliderInput(session, "start_res", value = parameters_Tanzania$prob_res*100)
      
      # Gametocyte
      updateSliderInput(session, "g_is", value = parameters_Tanzania$g_is[1])
      updateSliderInput(session, "g_ir", value = parameters_Tanzania$g_ir[1])
      updateSliderInput(session, "g_as", value = parameters_Tanzania$g_as[1])
      updateSliderInput(session, "g_ar", value = parameters_Tanzania$g_ar[1])
      updateSliderInput(session, "g_infs", value = parameters_Tanzania$g_infs[1])
      updateSliderInput(session, "g_ifs", value = parameters_Tanzania$g_ifs[1])
      updateSliderInput(session, "g_infr", value = parameters_Tanzania$g_infr[1])
      updateSliderInput(session, "g_ifr", value = parameters_Tanzania$g_ifr[1])
      
      # Treatment
      updateSliderInput(session, "D2_start", value = 2025)
      updateSliderInput(session, "D0_treatment", value = 90)
      updateSliderInput(session, "D1_treatment", value = 44.1)
      updateSliderInput(session, "D2_treatment", value = 4.88)
      updateSliderInput(session, "D0", value = 33.34)
      updateSliderInput(session, "D1", value = 33.33)
      updateSliderInput(session, "D2", value = 33.33)
      updateSliderInput(session, "fail_rate_s", value = parameters_Tanzania$Fail_rate_s*100)
      updateSliderInput(session, "fail_rate_r", value = parameters_Tanzania$Fail_rate_r*100)
      
      # Visualization
      updateSliderInput(session, "time_x", value = 2014)
    }
    
  })
  
  last_mod <- reactiveVal(list(id = "D0", value = 0))
  observeEvent(input$D0, { last_mod(list(id = "D0", value = input$D0)) })
  observeEvent(input$D1, { last_mod(list(id = "D1", value = input$D1)) })
  observeEvent(input$D2, { last_mod(list(id = "D2", value = input$D2)) })
  
  observeEvent(c(input$D0, input$D1, input$D2), {
    d0 <- isolate(input$D0); d1 <- isolate(input$D1); d2 <- isolate(input$D2)
    total <- d0 + d1 + d2
    last_mod_id <- last_mod()$id
    
    if (total > 100) {
      overvalue <- total - 100
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
    } else if (total < 100) {
      remainder <- 100 - total
      if (last_mod_id == "D0") {
        if (d1 + remainder <= 100) {
          updateSliderInput(session, "D1", value = d1 + remainder)
        } else {
          remainder_for_d2 <- remainder - (100 - d1)
          updateSliderInput(session, "D1", value = 100)
          updateSliderInput(session, "D2", value = d2 + remainder_for_d2)
        }
      } else if (last_mod_id == "D1") {
        if (d0 + remainder <= 100) {
          updateSliderInput(session, "D0", value = d0 + remainder)
        } else {
          remainder_for_d2 <- remainder - (100 - d0)
          updateSliderInput(session, "D0", value = 100)
          updateSliderInput(session, "D2", value = d2 + remainder_for_d2)
        }
      } else if (last_mod_id == "D2") {
        if (d0 + remainder <= 100) {
          updateSliderInput(session, "D0", value = d0 + remainder)
        } else {
          remainder_for_d1 <- remainder - (100 - d0)
          updateSliderInput(session, "D0", value = 100)
          updateSliderInput(session, "D1", value = d1 + remainder_for_d1)
        }
      }
    }
    
    values$parameters$propIs <- c(input$D0/100, input$D1/100, input$D2/100)
    values$parameters$propIr <- c(input$D0/100, input$D1/100, input$D2/100)
  })
  
  observeEvent(input$D2_start, {
    if(input$country == "Uganda"){
      year_d2 <- input$D2_start - 1999
      values$parameters$start_d <- 12 * year_d2
    }
    else{
      year_d2 <- input$D2_start - 2014
      values$parameters$start_d <- 12 * year_d2
    }

  })
  
  observeEvent(input$D0_treatment, {
    values$parameters$d0 <- (1/input$D0_treatment) * 30
  })
  
  observeEvent(input$D1_treatment, {
    values$parameters$d1 <- (1/input$D1_treatment) * 30
  })
  
  observeEvent(input$D2_treatment, {
    values$parameters$d2 <- (1/input$D2_treatment) * 30
  })
  
  observeEvent(input$start_res_year, {

    if(input$country == "Uganda"){
      year_r <- input$start_res_year - 2000
    }else if(input$country == "Tanzania"){
      year_r <- input$start_res_year - 2015
    }

    values$parameters$start_res_year <- 12 * year_r

    values$events <- list(
      func = function(time, state, parameters) {
        if (time == values$parameters$start_res_year) {
          state["Ir0"] <- state["Ir0"] + state["Is0"] * input$start_res/100
          state["Ir1"] <- state["Ir1"] + state["Is1"] * input$start_res/100
          state["Ir2"] <- state["Ir2"] + state["Is2"] * input$start_res/100
          state["Ar"]  <- state["Ar"]  + state["As"]  * input$start_res/100
          
          state["Is0"] <- state["Is0"] *(1-input$start_res/100)
          state["Is1"] <- state["Is1"] *(1-input$start_res/100)
          state["Is2"] <- state["Is2"] *(1-input$start_res/100)
          state["As"]  <- state["As"]  *(1-input$start_res/100)
        }
        return(state)
      },
      time = values$parameters$start_res_year
    )
    
  })
  observeEvent(input$start_res, {
    values$parameters$prob_res <- input$start_res / 100

    values$events <- list(
      func = function(time, state, parameters) {
        if (time == values$parameters$start_res_year) {
          state["Ir0"] <- state["Ir0"] + state["Is0"] * input$start_res/100
          state["Ir1"] <- state["Ir1"] + state["Is1"] * input$start_res/100
          state["Ir2"] <- state["Ir2"] + state["Is2"] * input$start_res/100
          state["Ar"]  <- state["Ar"]  + state["As"]  * input$start_res/100
          
          state["Is0"] <- state["Is0"] *(1-input$start_res/100)
          state["Is1"] <- state["Is1"] *(1-input$start_res/100)
          state["Is2"] <- state["Is2"] *(1-input$start_res/100)
          state["As"]  <- state["As"]  *(1-input$start_res/100)
        }
        return(state)
      },
      time = values$parameters$start_res_year
    )
  })
  observeEvent(input$fail_rate_s, {
    values$parameters$Fail_rate_s <- input$fail_rate_s/100
  })
  observeEvent(input$fail_rate_r, {
    values$parameters$Fail_rate_r <- input$fail_rate_r/100
  })
  
  observeEvent(input$Bss, {
    values$parameters$beta_s <- c(input$Bss,input$Bss,input$Bss)
  })
  observeEvent(input$Bas, {
    values$parameters$beta_as <- input$Bas
  })
  observeEvent(input$Bss2, {
    values$parameters$beta_s_2 <- c(input$Bss2,input$Bss2,input$Bss2)
  })
  observeEvent(input$Bas2, {
    values$parameters$beta_as_2 <- input$Bas2
  })
  observeEvent(input$Bsr, {
    values$parameters$beta_r <- c(input$Bsr,input$Bsr,input$Bsr)
  })
  observeEvent(input$Bar, {
    values$parameters$beta_ar <- input$Bar
  })
  
  # Data for Source Data tab
  source_data <- data.frame(
    Description = c(
      "Rate of loss of immunity",
      "Rate of recovery ACT & Primaquine 0.0625 mg/kg",
      "Population of Uganda",
      "Malaria Incidence of Uganda",
      "Population of Tanzania",
      "Malaria Incidence of Tanzania"
    ),
    Citation = c(
      "DOI: 10.1371/journal.pone.0001767",
      "DOI: 10.1002/cpt.2512",
      "https://statisticstimes.com/demographics/country/uganda-population.php",
      "https://data.worldbank.org/indicator/SH.MLR.INCD.P3?name_desc=true&locations=UG",
      "https://statisticstimes.com/demographics/country/tanzania-population.php",
      "https://www.who.int/teams/global-malaria-programme/reports/world-malaria-report-2024"
    ),
    stringsAsFactors = FALSE
  )
  
  output$source_table <- renderTable({
    source_data
  })
  
  observeEvent(input$g_is, {
    values$parameters$g_is <- c(input$g_is,input$g_is,input$g_is)
  })
  observeEvent(input$g_ir, {
    values$parameters$g_ir <- c(input$g_ir,input$g_ir,input$g_ir)
  })
  observeEvent(input$g_as, {
    values$parameters$g_as <- input$g_as
  })
  observeEvent(input$g_ar, {
    values$parameters$g_ar <- input$g_ar
  })
  observeEvent(input$g_infs, {
    values$parameters$g_infs <- c(input$g_infs,input$g_infs)
  })
  observeEvent(input$g_infr, {
    values$parameters$g_infr <- c(input$g_infr,input$g_infr)
  })
  observeEvent(input$g_ifr, {
    values$parameters$g_ifr <- c(input$g_ifr,input$g_ifr)
  })
  observeEvent(input$g_ifs, {
    values$parameters$g_ifs <- c(input$g_ifs,input$g_ifs)
  })
  
  
  observeEvent(input$prop_sym, {
    values$parameters$prob_sym_s <- input$prop_sym
    values$parameters$prob_sym_r <- input$prop_sym
  })
  
  ode_results <- reactive({
    req(input$D0+input$D1+input$D2 ==100 & values$parameters$start_res_year >=0)
    
    out <- ode(
      y = values$init_state,
      times = values$times,
      func = Malaria_model_with_Array,
      parms = values$parameters,
      events = values$events
    )
    vars_all <- colnames(out)
    if(input$country == "Uganda"){
      summarise_by_year(out, vars_all, 36)
    }else{
      summarise_by_year(out, vars_all, 36,start_year=2015)
    }

  })
  
  ode_baseline <- reactive({
    req(input$D0+input$D1+input$D2 ==100 & values$parameters$start_res_year >=0)
    parameters_baseline <- values$parameters
    parameters_baseline$start_d <- 6000
    out <- ode(
      y = values$init_state,
      times = values$times,
      func = Malaria_model_with_Array,
      parms = parameters_baseline,
      events = values$events
    )
    vars_all <- colnames(out)
    if(input$country == "Uganda"){
      summarise_by_year(out, vars_all, 36)
    }else{
      summarise_by_year(out, vars_all, 36,start_year=2015)
    }
    
  })
  
  output$distPlot <- renderPlot({
    out_in_year <- ode_results()
    out_in_year_baseline <- ode_baseline()
    max_y <- max(out_in_year[1:20, "inc"],out_in_year_baseline[, "inc"])

    
    plot(out_in_year[, 1], out_in_year[, "inc"], type = "l",
         xlab = "Years", ylab = "Incidence",
         main = "Total Incidence by Year",
         xlim = c(input$time_x, 2035), lwd = 2, lty = 2,
         col = 1, ylim = c(0, 1.2 * max_y))
    points(values$points_time, values$obs, col = "red", pch = 19)
    lines(out_in_year_baseline[, 1], out_in_year_baseline[, "inc"], lwd = 2)
    abline(v = input$D2_start, lty = 2)
    legend("right", c("Model", "Data"), col = c(1, 2), lty = c(1, NA), pch = c(NA, 19))
    legend("bottomright", c("Baseline","Treatment D2"), col = c(1), lty = c(1, 2),lwd=2)
    text(input$D2_start - 1, 2e07, "Start D2")
  })
  
  output$Plot_sym_asym <- renderPlot({
    out_in_year <- ode_results()
    out_in_year_baseline <- ode_baseline()
    max_y <- max(c(out_in_year[1:35, "inc_sym"], out_in_year[1:35, "inc_asym"],
                   out_in_year_baseline[1:35, "inc_sym"], out_in_year_baseline[1:35, "inc_asym"]
                   ))
    
    plot(out_in_year[, 1], out_in_year[, "inc_sym"], type = "l",
         xlab = "Years", ylab = "Incidence",
         main = "Symptomatic/Asymptomatic Incidence by Year",
         xlim = c(input$time_x, 2035), lwd = 2, lty = 2,
         col = "blue2", ylim = c(0, 1.2 * max_y))
    lines(out_in_year[, 1], out_in_year[, "inc_asym"], col = "green4", lty = 2, lwd = 2)
    lines(out_in_year_baseline[, 1], out_in_year_baseline[, "inc_sym"], col = "blue2", lwd = 2)
    lines(out_in_year_baseline[, 1], out_in_year_baseline[, "inc_asym"], col = "green4", lwd = 2)
    legend("bottomright", c("Symptomatic", "Asymptomatic"), col = c("blue2", "green4"), lty = c(1, 1), lwd = 2)
    abline(v = input$D2_start, lty = 2)
    text(input$D2_start - 1, 2e07, "Start D2")
    legend("topleft", c("Baseline","Treatment D2"), col = c(1), lty = c(1, 2),lwd=2)
  })
  
  output$Plot_sym_asym_ratio <- renderPlot({
    out_in_year <- ode_results()
    out_in_year_baseline <- ode_baseline()

    plot(out_in_year[, 1], out_in_year[, "inc_sym"]/out_in_year[, "inc"], type = "l",
         xlab = "Years", 
         ylab = "Proportion of incidence type",
         main = "Proportion of Symptomatic vs Asymptomatic Incidence Over Time",
         xlim = c(input$time_x, 2035), lwd = 2, lty = 2,
         col = "blue2", ylim = c(0,1))
    lines(out_in_year[, 1], out_in_year[, "inc_asym"]/out_in_year[, "inc"], col = "green4", lty = 2, lwd = 2)
    lines(out_in_year_baseline[, 1], out_in_year_baseline[, "inc_sym"]/out_in_year_baseline[, "inc"], col = "blue2", lwd = 2)
    lines(out_in_year_baseline[, 1], out_in_year_baseline[, "inc_asym"]/out_in_year_baseline[, "inc"], col = "green4", lwd = 2)
    abline(v = input$D2_start, lty = 2)
    legend("bottomright", c("Proportion Symptomatic", "Proportion Asymptomatic"), col = c("blue2", "green4"), lty = c(1, 1), lwd = 2)
    text(input$D2_start - 1, 0.8, "Start D2")
    legend("topright", c("Baseline","Treatment D2"), col = c(1), lty = c(1, 2),lwd=2)
  })
  
  output$Plot_s_r <- renderPlot({
    out_in_year <- ode_results()
    out_in_year_baseline <- ode_baseline()
    max_y <- max(c(out_in_year[1:35, "inc_s"], out_in_year[1:35, "inc_r"],
                   out_in_year_baseline[1:35, "inc_s"], out_in_year_baseline[1:35, "inc_r"]
                   ))
    
    plot(out_in_year[, 1], out_in_year[, "inc_s"], type = "l",
         xlab = "Years", ylab = "Incidence",
         main = "Sensitive/Resistant Incidence by Year",
         xlim = c(input$time_x, 2035), lwd = 2, lty = 2,
         col = "black", ylim = c(0, 1.2 * max_y))
    lines(out_in_year[, 1], out_in_year[, "inc_r"], col = "red", lty = 2, lwd = 2)
    lines(out_in_year_baseline[, 1], out_in_year_baseline[, "inc_r"], col = "red", lwd = 2)
    lines(out_in_year_baseline[, 1], out_in_year_baseline[, "inc_s"], lwd = 2)
    abline(v = input$D2_start, lty = 2)
    legend("bottomright", c("Sensitive", "Resistant"), col = c("black", "red"), lty = c(1, 1), lwd = 2)
    text(input$D2_start - 1, 2e07, "Start D2")
    legend("topleft", c("Baseline","Treatment D2"), col = c(1), lty = c(1, 2),lwd=2)
  })
  
  output$Plot_s_r_ratio <- renderPlot({
    out_in_year <- ode_results()
    out_in_year_baseline <- ode_baseline()
    plot(out_in_year[, 1], out_in_year[, "inc_s"]/out_in_year[, "inc"], type = "l",
         xlab = "Years", 
         ylab = "Proportion of parasite type",
         main = "Proportion of Sensitive vs Resistant Parasites Over Time",
         xlim = c(input$time_x, 2035), lwd = 2, lty = 2,
         col = "black", ylim = c(0,1 ))
    lines(out_in_year[, 1], out_in_year[, "inc_r"]/out_in_year[, "inc"], col = "red", lty = 2, lwd = 2)
    lines(out_in_year_baseline[, 1], out_in_year_baseline[, "inc_r"]/out_in_year_baseline[, "inc"], col = "red", lwd = 2)
    lines(out_in_year_baseline[, 1], out_in_year_baseline[, "inc_s"]/out_in_year_baseline[, "inc"], lwd = 2)
    abline(v = input$D2_start, lty = 2)
    legend("left", c("Proportion Sensitive", "Proportion Resistant"), col = c("black", "red"), lty = c(1, 1), lwd = 2)
    text(input$D2_start - 1, 0.8, "Start D2")
    legend("right", c("Baseline","Treatment D2"), col = c(1), lty = c(1, 2),lwd=2)
  })
  
  output$Plot_Gasym_sym <- renderPlot({
    out_in_year <- ode_results()
    out_in_year_baseline <- ode_baseline()
    max_y <- max(c(out_in_year[1:35, "Gasym_inc"], out_in_year[1:35, "Gsym_inc"],
                   out_in_year_baseline[1:35, "Gasym_inc"], out_in_year_baseline[1:35, "Gsym_inc"]
    ))
    
    plot(out_in_year[, 1], out_in_year[, "Gasym_inc"], type = "l",
         xlab = "Years", ylab = "Incidence",lty = 2,
         main = "Gametocyte Infections Asymptomatic/Symptomatic by Year",
         xlim = c(input$time_x, 2035), lwd = 2, col = 3,
         ylim = c(0, 1.2 * max(out_in_year_baseline[1:35, "Gasym_inc"],out_in_year_baseline[1:35, "Gsym_inc"])))
    lines(out_in_year_baseline[, 1], out_in_year_baseline[, "Gasym_inc"], lty = 1, lwd = 2, col = 3)
    lines(out_in_year[, 1], out_in_year[, "Gsym_inc"], lwd = 2,lty = 2, col = "blue")
    lines(out_in_year_baseline[, 1], out_in_year_baseline[, "Gsym_inc"], lwd = 2, col = "blue")
    abline(v = input$D2_start, lty = 2)
    legend("topleft", c("Asymptomatic", "Symptomatic"), col = c(3, 4), lty = c(1, 1))
    text(input$D2_start - 1, 4e05, "Start D2")
    legend("top", c("Baseline","Treatment D2"), col = c(1), lty = c(1, 2),lwd=2)
  })
  
  output$Plot_Gasym_sym_ratio <- renderPlot({
    out_in_year <- ode_results()
    out_in_year_baseline <- ode_baseline()
    plot(out_in_year[, 1], out_in_year[, "Gasym_inc"]/out_in_year[, "G_inc"], type = "l",
         xlab = "Years", ylab = "Proportion",lty = 2,
         main = "Asymptomatic/Symptomatic Gametocyte Infections Proportion Over Time",
         xlim = c(input$time_x, 2035), lwd = 2, col = 3,
         ylim = c(0, 1))
    lines(out_in_year_baseline[, 1], out_in_year_baseline[, "Gasym_inc"]/out_in_year_baseline[, "G_inc"], lty = 1, lwd = 2, col = 3)
    lines(out_in_year[, 1], out_in_year[, "Gsym_inc"]/out_in_year[, "G_inc"], lwd = 2,lty = 2, col = "blue")
    lines(out_in_year_baseline[, 1], out_in_year_baseline[, "Gsym_inc"]/out_in_year_baseline[, "G_inc"], lwd = 2, col = "blue")
    abline(v = input$D2_start, lty = 2)
    legend("left", c("Proportion Asymptomatic", "Proportion Symptomatic"), col = c(3, 4), lty = c(1, 1))
    text(input$D2_start - 1, 0.8, "Start D2")
    legend("right", c("Baseline","Treatment D2"), col = c(1), lty = c(1, 2),lwd=2)
  })
  
  output$Plot_GS_GR <- renderPlot({
    out_in_year <- ode_results()
    out_in_year_baseline <- ode_baseline()
    max_y <- max(c(out_in_year[1:35, "GR_inc"], out_in_year[1:35, "GR_inc"],
                   out_in_year_baseline[1:35, "GR_inc"], out_in_year_baseline[1:35, "GR_inc"]
    ))
    
    plot(out_in_year[, 1], out_in_year[, "GS_inc"], type = "l",
         xlab = "Years", ylab = "Incidence",
         main = "Gametocyte Infections Sensitive/Resistant by Year",
         xlim = c(input$time_x, 2035), lwd = 2, col = 1,lty = 2,
         ylim = c(0, 1.2 * max_y))
    lines(out_in_year_baseline[, 1], out_in_year_baseline[, "GS_inc"], lty = 1, lwd = 2)
    lines(out_in_year[, 1], out_in_year[, "GR_inc"], lwd = 2,lty = 2, col = "red")
    lines(out_in_year_baseline[, 1], out_in_year_baseline[, "GR_inc"], lwd = 2, col = "red")
    abline(v = input$D2_start, lty = 2)
    legend("topleft", c("Sensitive", "Resistant"), col = c(1, 2), lty = c(1, 1))
    text(input$D2_start - 1, 4e05, "Start D2")
    legend("top", c("Baseline","Treatment D2"), col = c(1), lty = c(1, 2),lwd=2)
  })
  
  output$Plot_GS_GR_ratio <- renderPlot({
    out_in_year <- ode_results()
    out_in_year_baseline <- ode_baseline()
    plot(out_in_year[, 1], out_in_year[, "GS_inc"]/out_in_year[, "G_inc"], type = "l",
         xlab = "Years", ylab = "Proportion",
         main = "Sensitive/Resistant Gametocyte Infections Proportion Over Time",
         xlim = c(input$time_x, 2035), lwd = 2, col = 1,lty = 2,
         ylim = c(0, 1))
    lines(out_in_year_baseline[, 1], out_in_year_baseline[, "GS_inc"]/out_in_year_baseline[, "G_inc"], lty = 1, lwd = 2)
    lines(out_in_year[, 1], out_in_year[, "GR_inc"]/out_in_year[, "G_inc"], lwd = 2,lty = 2, col = "red")
    lines(out_in_year_baseline[, 1], out_in_year_baseline[, "GR_inc"]/out_in_year_baseline[, "G_inc"], lwd = 2, col = "red")
    abline(v = input$D2_start, lty = 2)
    legend("left", c("Proportion Sensitive", "Proportion Resistant"), col = c(1, 2), lty = c(1, 1))
    text(input$D2_start - 1, 0.8, "Start D2")
    legend("right", c("Baseline","Treatment D2"), col = c(1), lty = c(1, 2),lwd=2)
  })
  
  output$Plot_ft <- renderPlot({
    out_in_year <- ode_results()
    out_in_year_baseline <- ode_baseline()
    max_y <- max(c(out_in_year[1:35, "F_T"]
    ))
    
    out_in_year <- ode_results()
    out_in_year_baseline <- ode_baseline()
    plot(out_in_year[, 1], out_in_year[, "F_T"], type = "l",
         xlab = "Years", ylab = "Proportion",
         main = "FT",
         xlim = c(input$time_x, 2035), lwd = 2, col = 1,lty = 2,
         ylim = c(0, max_y*1.2))
    lines(out_in_year_baseline[, 1], out_in_year_baseline[, "F_T"], lty = 1, lwd = 2)
    abline(v = input$D2_start, lty = 2)
    text(input$D2_start - 1, 0.8, "Start D2")
    legend("left", c("Baseline","Treatment D2"), col = c(1), lty = c(1, 2),lwd=2)
  })
  
  output$Plot_ft <- renderPlot({
    out_in_year <- ode_results()
    out_in_year_baseline <- ode_baseline()
    max_y <- max(c(out_in_year[1:35, "F_T"]
    ))
    
    out_in_year <- ode_results()
    out_in_year_baseline <- ode_baseline()
    plot(out_in_year[, 1], out_in_year[, "F_T"], type = "l",
         xlab = "Years", ylab = "FT",
         main = "FT",
         xlim = c(input$time_x, 2035), lwd = 2, col = 1,lty = 2,
         ylim = c(0, max_y*1.2))
    lines(out_in_year_baseline[, 1], out_in_year_baseline[, "F_T"], lty = 1, lwd = 2)
    abline(v = input$D2_start, lty = 2)
    text(input$D2_start - 1, 0.8, "Start D2")
    legend("left", c("Baseline","Treatment D2"), col = c(1), lty = c(1, 2),lwd=2)
  })
  
  output$Plot_p_res <- renderPlot({
    out_in_year <- ode_results()
    out_in_year_baseline <- ode_baseline()
    max_y <- max(c(out_in_year[1:35, "p_resistant_inc"]
    ))
    
    
    out_in_year <- ode_results()
    out_in_year_baseline <- ode_baseline()
    plot(out_in_year[, 1], out_in_year[, "p_resistant_inc"], type = "l",
         xlab = "Years", ylab = "Resistant Infection Fraction",
         main = "Resistant fraction of incident infections",
         xlim = c(input$time_x, 2035), lwd = 2, col = 1,lty = 2,
         ylim = c(0, 1))
    lines(out_in_year_baseline[, 1], out_in_year_baseline[, "p_resistant_inc"], lty = 1, lwd = 2)
    abline(v = input$D2_start, lty = 2)
    text(input$D2_start - 1, 0.8, "Start D2")
    legend("left", c("Baseline","Treatment D2"), col = c(1), lty = c(1, 2),lwd=2)
  })
  
  output$Plot_CT_total <- renderPlot({
    out_in_year <- ode_results()
    out_in_year_baseline <- ode_baseline()
    max_y <- max(c(out_in_year[1:35, "CT_total"],out_in_year_baseline[, "CT_total"]
    ))
    
    out_in_year <- ode_results()
    out_in_year_baseline <- ode_baseline()
    plot(out_in_year[, 1], out_in_year[, "CT_total"], type = "l",
         xlab = "Years", ylab = "Treated Contribution",
         main = "CT_total",
         xlim = c(input$time_x, 2035), lwd = 2, col = 1,lty = 2,
         ylim = c(0, max_y*1.2))
    lines(out_in_year_baseline[, 1], out_in_year_baseline[, "CT_total"], lty = 1, lwd = 2)
    abline(v = input$D2_start, lty = 2)
    text(input$D2_start - 1, 0.8, "Start D2")
    legend("bottomleft", c("Baseline","Treatment D2"), col = c(1), lty = c(1, 2),lwd=2)
  })
  
  output$Plot_CU_total <- renderPlot({
    out_in_year <- ode_results()
    out_in_year_baseline <- ode_baseline()
    max_y <- max(c(out_in_year[1:35, "CU_total"],out_in_year_baseline[, "CU_total"]
    ))
    
    out_in_year <- ode_results()
    out_in_year_baseline <- ode_baseline()
    plot(out_in_year[, 1], out_in_year[, "CU_total"], type = "l",
         xlab = "Years", ylab = "Untreated Contribution",
         main = "CU_total",
         xlim = c(input$time_x, 2035), lwd = 2, col = 1,lty = 2,
         ylim = c(0, max_y*1.2))
    lines(out_in_year_baseline[, 1], out_in_year_baseline[, "CU_total"], lty = 1, lwd = 2)
    abline(v = input$D2_start, lty = 2)
    text(input$D2_start - 1, 0.8, "Start D2")
    legend("bottomleft", c("Baseline","Treatment D2"), col = c(1), lty = c(1, 2),lwd=2)
  })
  
  output$kT_star_text <- renderText({
    out_in_year <- ode_results()
    keep <- out_in_year[, "F_T"] 
    FT_mean <- mean(keep, na.rm = TRUE)

    # implied treated-channel advantage
    kT_star <- 1 + input$s_hat / FT_mean
    paste("Treatment-associated transmission advantage:", round(kT_star, 2))
  })

}

shinyApp(ui = ui, server = server)
