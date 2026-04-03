library(shiny)
library(shinydashboard)
library(htmltools)
library(shinyBS)
library(deSolve)
library(shinyjs)
library(bsplus)
library(shinycssloaders)

# Read model
source("model_D0_33_D1_67_baseline_change_birth_death_2betaSets_SS.R")
init_state <- readRDS("init_state_d1_33.rds")
source("init_parameter_equilibrium_birth_2_rate_2beta_res_fitted_Tanzania.R")
tanzania_Incidence <- read.csv("data/Tanzania_incidence.csv")
tanzania_Population <- read.csv("data/Tanzania_pop_2000_2100.csv", header = TRUE)[1:51, 2]
init_state_tanzania <- readRDS("init_state_d1_33.rds")
parameters <- parameters_Tanzania

events <- list(
  func = function(time, state, parameters) {
    if (time == 1) {
      state["Is0"] <- state["Is0"] - 1
      state["Is1"] <- state["Is1"] - 1
      state["As"]  <- state["As"]  - 1
      
      state["Ir0"] <- state["Ir0"] + 1
      state["Ir1"] <- state["Ir1"] + 1
      state["Ar"]  <- state["Ar"]  + 1
      
      
    }
    return(state)
  },
  time = 1
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
        else if (v == "p_resistant_inc" |
                 v == "F_T" | v == "CT_total" | v == "CU_total") {
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
  dashboardHeader(
    title = "MalariaPQ",
    titleWidth = "300px",
    tags$li(
      actionLink("goto_intro", HTML("<b>Introduction</b>"), class = "btn btn-default"),
      class = "dropdown"
    ),
    tags$li(
      actionLink("goto_simulation", HTML("<b>Simulation</b>"), class = "btn btn-default"),
      class = "dropdown"
    ),
    tags$li(
      actionLink("goto_about", HTML("<b>About us</b>"), class = "btn btn-default"),
      class = "dropdown"
    ),
    tags$li(
      actionLink("goto_data", HTML("<b>Source Data</b>"), class = "btn btn-default"),
      class = "dropdown"
    )
  ),
  dashboardSidebar(
    width = "300px",
    sidebarMenu(
      # useShinyjs(),
      actionButton("reset_all", "Reset", icon = icon("undo")),
      div(style = "display:none;", radioButtons(
        "country" , "Data Country", choices = c("Tanzania")
      )),
      
      # Transmissibility group
      menuItem(
        "Transmissibility",
        tabName = "transmissibility",
        icon = icon("bug"),
        bsPopover(
          id = "Sensitive_2000_2014",
          title = NULL,
          content = "β is defined separately for sensitive/resistant strains. Sensitive infections are split into pre- and post-.",
          placement = "right",
          trigger = "hover"
        ),
        menuItem(
          tagList(tags$span(
            "Sensitive", tags$i(class = "fa fa-circle-info"), id = "Sensitive_2000_2014"
          )),
          tabName = "transmissibility_sensitive",
          bsPopover(
            id = "beta1",
            title = NULL,
            content = "Tanzania: β1 (2015–2020)",
            placement = "right",
            trigger = "hover"
          ),
          menuItem(
            tagList(tags$span(
              "Beta1", tags$i(class = "fa fa-circle-info"), id = "beta1"
            )),
            tabName = "Beta_1",
            sliderInput(
              "Bss",
              "β (Symptomatic)",
              min = 0,
              max = 10,
              step = 0.01,
              value = parameters$beta_s[1]
            ),
            sliderInput(
              "Bas",
              "β (Asymptomatic)",
              min = 0,
              max = 10,
              step = 0.01,
              value = parameters$beta_as
            )
          ),
          bsPopover(
            id = "beta2",
            title = NULL,
            content = "Tanzania: β2 (2020–2035)",
            placement = "right",
            trigger = "hover"
          ),
          menuItem(
            tagList(tags$span(
              "Beta2", tags$i(class = "fa fa-circle-info"), id = "beta2"
            )),
            tabName = "2014–2035",
            sliderInput(
              "Bss2",
              "β (Symptomatic)",
              min = 0,
              max = 10,
              step = 0.01,
              value = parameters$beta_s_2[1]
            ),
            sliderInput(
              "Bas2",
              "β (Asymptomatic)",
              min = 0,
              max = 10,
              step = 0.01,
              value = parameters$beta_as_2
            )
          )
        ),
        menuItem(
          "Resistant",
          tabName = "transmissibility_resistant",
          sliderInput(
            "M_beta_resistant",
            "Multiplier Beta Resistant",
            min = 0,
            max = 3,
            step = 0.01,
            value = 1.131
          ),
          sliderInput(
            "C_beta_resistant",
            "Resistant Fitness Cost",
            min = 0,
            max = 1,
            step = 0.01,
            value = parameters$c_beta_r
          )
        )
      ),
      
      # Symptomatic group
      menuItem(
        "Symptomatic",
        tabName = "symptomatic",
        icon = icon("user-md"),
        sliderInput(
          "prop_sym",
          "Proportion Symptomatic",
          min = 0.01,
          max = 1,
          step = 0.01,
          value = 0.25
        )
      ),
      
      # Resistance group
      menuItem(
        "Resistance",
        tabName = "resistance",
        icon = icon("dna"),
        sliderInput(
          "start_res_year",
          "Start Year of Resistance",
          min = 2000,
          max = 2025,
          step = 1,
          value = 2000,
          sep = ""
        )
      ),
      
      # Gametocyte group
      menuItem(
        "Gametocyte",
        tabName = "gametocyte",
        icon = icon("tint"),
        menuItem(
          "Sensitive",
          tabName = "gametocyte_sensitive",
          sliderInput(
            "g_is",
            "Gametocyte (Symptomatic)",
            min = 0,
            max = 10,
            step = 0.01,
            value = parameters$g_is[1]
          ),
          sliderInput(
            "g_as",
            "Gametocyte (Asymptomatic)",
            min = 0,
            max = 10,
            step = 0.01,
            value = parameters$g_as[1]
          ),
          sliderInput(
            "g_infs",
            "Gametocyte (Treatm0ent, Not Failed)",
            min = 0,
            max = 10,
            step = 0.01,
            value = parameters$g_infs[1]
          ),
          sliderInput(
            "g_ifs",
            "Gametocyte (Treatment, Failed)",
            min = 0,
            max = 10,
            step = 0.01,
            value = parameters$g_ifs[1]
          )
        ),
        menuItem(
          "Resistant",
          tabName = "gametocyte_resistant",
          sliderInput(
            "g_ir",
            "Gametocyte (Symptomatic)",
            min = 0,
            max = 10,
            step = 0.01,
            value = parameters$g_ir[1]
          ),
          sliderInput(
            "g_ar",
            "Gametocyte (Asymptomatic)",
            min = 0,
            max = 10,
            step = 0.01,
            value = parameters$g_ar[1]
          ),
          sliderInput(
            "g_infr",
            "Gametocyte (Treatment, Not Failed)",
            min = 0,
            max = 10,
            step = 0.01,
            value = parameters$g_infr[1]
          ),
          sliderInput(
            "g_ifr",
            "Gametocyte (Treatment, Failed)",
            min = 0,
            max = 10,
            step = 0.01,
            value = parameters$g_ifr[1]
          )
        )
      ),
      
      # Treatment (Drugs) group
      menuItem(
        "Treatment",
        tabName = "treatment",
        icon = icon("pills"),
        menuItem(
          "Start Time",
          tabName = "treatment_time",
          sliderInput(
            "D2_start",
            "Start Year of D2",
            min = 2000,
            max = 2033,
            step = 1,
            value = 2025,
            sep = ""
          )
        ),
        menuItem(
          "Clearance Gametocyte Time",
          tabName = "treatment_clearance",
          sliderInput(
            "D0_treatment",
            "Clearance Time (D0) [days]",
            min = 10,
            max = 120,
            value = 90,
            step = 0.01
          ),
          sliderInput(
            "D1_treatment",
            "Clearance Time (D1) [days]",
            min = 1,
            max = 90,
            value = 44.1,
            step = 0.01
          ),
          sliderInput(
            "D2_treatment",
            "Clearance Time (D2) [days]",
            min = 0.1,
            max = 60,
            value = 4.88,
            step = 0.01
          )
        ),
        menuItem(
          "Percentage",
          tabName = "treatment_percentage",
          sliderInput(
            "D0",
            "% of No Treatment (D0)",
            min = 0,
            max = 100,
            value = 33.34,
            step = 0.01
          ),
          sliderInput(
            "D1",
            "% of ACT (D1)",
            min = 0,
            max = 100,
            value = 33.33,
            step = 0.01
          ),
          sliderInput(
            "D2",
            "% of ACT + PQ (D2)",
            min = 0,
            max = 100,
            value = 33.33,
            step = 0.01
          )
        ),
        menuItem(
          "Treatment Failure",
          tabName = "treatment_fail",
          sliderInput(
            "fail_rate_s",
            "Sensitive ACT Treatment Failure (%)",
            min = 0,
            max = 50,
            step = 0.1,
            value = parameters$Fail_rate_s * 100
          ),
          sliderInput(
            "fail_rate_r",
            "Resistant ACT Treatment Failure (%)",
            min = 0,
            max = 50,
            step = 0.1,
            value = parameters$Fail_rate_r * 100
          )
        )
        
      ),
      div(
        style = "display:none;",
        # Visualization group
        menuItem(
          "Visualization",
          tabName = "visualization",
          icon = icon("chart-line"),
          sliderInput(
            "time_x",
            "X-axis Time Range (years)",
            min = 2000,
            max = 2035,
            step = 1,
            value = 2000,
            sep = ""
          )
        )
      )
    )
    
  ),
  
  dashboardBody(
    useShinyjs(),
    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "css/style.css")
    ),
    tags$head(
      tags$link(
        rel = "stylesheet", 
        type = "text/css", 
        href = "css/style.css"),
      tags$script(src = "js/toTheTop.js")
    ),
      div(class = "flag-accent"),
    div(
      id = "intro",
      div(
        class = "general-card",
        h3("Model Overview"),
        HTML(
          "<p>We developed a deterministic compartmental model to capture the transmission dynamics of Plasmodium falciparum in a human population, with the goal of evaluating the potential impact of artemisinin-based combination therapy (ACT) and ACT combined with a single low-dose of primaquine (PQ) on gametocyte prevalence and resistance propagation. The model stratifies infections by clinical presentation (asymptomatic or symptomatic), treatment status (no treatment, ACT, or ACT+PQ), and parasite phenotype (artemisinin-sensitive vs. artemisinin-resistant).</p>

             <p>The model is structured using ordinary differential equations (ODEs) and incorporates gametocyte-producing compartments to explicitly represent the infectious reservoir that contributes to onward transmission. Seasonal variation in transmission is included to approximate vector abundance changes, though vector dynamics are not modelled explicitly.</p>

             "
        ),
        
        
        h3("Model Compartment Descriptions"),
        tags$p(
          "The model divides the human population into several compartments to capture the infection dynamics of ",
          tags$i("Plasmodium falciparum"),
          ", stratified by clinical status, treatment pathway, and parasite resistance phenotype."
        ),
        
        tags$p(
          "Susceptible individuals (S) represent the uninfected population at risk of acquiring malaria through exposure to infectious mosquito bites. ",
          "Upon infection, individuals either develop asymptomatic or symptomatic infections, depending on their immune status and exposure history."
        ),
        
        tags$p(
          "Asymptomatic infections are represented by ",
          tags$b("A" , tags$sub("s")),
          " and ",
          tags$b("A" , tags$sub("r"), ),
          ", which denote individuals infected with artemisinin-sensitive and artemisinin-resistant parasites, respectively. ",
          "These individuals do not seek treatment but may still contribute to transmission by developing transmissible levels of gametocytes. ",
          "Their corresponding gametocyte-producing compartments are ",
          tags$b("GA" , tags$sub("s"), ),
          " and ",
          tags$b("GA" , tags$sub("r"), ),
          "."
        ),
        
        tags$p(
          "Symptomatic infections are further stratified by treatment received and parasite resistance status. ",
          "For artemisinin-sensitive strains, ",
          tags$b("I" , tags$sub("s0")),
          ", ",
          tags$b("I" , tags$sub("s1")),
          ", and ",
          tags$b("I" , tags$sub("s2")),
          " represent individuals receiving no treatment (D0), ACT only (D1), and ACT plus low-dose primaquine (D2), respectively. ",
          "Similarly, ",
          tags$b("I" , tags$sub("r0")),
          ", ",
          tags$b("I", tags$sub("r1")),
          ", and ",
          tags$b("I", tags$sub("r2")),
          " represent symptomatic infections caused by resistant parasites under the same treatment categories."
        ),
        
        tags$p(
          "Each treatment pathway is followed by one of two possible outcomes: success or failure. ",
          "Successful treatment outcomes are tracked in ",
          tags$b("ST", tags$sub("s1")),
          ", ",
          tags$b("ST", tags$sub("s2")),
          ", ",
          tags$b("ST", tags$sub("r1")),
          ", and ",
          tags$b("ST", tags$sub("r2")),
          ", while treatment failures are recorded in ",
          tags$b("F", tags$sub("s1")),
          ", ",
          tags$b("F", tags$sub("s2")),
          ", ",
          tags$b("F", tags$sub("r1")),
          ", and ",
          tags$b("F", tags$sub("r2")),
          "."
        ),
        
        tags$p(
          "Gametocyte production from symptomatic infections is represented by ",
          tags$b("GI", tags$sub("s0")),
          ", ",
          tags$b("GI", tags$sub("s1")),
          ", and ",
          tags$b("GI", tags$sub("s2")),
          " for sensitive infections, and ",
          tags$b("GI", tags$sub("r0")),
          ", ",
          tags$b("GI", tags$sub("r1")),
          ", and ",
          tags$b("GI", tags$sub("r2")),
          " for resistant infections."
        ),
        
        tags$p(
          "Finally, individuals who clear infection either naturally or through treatment enter the recovered compartments: ",
          tags$b("R", tags$sub("s")),
          " for sensitive infections and ",
          tags$b("R", tags$sub("r")),
          " for resistant infections."
        ),
        
        # Figure 1
        tags$figure(
          tags$div(
            class = "img-container",
            
            tags$img(
              src = "img/figure1.png",
              class = "zoom-img",
              onclick = "document.getElementById('imgModal').style.display='block'"
            ),
            
            tags$div(
              id = "imgModal",
              class = "img-modal",
              onclick = "this.style.display='none'",
              
              tags$img(src = "img/figure1.png", class = "img-modal-content")
            )
          ),
          tags$figcaption(
            tags$b(
              "Figure 1. Compartmental model of ",
              tags$em("Plasmodium falciparum"),
              " transmission incorporating treatment with ACT and low-dose primaquine."
            ),
            "The diagram illustrates the model structure for artemisinin-sensitive strains.
      The human population is divided into susceptible (S), asymptomatic (A), and symptomatic (I) compartments.
      Symptomatic infections are stratified by treatment status: untreated (D0), treated with ACT alone (D1),
      or treated with ACT plus low-dose primaquine (D2). Each infected class may produce gametocytes (G), which
      contribute to onward transmission. Treated individuals may either recover (ST: success) or experience
      treatment failure (F), with failed treatments associated with increased gametocyte carriage. Recovered
      individuals (R) can lose immunity and return to susceptibility. Seasonal variation influences the force
      of infection."
          )
        ),
      ),
      div(
        class = "general-card",
        h3("Transmission Dynamics"),
        tags$p(
          "Transmission in the model is governed by the force of infection arising from infectious gametocyte carriers in the human population. The total force of infection is decomposed into contributions from four sources: symptomatic and asymptomatic infections with artemisinin-sensitive parasites, and symptomatic and asymptomatic infections with artemisinin-resistant parasites."
        ),
        tags$p("Specifically, the force of infection components are defined as:"),
        withMathJax(),
        div(
          class = "eq-box",
        HTML(
          "$$season(t) = flo + amp \\sin\\left(\\frac{2\\pi}{period}(phase + t)\\right)$$"
        ),
        HTML(
          "
  $$\\begin{aligned}
  \\lambda_{is} &= season(t) \\frac{\\beta_s \\left(\\sum_{i=0}^{2} GI_{s,i}\\right)}{N(t)},
  &\\lambda_{as} &= season(t) \\frac{\\beta_{as} GA_s}{N(t)} \\\\
  \\lambda_{ir} &= season(t) \\frac{\\beta_r c_{\\beta r} \\left(\\sum_{i=0}^{2} GI_{r,i}\\right)}{N(t)},
  &\\lambda_{ar} &= season(t) \\frac{\\beta_{ar} c_{\\beta r} GA_r}{N(t)}
  \\end{aligned}$$
  "
        )
        ),
        tags$p("where:"),
        # Text + subscripts (NO MathJax)
        tags$ul(
          tags$li(
            tags$b("GI", tags$sub("s,i")),
            " and ",
            tags$b("GI", tags$sub("r,i")),
            " denote gametocyte carriers originating from symptomatic infections under treatment pathway i ∈ {0,1,2}"
          ),
          tags$li(
            tags$b("GA", tags$sub("s")),
            " and ",
            tags$b("GA", tags$sub("r")),
            " denote gametocyte carriers from asymptomatic infections"
          ),
          tags$li(
            tags$b("β", tags$sub("s")),
            " and ",
            tags$b("β", tags$sub("r")),
            " are transmission coefficients for sensitive and resistant parasites"
          ),
          tags$li(
            tags$b("c", tags$sub("βr")),
            " represents the transmission fitness cost of resistant parasites"
          ),
          tags$li(tags$b("N"), " is the total population size"),
          tags$li(
            tags$b("season(t)"),
            " is the time-varying seasonal forcing function"
          )
        ),
        tags$p("The total force of infection is given by:"),
        # Equation inline (MathJax)
        div(
          class = "eq-box",
        tags$p(
          HTML(
            "$$\\lambda = \\lambda_{is} + \\lambda_{as} + \\lambda_{ir} + \\lambda_{ar}$$"
          )
        )
        ),
        tags$p(
          "New infections occur at rate ",
          tags$b("λS"),
          ", and are partitioned into symptomatic and asymptomatic infections according to strain-specific probabilities (",
          tags$b("p", tags$sub("syms")),
          ", ",
          tags$b("p", tags$sub("syma")),
          ")."
        ),
        
        tags$p(
          "Symptomatic infections are further distributed into treatment pathways (D0, D1, D2), while asymptomatic infections remain untreated."
        ),
        
        tags$p(
          "This formulation ensures that transmission depends on both the infectious reservoir composition and treatment structure."
        ),
      ),
      tags$div(
        class = "general-card",
        tags$h3("Fitness cost of resistant parasites"),
        
        tags$p(
          "The model incorporates a transmission fitness cost for resistant parasites through the parameter ",
          tags$b("c", tags$sub("βr")),
          ", which scales the contribution of resistant infections to the force of infection. ",
          "This is implemented directly in the expressions for ",
          tags$b("λ", tags$sub("ir")),
          " and ",
          tags$b("λ", tags$sub("ar")),
          ":"
        ),
        div(
          class = "eq-box",
        tags$p(
          HTML(
            "$$\\lambda_{ir} \\propto \\beta_r c_{\\beta r}, \\quad \\lambda_{ar} \\propto \\beta_{ar} c_{\\beta r}$$"
          )
        )
        ),
        
        tags$p(
          "Thus, the effective transmission coefficients for resistant parasites become:"
        ),
        div(
          class = "eq-box",
        tags$p(
          HTML(
            "$$\\beta_r^{eff} = c_{\\beta r}\\beta_r, \\quad \\beta_{ar}^{eff} = c_{\\beta r}\\beta_{ar}$$"
          )
        )
        ),
        
        tags$p(
          "When ",
          tags$b("c", tags$sub("βr")),
          " < 1, resistant parasites have reduced transmissibility relative to sensitive parasites in the absence of treatment. ",
          "This reduction applies uniformly across both symptomatic and asymptomatic transmission pathways, reflecting the assumption that the fitness cost is an intrinsic property of the resistant parasite rather than dependent on clinical presentation."
        ),
        
        tags$p(
          "This parameter plays a key role in determining the rate of spread of resistance, as it counterbalances the transmission advantage of resistant parasites under drug pressure."
        )
      ),
      tags$div(
        class = "general-card",
        tags$h3("Recovery and loss of immunity"),
        
        tags$p(
          "Recovery in the model occurs through multiple pathways and contributes to the accumulation of individuals in the recovered compartments ",
          tags$b("R", tags$sub("s")),
          " and ",
          tags$b("R", tags$sub("r")),
          ". ",
          "These compartments represent individuals who have cleared infection with sensitive or resistant parasites, respectively."
        ),
        tags$p(
          "The general structure of the recovery dynamics can be expressed as:"
        ),
        div(
          class = "eq-box",
        tags$p(
          HTML(
            "$$\\frac{dR_s}{dt} = (\\text{recovery from sensitive infections}) - \\alpha R_s - \\mu_o R_s$$"
          )
        ),
        tags$p(
          HTML(
            "$$\\frac{dR_r}{dt} = (\\text{recovery from resistant infections}) - \\alpha R_r - \\mu_o R_r$$"
          )
        ),
        ),
        tags$p("Recovery terms include:"),
        
        tags$ul(
          tags$li(
            "clearance of asymptomatic infections (",
            tags$b("τ", tags$sub("as")),
            " ",
            tags$b("A", tags$sub("s")),
            ", ",
            tags$b("τ", tags$sub("ar")),
            " ",
            tags$b("A", tags$sub("r")),
            ")"
          ),
          tags$li(
            "clearance of gametocyte compartments (",
            tags$b("τ", tags$sub("is,i")),
            " ",
            tags$b("GI", tags$sub("s,i")),
            ", ",
            tags$b("τ", tags$sub("ir,i")),
            " ",
            tags$b("GI", tags$sub("r,i")),
            ", ",
            tags$b("τ", tags$sub("ags")),
            " ",
            tags$b("GA", tags$sub("s")),
            ", ",
            tags$b("τ", tags$sub("agr")),
            " ",
            tags$b("GA", tags$sub("r")),
            ")"
          ),
          tags$li(
            "recovery from untreated symptomatic infections (",
            tags$b("τ", tags$sub("ntsd0")),
            " ",
            tags$b("I", tags$sub("s,0")),
            ", ",
            tags$b("τ", tags$sub("ntrd0")),
            " ",
            tags$b("I", tags$sub("r,0")),
            ")"
          ),
          tags$li(
            "recovery following treatment success and failure (",
            tags$b("τ", tags$sub("nfs")),
            " ",
            tags$b("ST"),
            ", ",
            tags$b("τ", tags$sub("fs")),
            " ",
            tags$b("F"),
            ")"
          )
        ),
        
        tags$p(
          "Recovered individuals lose immunity at rate ",
          tags$b("α"),
          " and return to the susceptible compartment:"
        ),
        div(
          class = "eq-box",
        tags$p(
          HTML(
            "$$\\frac{dS}{dt} = \\mu_i N - \\mu_o S - \\lambda S + \\alpha (R_s + R_r)$$"
          )
        )
        ),
        
        tags$p(
          "Here ",
          tags$b("μ", tags$sub("i")),
          " and ",
          tags$b("μ", tags$sub("o")),
          " denote birth and death rates, respectively. ",
          "This SIRS structure allows individuals to experience repeated infections over time, reflecting the partial and waning immunity characteristic of malaria in endemic settings."
        )
      ),
      tags$div(
        class = "general-card",
        tags$h3("Model Equations"),
        
        tags$p(
          "Equations for each compartment (S, A, I, ST, F, GI, GA, R) are defined below."
        ),
        div(class = "row eq-box", div(
          class = "col-sm-12", tags$h4("Susceptible"), tags$p(
            HTML(
              "$$\\frac{dS}{dt} = \\mu_i N - \\mu_o S - \\lambda S + \\alpha R_s + \\alpha R_r$$"
            )
          ),
        ), ),
        div(
          class = "row eq-box",
          div(class = "col-md-6", tags$h4("Asymptomatic (Sensitive)"), tags$p(
            HTML(
              "$$\\frac{dA_s}{dt} = p\\lambda S (1 - p_s) - \\gamma_{as} A_s - \\tau_{as} A_s - \\mu_o A_s$$"
            )
          ), ),
          div(class = "col-md-6", tags$h4("Asymptomatic (Resistant)"), tags$p(
            HTML(
              "$$\\frac{dA_r}{dt} = (1-p)\\lambda S (1 - p_r) - \\gamma_{ar} A_r - \\tau_{ar} A_r - \\mu_o A_r$$"
            )
          ), ),
        ),
        div(
          class = "row eq-box",
          div(
            class = "col-md-6",
            tags$h4("Gametocytes from Asymptomatic (Sensitive)"),
            tags$p(
              HTML(
                "$$\\frac{dGA_s}{dt} = \\gamma_{as} A_s - \\tau_{ags} GA_s - \\mu_o GA_s$$"
              )
            ),
          ),
          div(
            class = "col-md-6",
            tags$h4("Gametocytes from Asymptomatic (Resistant)"),
            tags$p(
              HTML(
                "$$\\frac{dGA_r}{dt} = \\gamma_{ar} A_r - \\tau_{agr} GA_r - \\mu_o GA_r$$"
              )
            ),
          ),
        ),
        div(
          class = "row eq-box",
          div(
            class = "col-sm-12",
            tags$h4("Symptomatic (Sensitive), k = 0,1,2"),
            tags$p(
              HTML(
                "$$\\frac{dI_{s,k}}{dt} = (\\lambda_{is,k} S + \\lambda_{as} S)p_{sym_s}p_{is,k} - \\gamma_{is} I_{s,k} - \\mu_o I_{s,k}$$"
              )
            ),
          ),
          div(class = "col-sm-6", tags$p("For k = 0"), tags$p(
            HTML(
              "$$\\frac{dI_{s,0}}{dt} = \\frac{dI_{s,0}}{dt} - \\tau_{ntsd0} I_{s,0}$$"
            )
          ), ),
          div(class = "col-sm-6", tags$p("For k = 1,2"), tags$p(
            HTML("$$\\frac{dI_{s,k}}{dt} = \\frac{dI_{s,k}}{dt} - I_{s,k}$$")
          ), ),
        ),
        div(
          class = "row eq-box",
          div(
            class = "col-sm-12",
            tags$h4("Symptomatic (Resistant), k = 0,1,2"),
            tags$p(
              HTML(
                "$$\\frac{dI_{r,k}}{dt} = (\\lambda_{ir,k} S + \\lambda_{ar} S)p_{sym_r}p_{ir,k} - \\gamma_{ir,k} I_{r,k} - \\mu_o I_{r,k}$$"
              )
            ),
          ),
          div(class = "col-sm-6", tags$p("For k = 0"), tags$p(
            HTML(
              "$$\\frac{dI_{r,0}}{dt} = \\frac{dI_{r,0}}{dt} - \\tau_{ntrd0} I_{r,0}$$"
            )
          ), ),
          div(class = "col-sm-6", tags$p("For k = 1,2"), tags$p(
            HTML("$$\\frac{dI_{r,k}}{dt} = \\frac{dI_{r,k}}{dt} - I_{r,k}$$")
          ), ),
        ),
        div(
          class = "row eq-box",
        tags$h4("Treatment Success (Sensitive), k = 1,2"),
        tags$p(
          HTML(
            "$$\\frac{dST_{is,k}}{dt} = (1 - Failrate_s) I_{s,k} - (\\mu_o + \\gamma_{infs} + \\tau_{nfs}) ST_{is,k}$$"
          )
        )
        ),
        div(
          class = "row eq-box",
        tags$h4("Treatment Failure (Sensitive), k = 1,2"),
        tags$p(
          HTML(
            "$$\\frac{dF_{is,k}}{dt} = Failrate_s I_{s,k} - (\\mu_o + \\gamma_{ifs} + \\tau_{fs}) F_{is,k}$$"
          )
        )),
        div(
          class = "row eq-box",
        tags$h4("Treatment Success (Resistant)"),
        tags$p(
          HTML(
            "$$\\frac{dST_{ir,k}}{dt} = (1 - Failrate_r) I_{r,k} - (\\mu_o + \\gamma_{infr} + \\tau_{nfr}) ST_{ir,k}$$"
          )
        )),
        div(
          class = "row eq-box",
        tags$h4("Treatment Failure (Resistant)"),
        tags$p(
          HTML(
            "$$\\frac{dF_{ir,k}}{dt} = Failrate_r I_{r,k} - (\\mu_o + \\gamma_{ifr} + \\tau_{fr}) F_{ir,k}$$"
          )
        )),
        div(
          class = "row eq-box",
          
          tags$h4("Gametocytes from Symptomatic (Sensitive)"),
          div(class = "col-md-6", tags$p("For k = 0"), tags$p(
            HTML(
              "$$\\frac{dG_{is,0}}{dt} = I_{s,0}\\gamma_{is,0} - (\\mu_o + \\tau_{is,0}) G_{is,0}$$"
            )
          ), ),
          div(class = "col-md-6", tags$p("For k = 1,2"), tags$p(
            HTML(
              "$$\\frac{dG_{is,k}}{dt} = I_{s,k}\\gamma_{is,k} - (\\mu_o + \\tau_{is,k}) G_{is,k} + ST_{is,k}\\tau_{nfs} + F_{is,k}\\tau_{fs}$$"
            )
          ), ),
        ),
        div(
          class = "row eq-box",
          tags$h4("Gametocytes from Symptomatic (Resistant)"),
          div(class = "col-md-6", tags$p("For k = 0"), tags$p(
            HTML(
              "$$\\frac{dG_{ir,0}}{dt} = I_{r,0}\\gamma_{ir,0} - (\\mu_o + \\tau_{ir,0}) G_{ir,0}$$"
            )
          ), ),
          div(class = "col-md-6", tags$p("For k = 1,2"), tags$p(
            HTML(
              "$$\\frac{dG_{ir,k}}{dt} = I_{r,k}\\gamma_{ir,k} - (\\mu_o + \\tau_{ir,k}) G_{ir,k} + ST_{ir,k}\\tau_{nfr} + F_{ir,k}\\tau_{fr}$$"
            )
          ), ),
        ),
        div(
          class = "row eq-box",
        tags$h4("Recovered (Sensitive)"),
        tags$p(
          HTML(
            "$$\\frac{dR_s}{dt} = \\tau_{as}A_s + \\tau_{ags}GA_s + \\tau_{ntsd0}I_{s,0} + \\sum_{k=0}^{2} G_{is,k}\\tau_{is,k} + \\sum_{k=1}^{2}(ST_{is,k}\\tau_{nfs} + F_{is,k}\\tau_{fs}) - (\\mu_o + \\alpha)R_s$$"
          )
        )),
        div(
          class = "row eq-box",
        tags$h4("Recovered (Resistant)"),
        tags$p(
          HTML(
            "$$\\frac{dR_r}{dt} = \\tau_{ar}A_r + \\tau_{agr}GA_r + \\tau_{ntrd0}I_{r,0} + \\sum_{k=0}^{2} G_{ir,k}\\tau_{ir,k} + \\sum_{k=1}^{2}(ST_{ir,k}\\tau_{nfr} + F_{ir,k}\\tau_{fr}) - (\\mu_o + \\alpha)R_r$$"
          )
        )
      )),
      tags$div(
        class = "general-card",
        tags$h3("Model Parameters"),
        
        tags$table(
          class = "param-table",
          
          tags$thead(tags$tr(
            tags$th("Symbol"), tags$th("Description"), tags$th("Units")
          )),
          
          tags$tbody(
            tags$tr(
              tags$td(tags$b("μ", tags$sub("i")), ", ", tags$b("μ", tags$sub("o"))),
              tags$td("Birth and death rates"),
              tags$td("per month")
            ),
            
            tags$tr(
              tags$td(tags$b("α")),
              tags$td("Rate of waning immunity"),
              tags$td("per month")
            ),
            
            tags$tr(
              tags$td(tags$b("β", tags$sub("si")), ", ", tags$b("β", tags$sub("ri"))),
              tags$td("Transmission coefficients from gametocyte carriers"),
              tags$td("per month")
            ),
            
            tags$tr(
              tags$td(tags$b("β", tags$sub("as")), ", ", tags$b("β", tags$sub("ar"))),
              tags$td("Transmission from asymptomatic"),
              tags$td("per month")
            ),
            tags$tr(
              tags$td(
                tags$b("c", tags$sub("βr"))
              ),
              tags$td("Transmission fitness cost of resistant parasites"),
              tags$td("Unitless")
            ),
            tags$tr(
              tags$td(tags$b("γ", tags$sub("*"))),
              tags$td("Gametocyte production rates"),
              tags$td("per month")
            ),
            
            tags$tr(
              tags$td(tags$b("τ", tags$sub("*"))),
              tags$td("Recovery/removal rates"),
              tags$td("per month")
            ),
            
            tags$tr(
              tags$td(tags$b("p")),
              tags$td("Fraction of infections that are sensitive"),
              tags$td("Unitless")
            ),
            
            tags$tr(
              tags$td(tags$b("p", tags$sub("s")), ", ", tags$b("p", tags$sub("r"))),
              tags$td("Probability of becoming symptomatic"),
              tags$td("Unitless")
            ),
            
            tags$tr(
              tags$td(tags$b("f", tags$sub("s")), ", ", tags$b("f", tags$sub("r"))),
              tags$td("Probability of treatment failure"),
              tags$td("Unitless")
            ),
            
            tags$tr(
              tags$td(

                  tags$b("flo"), ", ",
                  tags$b("amp"), ", ",
                  tags$b("period"), ", ",
                  tags$b("phase")
                
              ),
              tags$td("Seasonal transmission parameters"),
              tags$td("Varies")
            )
            
          )
        )
      ),
      
    ),
    div(
      id = "simulation",
      tabsetPanel(
        id = "main",
        
        tabPanel(
          title = "Total Incidence & Total Population",
          withSpinner(plotOutput("distPlot")),
          br(),
          hr(),
          br(),
          withSpinner(plotOutput("Plot_Population"))
        ),
        tabPanel(
          title = "Fraction of transmission & Resistant fraction",
          withSpinner(plotOutput("Plot_ft")),
          br(),
          hr(),
          br(),
          withSpinner(plotOutput("Plot_p_res")),
        ),
        tabPanel(
          title = "Symptomatic/Asymptomatic Incidence",
          withSpinner(plotOutput("Plot_sym_asym")),
          br(),
          hr(),
          br(),
          withSpinner(plotOutput("Plot_sym_asym_ratio")),
        ),
        tabPanel(
          title = "Sensitive/Resistant Incidence",
          withSpinner(plotOutput("Plot_s_r")),
          br(),
          hr(),
          br(),
          withSpinner(plotOutput("Plot_s_r_ratio"))
        ),
        tabPanel(
          title = "Gametocyte Infections (Symptomatic vs Asymptomatic)",
          withSpinner(plotOutput("Plot_Gasym_sym")),
          br(),
          hr(),
          br(),
          withSpinner(plotOutput("Plot_Gasym_sym_ratio"))
        ),
        tabPanel(
          title = "Gametocyte Infections (Sensitive vs Resistant)",
          withSpinner(plotOutput("Plot_GS_GR")),
          br(),
          hr(),
          br(),
          withSpinner(plotOutput("Plot_GS_GR_ratio"))
        ),
        tabPanel(
          title = "Treated/untreated channel contribution",
          withSpinner(plotOutput("Plot_CT_total")),
          br(),
          hr(),
          br(),
          withSpinner(plotOutput("Plot_CU_total"))
        )
      ),
    ),
    div(
      id = "about_us",
      class = "source-card",
      tags$b("Model Development Team:"),
      br(),
      "Assoc. Prof. Wirichada Pan-ngum (email: ",
      tags$a(href = "mailto:pan@tropmedres.ac", "pan@tropmedres.ac"),
      ")",
      br(),
      "Dr. Sompob Saralamba (email: ",
      tags$a(href = "mailto:Sompob@tropmedres.ac", "Sompob@tropmedres.ac"),
      ")",
      br(),
      "Mr. Tanaphum Wichaita (email: ",
      tags$a(href = "mailto:Tanaphum@tropmedres.ac", "Tanaphum@tropmedres.ac"),
      ")",
      br(),
      br(),
      tags$b("Scientific and Clinical Advisors:"),
      br(),
      "Prof. Arjen Dondorp (email: ",
      tags$a(href = "mailto:arjen@tropmedres.ac", "arjen@tropmedres.ac"),
      ")",
      br(),
      "Dr. Thomas J Peto (email: ",
      tags$a(href = "mailto:Tom@tropmedres.ac", "Tom@tropmedres.ac"),
      ")",
    ),
    div(
      id = "Data",
      tags$div(
        class = "source-card",
        
        tags$h3("Source Data", class = "source-title"),
        tags$p("Below is the list of source data and their citations:", class = "source-subtitle"),
        
        tags$table(
          class = "source-table",
          
          tags$thead(tags$tr(
            tags$th("Description"), tags$th("Citation")
          )),
          
          tags$tbody(
            tags$tr(tags$td("Rate of loss of immunity"), tags$td(
              tags$a(
                href = "https://doi.org/10.1371/journal.pone.0001767",
                "DOI: 10.1371/journal.pone.0001767",
                target = "_blank"
              )
            )),
            
            tags$tr(
              tags$td("Rate of recovery ACT & Primaquine 0.0625 mg/kg"),
              tags$td(
                tags$a(href = "https://doi.org/10.1002/cpt.2512", "DOI: 10.1002/cpt.2512", target = "_blank")
              )
            ),
            
            tags$tr(tags$td("Population of Tanzania"), tags$td(
              tags$a(href = "https://data.who.int/countries/834", "WHO Data Tanzania", target = "_blank")
            )),
            
            tags$tr(tags$td("Malaria Incidence of Tanzania"), tags$td(
              tags$a(
                href = "https://www.who.int/teams/global-malaria-programme/reports/world-malaria-report-2024",
                "WHO World Malaria Report 2024",
                target = "_blank"
              )
            ))
            
          )
        )
      )
    ),
    tags$div(id = "goTopButton", "Go to Top"),
  )
)

server <- function(session, input, output) {
  useShinyjs()
  values <- reactiveValues(
    out_in_year = NULL,
    out_in_year_baseline = NULL,
    init_state = init_state,
    parameters = parameters_Tanzania,
    times = times,
    events = events,
    points_time = 2000:2022,
    obs = tanzania_Incidence[6:14, 4]
  )
  
  
  observeEvent(input$country, {
    if (input$country == "Tanzania") {
      values$init_state <- init_state_tanzania
      values$parameters <- parameters_Tanzania
      updateSliderInput(session, "start_res_year", value = 2000)
      
      year_r <- input$start_res_year - 1999
      values$points_time <- 2015:2023
      values$obs <- tanzania_Incidence[6:14, 4]
      year_d2 <- input$D2_start - 2014
      values$parameters$start_d <- 12 * year_d2
      
    }
    
    values$parameters$start_res_year <- 12 * year_r
    
    values$events <- list(
      func = function(time, state, parameters) {
        if (time == values$parameters$start_res_year) {
          state["Ir0"] <- state["Ir0"] + 1
          state["Ir1"] <- state["Ir1"] + 1
          state["Ar"]  <- state["Ar"]  + 1
          
          state["Is0"] <- state["Is0"] - 1
          state["Is1"] <- state["Is1"] - 1
          state["As"]  <- state["As"]  - 1
        }
        return(state)
      },
      time = 1
    )
    # Transmissibility
    updateSliderInput(session, "Bss", value = values$parameters$beta_s[1])
    updateSliderInput(session, "Bas", value = values$parameters$beta_as)
    updateSliderInput(session, "Bss2", value = values$parameters$beta_s_2[1])
    updateSliderInput(session, "Bas2", value = values$parameters$beta_as_2)
    updateSliderInput(session, "M_beta_resistant", value = 1.131)
    updateSliderInput(session, "C_beta_resistant", value = parameters$c_beta_r)
    
    # Symptomatic
    updateSliderInput(session, "prop_sym", value = 0.25)
    
    # Resistance
    # updateSliderInput(session, "start_res", value = values$parameters$prob_res*100)
    
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
    updateSliderInput(session, "D0", value = 66.67)
    updateSliderInput(session, "D1", value = 0)
    updateSliderInput(session, "D2", value = 33.33)
    updateSliderInput(session,
                      "fail_rate_s",
                      value = values$parameters$Fail_rate_s * 100)
    updateSliderInput(session,
                      "fail_rate_r",
                      value = values$parameters$Fail_rate_r * 100)
    
    # Visualization
    updateSliderInput(session, "time_x", value = 2000)
  })
  
  
  observeEvent("", {
    show("intro")
    hide("simulation")
    hide("about_us")
    hide("Data")
  }, once = TRUE)
  
  observeEvent(input$goto_intro, {
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
    # Transmissibility
    updateSliderInput(session, "Bss", value = parameters_Tanzania$beta_s[1])
    updateSliderInput(session, "Bas", value = parameters_Tanzania$beta_as)
    updateSliderInput(session, "Bss2", value = parameters_Tanzania$beta_s_2[1])
    updateSliderInput(session, "Bas2", value = parameters_Tanzania$beta_as_2)
    updateSliderInput(session, "M_beta_resistant", value = 1.131)
    updateSliderInput(session, "C_beta_resistant", value = parameters$c_beta_r)
    
    # Symptomatic
    updateSliderInput(session, "prop_sym", value = 0.25)
    
    # Resistance
    updateSliderInput(session, "start_res_year", value = 2000)
    # updateSliderInput(session, "start_res", value = parameters_Tanzania$prob_res*100)
    
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
    updateSliderInput(session, "D0", value = 66.67)
    updateSliderInput(session, "D1", value = 0)
    updateSliderInput(session, "D2", value = 33.33)
    updateSliderInput(session, "fail_rate_s", value = parameters_Tanzania$Fail_rate_s *
                        100)
    updateSliderInput(session, "fail_rate_r", value = parameters_Tanzania$Fail_rate_r *
                        100)
    
    # Visualization
    updateSliderInput(session, "time_x", value = 2000)
  })
  
  last_mod <- reactiveVal(list(id = "D0", value = 0))
  observeEvent(input$D0, {
    last_mod(list(id = "D0", value = input$D0))
  })
  observeEvent(input$D1, {
    last_mod(list(id = "D1", value = input$D1))
  })
  observeEvent(input$D2, {
    last_mod(list(id = "D2", value = input$D2))
  })
  
  observeEvent(c(input$D0, input$D1, input$D2), {
    d0 <- isolate(input$D0)
    d1 <- isolate(input$D1)
    d2 <- isolate(input$D2)
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
    
    values$parameters$propIs_new <- c(input$D0 / 100, input$D1 / 100, input$D2 /
                                        100)
    values$parameters$propIr_new <- c(input$D0 / 100, input$D1 / 100, input$D2 /
                                        100)
  })
  
  observeEvent(input$D2_start, {
    year_d2 <- input$D2_start - 1999
    values$parameters$start_d <- 12 * year_d2
    
    
  })
  
  observeEvent(input$D0_treatment, {
    values$parameters$d0 <- (1 / input$D0_treatment) * 30
  })
  
  observeEvent(input$D1_treatment, {
    values$parameters$d1 <- (1 / input$D1_treatment) * 30
  })
  
  observeEvent(input$D2_treatment, {
    values$parameters$d2 <- (1 / input$D2_treatment) * 30
  })
  
  observeEvent(input$start_res_year, {
    if (input$country == "Tanzania") {
      year_r <- input$start_res_year - 2000
    }
    
    values$parameters$start_res_year <- (12 * year_r) + 1
    
    values$events <- list(
      func = function(time, state, parameters) {
        if (time == values$parameters$start_res_year) {
          state["Is0"] <- state["Is0"] - 1
          state["Is1"] <- state["Is1"] - 1
          state["As"]  <- state["As"]  - 1
          
          state["Ir0"] <- state["Ir0"] + 1
          state["Ir1"] <- state["Ir1"] + 1
          state["Ar"]  <- state["Ar"]  + 1
          
          
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
        if (time == 1) {
          state["Is0"] <- state["Is0"] - 1
          state["Is1"] <- state["Is1"] - 1
          state["As"]  <- state["As"]  - 1
          
          state["Ir0"] <- state["Ir0"] + 1
          state["Ir1"] <- state["Ir1"] + 1
          state["Ar"]  <- state["Ar"]  + 1
          
          
        }
        return(state)
      },
      time = 1
    )
  })
  observeEvent(input$fail_rate_s, {
    values$parameters$Fail_rate_s <- input$fail_rate_s / 100
  })
  observeEvent(input$fail_rate_r, {
    values$parameters$Fail_rate_r <- input$fail_rate_r / 100
  })
  
  observeEvent(input$Bss, {
    values$parameters$beta_s <- c(input$Bss, input$Bss, input$Bss)
  })
  observeEvent(input$Bas, {
    values$parameters$beta_as <- input$Bas
  })
  observeEvent(input$Bss2, {
    values$parameters$beta_s_2 <- c(input$Bss2, input$Bss2, input$Bss2)
  })
  observeEvent(input$Bas2, {
    values$parameters$beta_as_2 <- input$Bas2
  })
  observeEvent(input$M_beta_resistant, {
    values$parameters$beta_r <- c(input$Bss, input$Bss, input$Bss) * input$M_beta_resistant
    values$parameters$beta_ar <- input$Bas * input$M_beta_resistant
    values$parameters$beta_r_2 <- c(input$Bss2, input$Bss2, input$Bss2) *
      input$M_beta_resistant
    values$parameters$beta_ar_2 <- input$Bas2 * input$M_beta_resistant
  })
  observeEvent(input$C_beta_resistant, {
    values$parameters$c_beta_r <- input$C_beta_resistant
  })
  
  # Data for Source Data tab
  source_data <- data.frame(
    Description = c(
      "Rate of loss of immunity",
      "Rate of recovery ACT & Primaquine 0.0625 mg/kg",
      "Population of Tanzania",
      "Malaria Incidence of Tanzania"
    ),
    Citation = c(
      "DOI: 10.1371/journal.pone.0001767",
      "DOI: 10.1002/cpt.2512",
      "https://statisticstimes.com/demographics/country/tanzania-population.php",
      "https://www.who.int/teams/global-malaria-programme/reports/world-malaria-report-2024"
    ),
    stringsAsFactors = FALSE
  )
  
  output$source_table <- renderTable({
    source_data
  })
  
  observeEvent(input$g_is, {
    values$parameters$g_is <- c(input$g_is, input$g_is, input$g_is)
  })
  observeEvent(input$g_ir, {
    values$parameters$g_ir <- c(input$g_ir, input$g_ir, input$g_ir)
  })
  observeEvent(input$g_as, {
    values$parameters$g_as <- input$g_as
  })
  observeEvent(input$g_ar, {
    values$parameters$g_ar <- input$g_ar
  })
  observeEvent(input$g_infs, {
    values$parameters$g_infs <- c(input$g_infs, input$g_infs)
  })
  observeEvent(input$g_infr, {
    values$parameters$g_infr <- c(input$g_infr, input$g_infr)
  })
  observeEvent(input$g_ifr, {
    values$parameters$g_ifr <- c(input$g_ifr, input$g_ifr)
  })
  observeEvent(input$g_ifs, {
    values$parameters$g_ifs <- c(input$g_ifs, input$g_ifs)
  })
  
  
  observeEvent(input$prop_sym, {
    values$parameters$prob_sym_s <- input$prop_sym
    values$parameters$prob_sym_r <- input$prop_sym
  })
  
  ode_results <- reactive({
    req(input$D0 + input$D1 + input$D2 == 100 &
          values$parameters$start_res_year >= 0)
    
    out <- ode(
      y = values$init_state,
      times = values$times,
      func = Malaria_model_with_Array,
      parms = values$parameters,
      events = values$events
    )
    vars_all <- colnames(out)
    summarise_by_year(out, vars_all, 36, start_year = 2000)
    
  })
  
  ode_baseline <- reactive({
    req(input$D0 + input$D1 + input$D2 == 100 &
          values$parameters$start_res_year >= 0)
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
    summarise_by_year(out, vars_all, 36, start_year = 2000)
    
    
  })
  
  output$distPlot <- renderPlot({
    out_in_year <- ode_results()
    out_in_year_baseline <- ode_baseline()
    max_y <- max(out_in_year[1:20, "inc"], out_in_year_baseline[, "inc"])
    
    plot(
      out_in_year[, 1],
      out_in_year[, "inc"],
      type = "l",
      xlab = "Years",
      ylab = "Incidence",
      main = "Total Incidence by Year",
      xlim = c(input$time_x, 2035),
      lwd = 2,
      lty = 2,
      col = 1,
      ylim = c(0, 1.2 * max_y)
    )
    points(values$points_time,
           values$obs,
           col = "red",
           pch = 19)
    lines(out_in_year_baseline[, 1], out_in_year_baseline[, "inc"], lwd = 2)
    abline(v = input$D2_start, lty = 2)
    legend(
      "topleft",
      c("Model", "Data"),
      col = c(1, 2),
      lty = c(1, NA),
      pch = c(NA, 19)
    )
    legend(
      "bottomright",
      c("Baseline", "Treatment D2"),
      col = c(1),
      lty = c(1, 2),
      lwd = 2
    )
    text(input$D2_start - 1, 2e07, "Start D2")
  })
  
  output$Plot_Population <- renderPlot({
    out_in_year <- ode_results()
    out_in_year_baseline <- ode_baseline()
    max_y <- max(out_in_year[1:20, "N"], out_in_year_baseline[, "N"])
    
    plot(
      out_in_year[, 1],
      out_in_year[, "N"],
      type = "l",
      xlab = "Years",
      ylab = "Population",
      main = "Total Population by Year",
      xlim = c(input$time_x, 2035),
      lwd = 2,
      lty = 1,
      col = "green2",
      ylim = c(0, 1.2 * max_y)
    )
    points(2000:2050,
           tanzania_Population,
           col = 1,
           pch = 19)
    legend(
      "right",
      c("Model", "Data"),
      col = c("green2", 1),
      lty = c(1, NA),
      pch = c(NA, 19)
    )
  })
  
  output$Plot_sym_asym <- renderPlot({
    out_in_year <- ode_results()
    out_in_year_baseline <- ode_baseline()
    max_y <- max(c(
      out_in_year[1:35, "inc_sym"],
      out_in_year[1:35, "inc_asym"],
      out_in_year_baseline[1:35, "inc_sym"],
      out_in_year_baseline[1:35, "inc_asym"]
    ))
    
    plot(
      out_in_year[, 1],
      out_in_year[, "inc_sym"],
      type = "l",
      xlab = "Years",
      ylab = "Incidence",
      main = "Symptomatic/Asymptomatic Incidence by Year",
      xlim = c(input$time_x, 2035),
      lwd = 2,
      lty = 2,
      col = "blue2",
      ylim = c(0, 1.2 * max_y)
    )
    lines(
      out_in_year[, 1],
      out_in_year[, "inc_asym"],
      col = "green4",
      lty = 2,
      lwd = 2
    )
    lines(
      out_in_year_baseline[, 1],
      out_in_year_baseline[, "inc_sym"],
      col = "blue2",
      lwd = 2
    )
    lines(
      out_in_year_baseline[, 1],
      out_in_year_baseline[, "inc_asym"],
      col = "green4",
      lwd = 2
    )
    legend(
      "bottomright",
      c("Symptomatic", "Asymptomatic"),
      col = c("blue2", "green4"),
      lty = c(1, 1),
      lwd = 2
    )
    abline(v = input$D2_start, lty = 2)
    text(input$D2_start - 1, 2e07, "Start D2")
    legend(
      "topleft",
      c("Baseline", "Treatment D2"),
      col = c(1),
      lty = c(1, 2),
      lwd = 2
    )
  })
  
  output$Plot_sym_asym_ratio <- renderPlot({
    out_in_year <- ode_results()
    out_in_year_baseline <- ode_baseline()
    
    plot(
      out_in_year[, 1],
      out_in_year[, "inc_sym"] / out_in_year[, "inc"],
      type = "l",
      xlab = "Years",
      ylab = "Proportion of incidence type",
      main = "Proportion of Symptomatic vs Asymptomatic Incidence Over Time",
      xlim = c(input$time_x, 2035),
      lwd = 2,
      lty = 2,
      col = "blue2",
      ylim = c(0, 1)
    )
    lines(
      out_in_year[, 1],
      out_in_year[, "inc_asym"] / out_in_year[, "inc"],
      col = "green4",
      lty = 2,
      lwd = 2
    )
    lines(
      out_in_year_baseline[, 1],
      out_in_year_baseline[, "inc_sym"] / out_in_year_baseline[, "inc"],
      col = "blue2",
      lwd = 2
    )
    lines(
      out_in_year_baseline[, 1],
      out_in_year_baseline[, "inc_asym"] / out_in_year_baseline[, "inc"],
      col = "green4",
      lwd = 2
    )
    abline(v = input$D2_start, lty = 2)
    legend(
      "bottomright",
      c("Proportion Symptomatic", "Proportion Asymptomatic"),
      col = c("blue2", "green4"),
      lty = c(1, 1),
      lwd = 2
    )
    text(input$D2_start - 1, 0.8, "Start D2")
    legend(
      "topright",
      c("Baseline", "Treatment D2"),
      col = c(1),
      lty = c(1, 2),
      lwd = 2
    )
  })
  
  output$Plot_s_r <- renderPlot({
    out_in_year <- ode_results()
    out_in_year_baseline <- ode_baseline()
    max_y <- max(c(
      out_in_year[1:35, "inc_s"],
      out_in_year[1:35, "inc_r"],
      out_in_year_baseline[1:35, "inc_s"],
      out_in_year_baseline[1:35, "inc_r"]
    ))
    
    plot(
      out_in_year[, 1],
      out_in_year[, "inc_s"],
      type = "l",
      xlab = "Years",
      ylab = "Incidence",
      main = "Sensitive/Resistant Incidence by Year",
      xlim = c(input$time_x, 2035),
      lwd = 2,
      lty = 2,
      col = "black",
      ylim = c(0, 1.2 * max_y)
    )
    lines(
      out_in_year[, 1],
      out_in_year[, "inc_r"],
      col = "red",
      lty = 2,
      lwd = 2
    )
    lines(
      out_in_year_baseline[, 1],
      out_in_year_baseline[, "inc_r"],
      col = "red",
      lwd = 2
    )
    lines(out_in_year_baseline[, 1], out_in_year_baseline[, "inc_s"], lwd = 2)
    abline(v = input$D2_start, lty = 2)
    legend(
      "bottomright",
      c("Sensitive", "Resistant"),
      col = c("black", "red"),
      lty = c(1, 1),
      lwd = 2
    )
    text(input$D2_start - 1, 2e07, "Start D2")
    legend(
      "topleft",
      c("Baseline", "Treatment D2"),
      col = c(1),
      lty = c(1, 2),
      lwd = 2
    )
  })
  
  output$Plot_s_r_ratio <- renderPlot({
    out_in_year <- ode_results()
    out_in_year_baseline <- ode_baseline()
    plot(
      out_in_year[, 1],
      out_in_year[, "inc_s"] / out_in_year[, "inc"],
      type = "l",
      xlab = "Years",
      ylab = "Proportion of parasite type",
      main = "Proportion of Sensitive vs Resistant Parasites Over Time",
      xlim = c(input$time_x, 2035),
      lwd = 2,
      lty = 2,
      col = "black",
      ylim = c(0, 1)
    )
    lines(
      out_in_year[, 1],
      out_in_year[, "inc_r"] / out_in_year[, "inc"],
      col = "red",
      lty = 2,
      lwd = 2
    )
    lines(
      out_in_year_baseline[, 1],
      out_in_year_baseline[, "inc_r"] / out_in_year_baseline[, "inc"],
      col = "red",
      lwd = 2
    )
    lines(out_in_year_baseline[, 1],
          out_in_year_baseline[, "inc_s"] / out_in_year_baseline[, "inc"],
          lwd = 2)
    abline(v = input$D2_start, lty = 2)
    legend(
      "left",
      c("Proportion Sensitive", "Proportion Resistant"),
      col = c("black", "red"),
      lty = c(1, 1),
      lwd = 2
    )
    text(input$D2_start - 1, 0.8, "Start D2")
    legend(
      "right",
      c("Baseline", "Treatment D2"),
      col = c(1),
      lty = c(1, 2),
      lwd = 2
    )
  })
  
  output$Plot_Gasym_sym <- renderPlot({
    out_in_year <- ode_results()
    out_in_year_baseline <- ode_baseline()
    max_y <- max(c(
      out_in_year[1:35, "Gasym_inc"],
      out_in_year[1:35, "Gsym_inc"],
      out_in_year_baseline[1:35, "Gasym_inc"],
      out_in_year_baseline[1:35, "Gsym_inc"]
    ))
    
    plot(
      out_in_year[, 1],
      out_in_year[, "Gasym_inc"],
      type = "l",
      xlab = "Years",
      ylab = "Incidence",
      lty = 2,
      main = "Gametocyte Infections Asymptomatic/Symptomatic by Year",
      xlim = c(input$time_x, 2035),
      lwd = 2,
      col = 3,
      ylim = c(0, 1.2 * max(
        out_in_year_baseline[1:35, "Gasym_inc"], out_in_year_baseline[1:35, "Gsym_inc"]
      ))
    )
    lines(
      out_in_year_baseline[, 1],
      out_in_year_baseline[, "Gasym_inc"],
      lty = 1,
      lwd = 2,
      col = 3
    )
    lines(
      out_in_year[, 1],
      out_in_year[, "Gsym_inc"],
      lwd = 2,
      lty = 2,
      col = "blue"
    )
    lines(
      out_in_year_baseline[, 1],
      out_in_year_baseline[, "Gsym_inc"],
      lwd = 2,
      col = "blue"
    )
    abline(v = input$D2_start, lty = 2)
    legend(
      "topleft",
      c("Asymptomatic", "Symptomatic"),
      col = c(3, 4),
      lty = c(1, 1)
    )
    text(input$D2_start - 1, 4e05, "Start D2")
    legend(
      "top",
      c("Baseline", "Treatment D2"),
      col = c(1),
      lty = c(1, 2),
      lwd = 2
    )
  })
  
  output$Plot_Gasym_sym_ratio <- renderPlot({
    out_in_year <- ode_results()
    out_in_year_baseline <- ode_baseline()
    plot(
      out_in_year[, 1],
      out_in_year[, "Gasym_inc"] / out_in_year[, "G_inc"],
      type = "l",
      xlab = "Years",
      ylab = "Proportion",
      lty = 2,
      main = "Asymptomatic/Symptomatic Gametocyte Infections Proportion Over Time",
      xlim = c(input$time_x, 2035),
      lwd = 2,
      col = 3,
      ylim = c(0, 1)
    )
    lines(
      out_in_year_baseline[, 1],
      out_in_year_baseline[, "Gasym_inc"] / out_in_year_baseline[, "G_inc"],
      lty = 1,
      lwd = 2,
      col = 3
    )
    lines(
      out_in_year[, 1],
      out_in_year[, "Gsym_inc"] / out_in_year[, "G_inc"],
      lwd = 2,
      lty = 2,
      col = "blue"
    )
    lines(
      out_in_year_baseline[, 1],
      out_in_year_baseline[, "Gsym_inc"] / out_in_year_baseline[, "G_inc"],
      lwd = 2,
      col = "blue"
    )
    abline(v = input$D2_start, lty = 2)
    legend(
      "left",
      c("Proportion Asymptomatic", "Proportion Symptomatic"),
      col = c(3, 4),
      lty = c(1, 1)
    )
    text(input$D2_start - 1, 0.8, "Start D2")
    legend(
      "right",
      c("Baseline", "Treatment D2"),
      col = c(1),
      lty = c(1, 2),
      lwd = 2
    )
  })
  
  output$Plot_GS_GR <- renderPlot({
    out_in_year <- ode_results()
    out_in_year_baseline <- ode_baseline()
    max_y <- max(c(
      out_in_year[1:35, "GR_inc"],
      out_in_year[1:35, "GR_inc"],
      out_in_year_baseline[1:35, "GR_inc"],
      out_in_year_baseline[1:35, "GR_inc"]
    ))
    
    plot(
      out_in_year[, 1],
      out_in_year[, "GS_inc"],
      type = "l",
      xlab = "Years",
      ylab = "Incidence",
      main = "Gametocyte Infections Sensitive/Resistant by Year",
      xlim = c(input$time_x, 2035),
      lwd = 2,
      col = 1,
      lty = 2,
      ylim = c(0, 1.2 * max_y)
    )
    lines(out_in_year_baseline[, 1],
          out_in_year_baseline[, "GS_inc"],
          lty = 1,
          lwd = 2)
    lines(
      out_in_year[, 1],
      out_in_year[, "GR_inc"],
      lwd = 2,
      lty = 2,
      col = "red"
    )
    lines(
      out_in_year_baseline[, 1],
      out_in_year_baseline[, "GR_inc"],
      lwd = 2,
      col = "red"
    )
    abline(v = input$D2_start, lty = 2)
    legend(
      "topleft",
      c("Sensitive", "Resistant"),
      col = c(1, 2),
      lty = c(1, 1)
    )
    text(input$D2_start - 1, 4e05, "Start D2")
    legend(
      "top",
      c("Baseline", "Treatment D2"),
      col = c(1),
      lty = c(1, 2),
      lwd = 2
    )
  })
  
  output$Plot_GS_GR_ratio <- renderPlot({
    out_in_year <- ode_results()
    out_in_year_baseline <- ode_baseline()
    plot(
      out_in_year[, 1],
      out_in_year[, "GS_inc"] / out_in_year[, "G_inc"],
      type = "l",
      xlab = "Years",
      ylab = "Proportion",
      main = "Sensitive/Resistant Gametocyte Infections Proportion Over Time",
      xlim = c(input$time_x, 2035),
      lwd = 2,
      col = 1,
      lty = 2,
      ylim = c(0, 1)
    )
    lines(
      out_in_year_baseline[, 1],
      out_in_year_baseline[, "GS_inc"] / out_in_year_baseline[, "G_inc"],
      lty = 1,
      lwd = 2
    )
    lines(
      out_in_year[, 1],
      out_in_year[, "GR_inc"] / out_in_year[, "G_inc"],
      lwd = 2,
      lty = 2,
      col = "red"
    )
    lines(
      out_in_year_baseline[, 1],
      out_in_year_baseline[, "GR_inc"] / out_in_year_baseline[, "G_inc"],
      lwd = 2,
      col = "red"
    )
    abline(v = input$D2_start, lty = 2)
    legend(
      "left",
      c("Proportion Sensitive", "Proportion Resistant"),
      col = c(1, 2),
      lty = c(1, 1)
    )
    text(input$D2_start - 1, 0.8, "Start D2")
    legend(
      "right",
      c("Baseline", "Treatment D2"),
      col = c(1),
      lty = c(1, 2),
      lwd = 2
    )
  })
  
  output$Plot_ft <- renderPlot({
    out_in_year <- ode_results()
    out_in_year_baseline <- ode_baseline()
    max_y <- max(c(out_in_year[1:35, "F_T"]))
    
    out_in_year <- ode_results()
    out_in_year_baseline <- ode_baseline()
    plot(
      out_in_year[, 1],
      out_in_year[, "F_T"],
      type = "l",
      xlab = "Years",
      ylab = "Proportion",
      main = "FT",
      xlim = c(input$time_x, 2035),
      lwd = 2,
      col = 1,
      lty = 2,
      ylim = c(0, max_y * 1.2)
    )
    lines(out_in_year_baseline[, 1],
          out_in_year_baseline[, "F_T"],
          lty = 1,
          lwd = 2)
    abline(v = input$D2_start, lty = 2)
    text(input$D2_start - 1, 0.8, "Start D2")
    legend(
      "left",
      c("Baseline", "Treatment D2"),
      col = c(1),
      lty = c(1, 2),
      lwd = 2
    )
  })
  
  output$Plot_ft <- renderPlot({
    out_in_year <- ode_results()
    out_in_year_baseline <- ode_baseline()
    max_y <- max(c(out_in_year[1:35, "F_T"]))
    
    out_in_year <- ode_results()
    out_in_year_baseline <- ode_baseline()
    plot(
      out_in_year[, 1],
      out_in_year[, "F_T"],
      type = "l",
      xlab = "Years",
      ylab = "FT",
      main = "Fraction of transmission",
      xlim = c(input$time_x, 2035),
      lwd = 2,
      col = 1,
      lty = 2,
      ylim = c(0, max_y * 1.2)
    )
    lines(out_in_year_baseline[, 1],
          out_in_year_baseline[, "F_T"],
          lty = 1,
          lwd = 2)
    abline(v = input$D2_start, lty = 2)
    text(input$D2_start - 1, 0.8, "Start D2")
    legend(
      "left",
      c("Baseline", "Treatment D2"),
      col = c(1),
      lty = c(1, 2),
      lwd = 2
    )
  })
  
  output$Plot_p_res <- renderPlot({
    out_in_year <- ode_results()
    out_in_year_baseline <- ode_baseline()
    max_y <- max(c(out_in_year[1:35, "p_resistant_inc"]))
    
    
    out_in_year <- ode_results()
    out_in_year_baseline <- ode_baseline()
    plot(
      out_in_year[, 1],
      out_in_year[, "p_resistant_inc"],
      type = "l",
      xlab = "Years",
      ylab = "Resistant Infection Fraction",
      main = "Resistant fraction of incident infections",
      xlim = c(input$time_x, 2035),
      lwd = 2,
      col = 1,
      lty = 2,
      ylim = c(0, 1)
    )
    lines(out_in_year_baseline[, 1],
          out_in_year_baseline[, "p_resistant_inc"],
          lty = 1,
          lwd = 2)
    abline(v = input$D2_start, lty = 2)
    text(input$D2_start - 1, 0.8, "Start D2")
    legend(
      "left",
      c("Baseline", "Treatment D2"),
      col = c(1),
      lty = c(1, 2),
      lwd = 2
    )
  })
  
  output$Plot_CT_total <- renderPlot({
    out_in_year <- ode_results()
    out_in_year_baseline <- ode_baseline()
    max_y <- max(c(out_in_year[1:35, "CT_total"], out_in_year_baseline[, "CT_total"]))
    
    out_in_year <- ode_results()
    out_in_year_baseline <- ode_baseline()
    plot(
      out_in_year[, 1],
      out_in_year[, "CT_total"],
      type = "l",
      xlab = "Years",
      ylab = "Treated Contribution",
      main = "Treated Channel Contribution",
      xlim = c(input$time_x, 2035),
      lwd = 2,
      col = 1,
      lty = 2,
      ylim = c(0, max_y * 1.2)
    )
    lines(out_in_year_baseline[, 1],
          out_in_year_baseline[, "CT_total"],
          lty = 1,
          lwd = 2)
    abline(v = input$D2_start, lty = 2)
    text(input$D2_start - 1, 0.8, "Start D2")
    legend(
      "bottomleft",
      c("Baseline", "Treatment D2"),
      col = c(1),
      lty = c(1, 2),
      lwd = 2
    )
  })
  
  output$Plot_CU_total <- renderPlot({
    out_in_year <- ode_results()
    out_in_year_baseline <- ode_baseline()
    max_y <- max(c(out_in_year[1:35, "CU_total"], out_in_year_baseline[, "CU_total"]))
    
    out_in_year <- ode_results()
    out_in_year_baseline <- ode_baseline()
    plot(
      out_in_year[, 1],
      out_in_year[, "CU_total"],
      type = "l",
      xlab = "Years",
      ylab = "Untreated Contribution",
      main = "Untreated Channel Contribution",
      xlim = c(input$time_x, 2035),
      lwd = 2,
      col = 1,
      lty = 2,
      ylim = c(0, max_y * 1.2)
    )
    lines(out_in_year_baseline[, 1],
          out_in_year_baseline[, "CU_total"],
          lty = 1,
          lwd = 2)
    abline(v = input$D2_start, lty = 2)
    text(input$D2_start - 1, 0.8, "Start D2")
    legend(
      "bottomleft",
      c("Baseline", "Treatment D2"),
      col = c(1),
      lty = c(1, 2),
      lwd = 2
    )
  })
  
  
}

shinyApp(ui = ui, server = server)
