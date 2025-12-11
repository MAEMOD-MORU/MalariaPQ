# Run Simulation
library(deSolve)
library(manipulate)
#Read model
source("model_D0_33_D1_67_baseline_change_birth_death_2betaSets_SS.R")

source("init_parameter_after_1st_run_2015_D1_33_SS.R")

Tanzania_data <- read.csv("data/Reported malaria cases by method of confirmation.csv")
Tanzania_Incidence <- Tanzania_data[6:14,c(1,4)]
Tanzania_k13_Allele_frequency <- read.csv("data/Tanzania_K13_Allele_frequency.csv")[,2]

init_state <- readRDS("init_state_d1_33.rds")

times <- seq(0,12*81)

# Poputation
totalpop_2010_2100 <- read.csv("data/Tanzania_pop_2010_2100.csv", header = TRUE)

#test model ode
out <- ode(y = init_state, times = times, func = Malaria_model_with_Array, parms = parameters)

#
plot(out[,"time"],out[,"N"], type = "l", xlab = "Months", ylab = "Incidence Sensitive",
     main = "Incidence Sensitive", col = "purple", lwd = 2, ylim = c(0, 1.1*max(out[,"N"])))

events <- list(
  func = function(time, state, parameters) {
    if (time == 12 * 1) {
      state["Ir0"] <- state["Ir0"] + 1
      state["Ir1"] <- state["Ir1"] + 1
      # state["Ir2"] <- state["Ir2"] + 1
      state["Ar"] <- state["Ar"] + 1
      
      # reduce Is
      state["Is0"] <- state["Is0"] -1
      state["Is1"] <- state["Is1"] -1
      # state["Is2"] <- state["Is2"] -1
      state["As"] <- state["As"] -1
      
    }
    return(state)
  },
  time = 12 * 1  # only one time point
)

out_event <- ode(y = init_state, times = times, 
                 func = Malaria_model_with_Array, 
                 parms = parameters, events = events)


# Plot the results inc_s and inc_r 0 : 6000
plot(out_event[,"inc"], type = "l", xlab = "Months", ylab = "Incidence Sensitive",
     main = "Incidence Sensitive", col = "purple", lwd = 2, ylim = c(0, 1.1*max(out[,"inc_s"])))
# lines(out[,"inc_r"], col = "red", lwd = 2)

plot(out_event[,"inc_s"], type = "l", xlab = "Months", ylab = "Incidence Sensitive",
     main = "Incidence Sensitive", col = "purple", lwd = 2, ylim = c(0, 1.1*max(out[,"inc_s"])))
lines(out_event[,"inc_r"], col = "red", lwd = 2)
legend("topleft", legend = c("Incidence Sensitive", "Incidence Resistant"),
       col = c("purple", "red"), lwd = 2, cex = 0.8)

parameters_manipulate <- parameters

time_manipulate <- seq(1, 12*51, by = 1)  # Define the time points for optimized parameters

manipulate(
  {
    parameters_manipulate$beta_s  <- rep(beta_s, 3)   # assuming same for all Is subgroups
    parameters_manipulate$beta_as <- beta_as          # assuming same as beta_s0
    parameters_manipulate$beta_r  <- m * rep(beta_s, 3)
    parameters_manipulate$beta_ar <- m * beta_as
    
    parameters_manipulate$c_beta_r <- c_beta_r 
    
    ## Transmission level from radio (75% or 50%) ----
    if (Tra_level == "High Transmission"){
      parameters_manipulate$prob_sym_s <- 0.25
      parameters_manipulate$prob_sym_r <- 0.25
    }else{
      parameters_manipulate$prob_sym_s <- 0.5
      parameters_manipulate$prob_sym_r <- 0.5
    }
    
    events <- list(
      func = function(time, state, parameters) {
        if (time == (12 * time_start_res_after) + 1) {
          state["Ir0"] <- state["Ir0"] + 1
          state["Ir1"] <- state["Ir1"] + 1
          # state["Ir2"] <- state["Ir2"] + 1
          state["Ar"]  <- state["Ar"]  + 1
          
          # reduce Is
          state["Is0"] <- state["Is0"] - 1
          state["Is1"] <- state["Is1"] - 1
          # state["Is2"] <- state["Is2"] - 1
          state["As"]  <- state["As"]  - 1
        }
        return(state)
      },
      time = (12 * time_start_res_after) + 1  # only one time point
    )
    
    init_state_test <- init_state * i
    out <- ode(
      y     = init_state_test,
      times = time_manipulate,
      func  = Malaria_model_with_Array,
      parms = parameters_manipulate,
      events = events
    )
    
    ## Pre-calc all things you might plot ----
    mat_inc  <- matrix(out[,"inc"],   nrow = 12)
    inc_model_2000_2050   <- colSums(mat_inc)
    
    mat_inc  <- matrix(out[,"inc_r"], nrow = 12)
    inc_model_2000_2050_r <- colSums(mat_inc)
    
    mat_inc  <- matrix(out[,"N"],     nrow = 12)
    inc_model_2000_2050_n <- colMeans(mat_inc)
    
    mat_ratio <- matrix(out[,"inc_r"] / out[,"inc"], nrow = 12)
    ratio_2000_2050 <- colMeans(mat_ratio)
    
    # symptomatic / asymptomatic
    mat_inc  <- matrix(out[,"inc_sym_r"], nrow = 12)
    inc_model_2000_2050_sym_r <- colSums(mat_inc)
    
    mat_inc  <- matrix(out[,"inc_sym"], nrow = 12)
    inc_model_2000_2050_sym   <- colSums(mat_inc)
    
    mat_inc  <- matrix(out[,"inc_asym_r"], nrow = 12)
    inc_model_2000_2050_asym_r <- colSums(mat_inc)
    
    mat_inc  <- matrix(out[,"inc_asym"], nrow = 12)
    inc_model_2000_2050_asym   <- colSums(mat_inc)
    
    years <- 2000:2050
    
    ## Choose which plot(s) to draw based on radio/picker ----
    if (plot_type == "Total Inc + Res") {
      par(mfrow = c(2, 1))
      
      # Top: total incidence
      plot(
        years, inc_model_2000_2050,
        xlab = "Year", ylab = "Incidence",
        lwd  = 2, type = "l",
        main = "Total Incidence",
        xlim = c(year_plot, 2050),
        ylim = c(0, max(inc_model_2000_2050) * 1.1)
      )
      lines(years, inc_model_2000_2050_r, col = "red")
      lines(years, inc_model_2000_2050_n, col = "green3")
      points(2015:2023, Tanzania_Incidence[, 2], col = "blue", lwd = 2)
      
      # Bottom: resistance / total ratio
      plot(
        years, ratio_2000_2050,
        type = "l", col = "red",
        xlab = "Year", ylab = "Ratio",
        xlim = c(year_plot, 2050),
        main = "Ratio Resistant Incidence / Total Incidence"
      )
      # abline(h = ref_val, lty = 2)  # 75% or 50% line
      points(
        c(2016:2022, 2025),
        c(Tanzania_k13_Allele_frequency, 0.4),
        col = 2, lwd = 2, pch = 3
      )
      abline(v = 2025, lty = 2)
      text(2027,0.4,round(ratio_2000_2050[26],2))
      
    } else if (plot_type == "Symptomatic") {
      par(mfrow = c(2, 1))
      
      # Symptomatic R / Total Inc
      plot(
        years,
        inc_model_2000_2050_sym_r / inc_model_2000_2050,
        type = "l", col = "red",
        xlab = "Year", ylab = "Ratio",
        main = "Symptomatic Resistant Incidence / Total Incidence",
        ylim=c(0,0.5)
      )
      points(c(2016:2022,2025),c(Tanzania_k13_Allele_frequency,0.4),col=2,lwd=2,pch=3)
      abline(v = 2025, lty = 2)
      text(2027,0.4,round(inc_model_2000_2050_sym_r[26] / inc_model_2000_2050[26],2),pos=4)
      
      # Symptomatic R / Symptomatic Inc
      plot(
        years,
        inc_model_2000_2050_sym_r / inc_model_2000_2050_sym,
        type = "l", col = "red",
        xlab = "Year", ylab = "Ratio",
        main = "Symptomatic Resistant Incidence / Symptomatic Incidence",
        ylim = c(0, 0.5)
      )
      points(c(2016:2022,2025),c(Tanzania_k13_Allele_frequency,0.4),col=2,lwd=2,pch=3)
      abline(v = 2025, lty = 2)
      text(2027,0.4,round(inc_model_2000_2050_sym_r[26] / inc_model_2000_2050_sym[26],2),pos=4)
      # abline(h = ref_val, lty = 2)
      
    } else if (plot_type == "Asymptomatic") {
      par(mfrow = c(2, 1))
      
      # Asymptomatic R / Total Inc
      plot(
        years,
        inc_model_2000_2050_asym_r / inc_model_2000_2050,
        type = "l", col = "red",
        xlab = "Year", ylab = "Ratio",
        main = "Asymptomatic Resistant Incidence / Total Incidence",
        ylim=c(0,0.5)
      )
      points(c(2016:2022,2025),c(Tanzania_k13_Allele_frequency,0.4),col=2,lwd=2,pch=3)
      abline(v = 2025, lty = 2)
      text(2027,0.4,round(inc_model_2000_2050_asym_r[26] / inc_model_2000_2050[26],2),pos=4)
      # abline(h = ref_val, lty = 2)
      
      # Asymptomatic R / Asymptomatic Inc
      plot(
        years,
        inc_model_2000_2050_asym_r / inc_model_2000_2050_asym,
        type = "l", col = "red",
        xlab = "Year", ylab = "Ratio",
        main = "Asymptomatic Resistant Incidence / Asymptomatic Incidence",
        ylim = c(0, 0.5)
      )
      points(c(2016:2022,2025),c(Tanzania_k13_Allele_frequency,0.4),col=2,lwd=2,pch=3)
      abline(v = 2025, lty = 2)
      text(2027,0.4,round(inc_model_2000_2050_asym_r[26] / inc_model_2000_2050_asym[26],2),pos=4)
      # abline(h = ref_val, lty = 2)
    }
  },
  
  ## 3-type radio for plot group
  plot_type = picker(
    "Total Inc + Res",
    "Symptomatic",
    "Asymptomatic",
    initial = "Total Inc + Res"
  ),
  
  ## radio for 75% or 50% Transmission level
  Tra_level = picker(
    "High Transmission",
    "Moderate Transmission",
    initial = "High Transmission"
  ),
  
  beta_s  = slider(0, 10, step = 0.00001, initial = 0.671),
  beta_as = slider(0, 10, step = 0.00001, initial = 5.992),
  m       = slider(0.1, 3, step = 0.01, initial = 1.05),
  i       = slider(0.00001, 2, step = 0.00001, initial = 0.092),
  c_beta_r = slider(0.01, 2, step = 0.001, initial = 0.745),
  year_plot = slider(2000, 2035, step = 1, initial = 2000),
  time_start_res_after = slider(0, 25, step = 1, initial = 0)
)

