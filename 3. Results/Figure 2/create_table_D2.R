library(deSolve)
library(openxlsx)

source("../../model/model_D0_33_D1_67_baseline_change_birth_death_2betaSets_SS.R")
source("../init_parameter_MCMC.R")

init_state <- readRDS("../init_state_2000.rds")

# --- Settings ---
yrs <- 2000:2050
sim_time <- seq(1, 12 * 51)
asym_seq <- seq(0, 100, by = 10)

start_year_d2 <- 2025
idx1 <- which(yrs == start_year_d2 + 1) # 1-year window end (year+1) (2026)
idx3 <- which(yrs == (start_year_d2 + 3))  # 3-year window end (year+2) (2028)

events <- list(
  func = function(time, state, parameters) {
    if (time == 1) {
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
  time =1  # only one time point
)

pct_red <- function(base, scen) {
  # base should be scalar, scen can be vector
  if (is.na(base) || any(is.na(scen))) return(rep(NA_real_, length(scen)))
  100 * (base - scen) / base
}

# --- Baseline row (D0=0.67, D1=0.33, D2=0) (once) ---
parameters_bl <- parameters
parameters_bl$propIs_new <- c(0.67, 0.33, 0) # no D2
parameters_bl$propIr_new <- c(0.67, 0.33, 0)

out_bl <- ode(y = init_state, times = sim_time,
              func = Malaria_model_with_Array, parms = parameters_bl,
              events = events)
annual_inc_r_base  <- colSums(matrix(out_bl[, "inc_r"],  nrow = 12))
annual_GR_inc_base <- colSums(matrix(out_bl[, "GR_inc"], nrow = 12))

row_baseline <- data.frame(
  Asym                 = "Baseline",
  c_beta               = if (length(parameters$c_beta) == 1) parameters$c_beta else NA_real_,
  D0 = 0.67, D1 = 0.33, D2 = 0,
  Incidence_Res_year_1 = annual_inc_r_base[idx1],
  Incidence_Res_year_3 = annual_inc_r_base[idx3],
  red_inc_r_1y         =  0,
  red_inc_r_3y         =  0,
  Gametocyte_Res_year_1    = annual_GR_inc_base[idx1],
  Gametocyte_Res_year_3    = annual_GR_inc_base[idx3],
  red_GR_inc_1y        = 0,
  red_GR_inc_3y        = 0
)

# --- Workbook ---
wb <- createWorkbook()

# Grid for D2 from 0 to 1
D2_seq <- seq(0, 1, by = 0.1)
pb <- txtProgressBar(min = 0, max = length(D2_seq), style = 3)

# --- Loop over D2 values ---
for (D2 in D2_seq) {
  
  # Choose how to set D0, D1, D2 (must sum to 1)
  D1 <- 0
  D0 <- 1 - D2
  if (D0 < 0) next  # safety, though not needed with 0..1
  
  # parameters for this sheet
  parameters_plot <- parameters
  parameters_plot$start_d <- 12 * 26
  
  # Apply to BOTH propIs_new and propIr_new
  parameters_plot$propIs_new <- c(D0, D1, D2)
  parameters_plot$propIr_new <- c(D0, D1, D2)
  
  # Storage per D2
  Y_inc_r  <- matrix(NA_real_, nrow = length(yrs), ncol = length(asym_seq))
  Y_GR_inc <- matrix(NA_real_, nrow = length(yrs), ncol = length(asym_seq))
  
  for (i in seq_along(asym_seq)) {
    asym <- asym_seq[i]
    prob_sym <- 1 - (asym / 100)
    
    parameters_plot$prob_sym_s <- prob_sym
    parameters_plot$prob_sym_r <- prob_sym
    
    out <- ode(y = init_state, times = sim_time,
               func = Malaria_model_with_Array, parms = parameters_plot,
               events = events)
    
    annual_inc_r <- colSums(matrix(out[, "inc_r"], nrow = 12))
    annual_GR_inc <- colSums(matrix(out[, "GR_inc"], nrow = 12))
    
    Y_inc_r[, i]  <- annual_inc_r
    Y_GR_inc[, i] <- annual_GR_inc
  }
  
  # Table for this D2
  table_like <- data.frame(
    Asym = asym_seq,
    c_beta = if (length(parameters_plot$c_beta) == 1) parameters_plot$c_beta else NA_real_,
    D0 = D0, D1 = D1, D2 = D2,
    Incidence_Res_year_1 = Y_inc_r[idx1, ],
    Incidence_Res_year_3 = Y_inc_r[idx3, ],
    red_inc_r_1y  = pct_red(annual_inc_r_base[idx1], Y_inc_r[idx1, ]),
    red_inc_r_3y  = pct_red(annual_inc_r_base[idx3], Y_inc_r[idx3, ]),
    Gametocyte_Res_year_1 = Y_GR_inc[idx1, ],
    Gametocyte_Res_year_3 = Y_GR_inc[idx3, ],
    red_GR_inc_1y = pct_red(annual_GR_inc_base[idx1], Y_GR_inc[idx1, ]),
    red_GR_inc_3y = pct_red(annual_GR_inc_base[idx3], Y_GR_inc[idx3, ])
  )
  
  # ใส่ baseline เป็นแถวแรก
  table_like <- rbind(row_baseline, table_like)
  setTxtProgressBar(pb, (D2+1)*10)
  # Write sheet
  sheet_name <- sprintf("D2_%0.1f", D2)
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, table_like)

}

# Save
saveWorkbook(wb, "Scenario_tables_by_D2.xlsx", overwrite = TRUE)
