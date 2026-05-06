# ============================================================
# Run Simulation
# ============================================================
library(deSolve)
library(scales)
library(openxlsx)

# Read model
source("../../model/model_D0_33_D1_67_baseline_change_birth_death_2betaSets_SS.R")
source("../init_parameter_MCMC.R")

Tanzania_k13_Allele_frequency <- read.csv("../../data/Tanzania_K13_Allele_frequency.csv")[, 2]
init_state <- readRDS("../init_state_2000.rds")
times <- seq(1, 12 * 51)

# ============================================================
# Function: extract initial state from ODE output at a given position
# ============================================================
get_init_state <- function(ode_output, position) {
  select_values <- ode_output[position, , drop = FALSE]
  
  new_init <- c(
    S = unname(select_values[,"S"]),
    Is0 = unname(select_values[,"Is0"]),
    Is1 = unname(select_values[,"Is1"]),
    Is2 = unname(select_values[,"Is2"]),
    Stis1 = unname(select_values[,"Stis1"]),
    Stis2 = unname(select_values[,"Stis2"]),
    Fis1 = unname(select_values[,"Fis1"]),
    Fis2 = unname(select_values[,"Fis2"]),
    GIs0 = unname(select_values[,"GIs0"]),
    GIs1 = unname(select_values[,"GIs1"]),
    GIs2 = unname(select_values[,"GIs2"]),
    GAs = unname(select_values[,"GAs"]),
    As = unname(select_values[,"As"]),
    Rs = unname(select_values[,"Rs"]),
    Ir0 = unname(select_values[,"Ir0"]),
    Ir1 = unname(select_values[,"Ir1"]),
    Ir2 = unname(select_values[,"Ir2"]),
    Stir1 = unname(select_values[,"Stir1"]),
    Stir2 = unname(select_values[,"Stir2"]),
    Fir1 = unname(select_values[,"Fir1"]),
    Fir2 = unname(select_values[,"Fir2"]),
    GIr0 = unname(select_values[,"GIr0"]),
    GIr1 = unname(select_values[,"GIr1"]),
    GIr2 = unname(select_values[,"GIr2"]),
    GAr = unname(select_values[,"GAr"]),
    Ar = unname(select_values[,"Ar"]),
    Rr = unname(select_values[,"Rr"])
  )
}

# redistribute move all asymptomatic resistant to symptomatic resistant, 
# while keeping the same total number of resistant infections. 
redistribute_init <- function(init,text_show=F){
  
  # ratios before
  r_before  <- init["Ir0"] / init["Ir1"]
  gr_before <- init["GIr0"] / init["GIr1"]
  
  # redistribute
  r <- init["Ir0"] / (init["Ir0"] + init["Ir1"])
  init["Ir0"] <- init["Ir0"] + init["Ar"] * r
  init["Ir1"] <- init["Ir1"] + init["Ar"] * (1 - r)
  init["Ar"]  <- 0
  
  gr <- init["GIr0"] / (init["GIr0"] + init["GIr1"])
  init["GIr0"] <- init["GIr0"] + init["GAr"] * gr
  init["GIr1"] <- init["GIr1"] + init["GAr"] * (1 - gr)
  init["GAr"]  <- 0
  
  # ratios after
  r_after  <- init["Ir0"] / init["Ir1"]
  gr_after <- init["GIr0"] / init["GIr1"]
  if(text_show){
  cat("Ir ratio before:", r_before, " after:", r_after, "\n")
  cat("GIr ratio before:", gr_before, " after:", gr_after, "\n")
  }
  return(init)
}

# ============================================================
# Run ODE with resistance introduction event at time = 1
# ============================================================
events <- list(
  func = function(time, state, parameters) {
    if (time == 1) {
      state["Ir1"] <- state["Ir1"] + 1
      state["Is1"] <- state["Is1"] - 1
    }
    return(state)
  },
  time = 1
)

out_event <- ode(y = init_state, times = times,
                 func = Malaria_model_with_Array,
                 parms = parameters, events = events)

# ============================================================
# Compute incidence and resistance ratio
# ============================================================
mat_inc    <- matrix(out_event[, "inc"], nrow = 12)
inc_model  <- colSums(mat_inc)

inc_r_ratio <- out_event[, "inc_r"] / out_event[, "inc"]

# Rolling 12-month mean of resistance ratio
roll_mean <- sapply(1:(length(inc_r_ratio) - 11),
                    function(i) mean(inc_r_ratio[i:(i + 11)]))

# ============================================================
# Find positions closest to target resistance ratios
# ============================================================
target_ratios <- c(0.01, 0.05, 0.10, 0.15)
positions     <- sapply(target_ratios, function(t) which.min(abs(roll_mean - t)))
names(positions) <- paste0("pos_", gsub("\\.", "_", target_ratios))

# Preview roll_mean at each position
position_summary <- data.frame(
  target    = target_ratios,
  position  = positions,
  roll_mean = roll_mean[positions],
  inc_r_ratio_at_pos = inc_r_ratio[positions]
)
print(position_summary)

# ============================================================
# Function: re-run ODE from a given position
# ============================================================
run_recheck_position <- function(target,
                                 range_search = 0) {
  
  # Phase 1: Find initial position closest to target via roll_mean
  best_pos <- which.min(abs(roll_mean - target))
  message("Phase 1 | Initial position: ", best_pos,
          " | roll_mean = ", round(roll_mean[best_pos], 6),
          " | target = ", target)
  
  # Phase 2: Search within [best_pos - range_search, best_pos + range_search]
  search_positions <- (best_pos - range_search):(best_pos + range_search)
  search_positions <- search_positions[search_positions >= 1 &
                                         search_positions <= nrow(out_event)]
  
  run_single <- function(pos) {
    parameters_phase2         <- parameters
    parameters_phase2$start_d <- 1000
    parameters_phase2$start_b <- 1
    
    # Apply redistribute_init before evaluating ratio
    init_redistributed <- redistribute_init(get_init_state(out_event, pos))
    
    out_phase2 <- ode(y     = init_redistributed,
                      times = 0:11,
                      func  = Malaria_model_with_Array,
                      parms = parameters_phase2)
    
    mat_ratio <- matrix(out_phase2[, "inc_r"] / out_phase2[, "inc"], nrow = 12)
    
    list(
      mean       = mean(colMeans(mat_ratio)),
      init_state = get_init_state(out_phase2, nrow(out_phase2))  # ← last row ของ out_phase2
    )
  }
  
  results      <- lapply(search_positions, run_single)
  result_means <- sapply(results, function(x) x$mean)
  
  mean_start_pos <- result_means[search_positions == best_pos]
  best_in_range  <- which.min(abs(result_means - target))
  final_pos      <- search_positions[best_in_range]
  
  message("Phase 2 | Search window: [", min(search_positions), ", ", max(search_positions), "]",
          " | mean at start pos (", best_pos, ") = ", round(mean_start_pos, 3),
          " | Best position: ", final_pos,
          " | result mean = ", round(result_means[best_in_range], 3))
  
  return(get_init_state(out_event, final_pos))
}

#### Find ±5 positions around the initial position.####
init_state_0_01 <- run_recheck_position(target = 0.01, range_search = 5)
init_state_0_05 <- run_recheck_position(target = 0.05, range_search = 5)
init_state_0_1 <- run_recheck_position(target = 0.1, range_search = 5)
init_state_0_15 <- run_recheck_position(target = 0.15, range_search = 5)

init_state_0_01_symR <- redistribute_init(init_state_0_01)
init_state_0_05_symR <- redistribute_init(init_state_0_05)
init_state_0_1_symR <- redistribute_init(init_state_0_1)
init_state_0_15_symR <- redistribute_init(init_state_0_15)

# ============================================================
# create results table for different D2 values,
#### starting from the same initial state ####
# ============================================================
# =========================
# Parameters
# =========================

parameters_new_init <- parameters
parameters_new_init$start_d <- 1
parameters_new_init$start_b <- 1

time_new <- seq(0, 12 * 10 - 1)

# =========================
# Function: summarize model output
# =========================

summarize_output <- function(out_event) {
  
  mat_inc <- matrix(out_event[, "inc"], nrow = 12)
  inc_model <- colSums(mat_inc)
  
  mat_inc <- matrix(out_event[, "inc_s"], nrow = 12)
  inc_s_model <- colSums(mat_inc)
  
  mat_inc <- matrix(out_event[, "inc_r"], nrow = 12)
  inc_r_model <- colSums(mat_inc)
  
  mat_inc <- matrix(out_event[, "inc_sym_r"], nrow = 12)
  inc_sym_r_model <- colSums(mat_inc)
  
  mat_inc <- matrix(out_event[, "inc_asym_r"], nrow = 12)
  inc_asym_r_model <- colSums(mat_inc)
  
  mat_inc <- matrix(out_event[, "GR_inc"], nrow = 12)
  inc_GR_model <- colSums(mat_inc)
  
  mat_ratio <- matrix(out_event[, "inc_r"] / out_event[, "inc"], nrow = 12)
  ratio_model <- colMeans(mat_ratio)
  
  list(
    inc_model = inc_model,
    inc_s_model = inc_s_model,
    inc_r_model = inc_r_model,
    inc_sym_r_model = inc_sym_r_model,
    inc_asym_r_model = inc_asym_r_model,
    inc_GR_model = inc_GR_model,
    ratio_model = ratio_model
  )
}

# =========================
#### Function: change D2 ####
# =========================

update_prop_params <- function(parameters, d2) {
  parameters$propIs_new <- c(1 - d2, 0, d2)
  parameters$propIr_new <- c(1 - d2, 0, d2)
  return(parameters)
}

# =========================
#### Function: run one full scenario####
# =========================

run_simulation_set <- function(init_state_object) {
  
  results <- list()
  
  # Baseline: start_d = 6000
  parameters_baseline <- parameters
  parameters_baseline$start_d <- 6000
  parameters_baseline$start_b <- 1
  
  out_baseline <- ode(
    y = init_state_object,
    times = time_new,
    func = Malaria_model_with_Array,
    parms = parameters_baseline
  )
  
  results[["baseline"]] <- summarize_output(out_baseline)
  
  # D2 values
  d2_values <- c(0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80)
  
  for (d2 in d2_values) {
    
    parameters_tmp <- update_prop_params(parameters_new_init, d2)
    
    out_event <- ode(
      y = init_state_object,
      times = time_new,
      func = Malaria_model_with_Array,
      parms = parameters_tmp
    )
    
    results[[paste0("d2_", d2)]] <- summarize_output(out_event)
  }
  
  return(results)
}

# =========================
#### Run all initial states ####
# =========================
results_0_01 <- run_simulation_set(init_state_0_01_symR)
results_0_05 <- run_simulation_set(init_state_0_05_symR)
results_0_1  <- run_simulation_set(init_state_0_1_symR)
results_0_15 <- run_simulation_set(init_state_0_15_symR)

# =========================
#### Plotting ####
# =========================
i <- 1
years <- 1:length(results_0_05$baseline$ratio_model)
cols <- rainbow(length(results_0_05) - 1)
plot(years,
     rep(NA, length(years)),
     type = "l",
     col = "black",
     lwd = 3,
     ylim = c(0, 1),
     xlab = "Year",
     ylab = "Resistance Ratio",
     main = "Resistance ratio under different scenarios")
for (name in c("d2_0.3","d2_0.8")) {
  if (name != "baseline") {
    lines(years, results_0_01[[name]]$ratio_model,
          col = 1, lwd = 2, lty = i)
    lines(years, results_0_05[[name]]$ratio_model,
          col = 2, lwd = 2, lty = i)
    
    lines(years, results_0_1[[name]]$ratio_model,
          col = 3, lwd = 2, lty = i)
    
    lines(years, results_0_15[[name]]$ratio_model,
          col = 4, lwd = 2, lty = i)
    i <- i + 1
  }
}
abline(h = 0.4, col = "gray", lty = 2)
legend("topleft",
       legend = c(
                  paste0("D2 30%, init 0.01"),
                  paste0("D2 80%, init 0.01"),
                  paste0("D2 30%, init 0.05"),
                  paste0("D2 80%, init 0.05"),
                  paste0("D2 30%, init 0.10"),
                  paste0("D2 80%, init 0.10"),
                  paste0("D2 30%, init 0.15"),
                  paste0("D2 80%, init 0.15")),
       col = c(rep(c(1:4), each = 2)),
       lty = c(1:2),
       lwd = 2,
       cex = 0.7)

# =========================
### Save separate Excel files ####
# =========================

vars <- c(
  "inc_r_model",
  "inc_asym_r_model",
  "inc_GR_model",
  "ratio_model")
all_results <- list(
  "results_0_01" = results_0_01,
  "results_0_05" = results_0_05,
  "results_0_1"  = results_0_1,
  "results_0_15" = results_0_15
)

# If you have results_0_01, add it:
# all_results[["results_0_01"]] <- results_0_01

rename_scenario <- function(name) {
  
  if (name == "baseline") {
    return("Baseline (D1 33%)")
  }
  
  val <- as.numeric(sub("d2_", "", name))
  percent <- val * 100
  
  paste0("D2 ", percent, "%")
}

for (res_name in names(all_results)) {
  
  wb <- createWorkbook()
  res <- all_results[[res_name]]
  
  for (v in vars) {
    
    n <- length(res[[1]][[v]])
    
    df <- data.frame(
      year = 1:n,
      c_beta = rep(0.745, n)
    )
    
    for (scenario in names(res)) {
      
      new_name <- rename_scenario(scenario)
      df[[new_name]] <- res[[scenario]][[v]]
    }
    
    addWorksheet(wb, sheetName = v)
    writeData(wb, sheet = v, x = df)
  }
  
  saveWorkbook(
    wb,
    file = paste0(res_name, ".xlsx"),
    overwrite = TRUE
  )
}