library(deSolve)

# --- Load model and parameters ---
source("../model/model_D0_33_D1_67_baseline_change_birth_death_2betaSets_SS.R")
source("init_parameter_MCMC.R")

init_state <- readRDS("init_state_2000.rds")
times_fit  <- seq(1, 12*51, 1)

events <- list(
  func = function(time, state, parameters) {
    if (time == 1) {
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

run_annual <- function(pars) {
  out <- ode(y = init_state, times = times_fit,
             func = Malaria_model_with_Array,
             parms = pars, events = events)
  list(
    inc_r = colSums(matrix(out[, "inc_r"], nrow = 12)),
    GR    = colSums(matrix(out[, "GR_inc"], nrow = 12))
  )
}

pct_red <- function(base_val, scen_val) {
  if (is.na(base_val) || base_val == 0) return(NA_real_)
  (base_val - scen_val) / base_val * 100
}

i26 <- 27; i28 <- 29  # 2026=index27, 2028=index29

asym_levels     <- seq(0, 100, by = 10)
coverage_levels <- seq(0, 100, by = 10)

results <- list()
pb <- txtProgressBar(min = 0, max = length(asym_levels) * length(coverage_levels), style = 3)
k  <- 0

for (asym in asym_levels) {
  prob_sym <- 1 - asym / 100
  
  # ✅ Baseline = scenario ที่ cov=0% (D2=0, switches at 2025 แต่ไม่มี D2)
  pars_cov0 <- parameters
  pars_cov0$prob_sym_s <- prob_sym
  pars_cov0$prob_sym_r <- prob_sym
  pars_cov0$propIs_new <- c(1, 0, 0)   # D0=100%, no D1, no D2
  pars_cov0$propIr_new <- c(1, 0, 0)
  pars_cov0$start_d    <- (2025 - 2000) * 12
  
  base <- run_annual(pars_cov0)
  
  for (cov in coverage_levels) {
    k <- k + 1
    
    D2_frac <- cov / 100
    D0_frac <- 1 - D2_frac
    
    pars_scen <- parameters
    pars_scen$prob_sym_s <- prob_sym
    pars_scen$prob_sym_r <- prob_sym
    pars_scen$propIs_new <- c(D0_frac, 0, D2_frac)
    pars_scen$propIr_new <- c(D0_frac, 0, D2_frac)
    pars_scen$start_d    <- (2025 - 2000) * 12
    
    scen <- run_annual(pars_scen)
    
    results[[k]] <- data.frame(
      Asym          = asym,
      D2_coverage   = cov,
      red_inc_r_1y  = round(pct_red(base$inc_r[i26], scen$inc_r[i26]), 2),
      red_inc_r_3y  = round(pct_red(base$inc_r[i28], scen$inc_r[i28]), 2),
      red_GR_inc_1y = round(pct_red(base$GR[i26],    scen$GR[i26]),    2),
      red_GR_inc_3y = round(pct_red(base$GR[i28],    scen$GR[i28]),    2)
    )
    
    setTxtProgressBar(pb, k)
  }
}
close(pb)

df_results <- do.call(rbind, results)

# cov=0 ต้องได้ 0% เสมอ (เทียบกับตัวเอง)
df_results <- df_results %>%
  mutate(across(starts_with("red_"), ~ ifelse(D2_coverage == 0, 0, .)))

# =============================================
# Wide tables
# =============================================
library(tidyr)
library(dplyr)
library(openxlsx)

make_wide <- function(data, var) {
  data %>%
    select(Asym, D2_coverage, all_of(var)) %>%
    pivot_wider(names_from  = D2_coverage,
                values_from = all_of(var),
                names_prefix = "Cov_") %>%
    arrange(Asym) %>%
    rename(`Asym (%)` = Asym)
}

tbl_inc_r_1y  <- make_wide(df_results, "red_inc_r_1y")
tbl_inc_r_3y  <- make_wide(df_results, "red_inc_r_3y")
tbl_GR_inc_1y <- make_wide(df_results, "red_GR_inc_1y")
tbl_GR_inc_3y <- make_wide(df_results, "red_GR_inc_3y")

cat("\n=== % reduction in resistant incidence (1-year) ===\n");  print(tbl_inc_r_1y)
cat("\n=== % reduction in resistant incidence (3-year) ===\n");  print(tbl_inc_r_3y)
cat("\n=== % reduction in GR incidence (1-year) ===\n");         print(tbl_GR_inc_1y)
cat("\n=== % reduction in GR incidence (3-year) ===\n");         print(tbl_GR_inc_3y)

wb <- createWorkbook()
for (nm in c("red_inc_r_1y","red_inc_r_3y","red_GR_inc_1y","red_GR_inc_3y")) {
  addWorksheet(wb, nm)
  writeData(wb, nm, make_wide(df_results, nm))
}
saveWorkbook(wb, "model_wide_tables.xlsx", overwrite = TRUE)
