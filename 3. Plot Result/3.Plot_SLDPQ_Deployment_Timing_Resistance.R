library(deSolve)
library(BayesianTools)
library(scales)
library(ggplot2)
library(dplyr)
library(tidyr)
library(openxlsx)
# Load data
source("model_D0_33_D1_67_baseline_change_birth_death_2betaSets_SS.R")
source("init_parameter_MCMC.R")
init_state           <- readRDS("init_state_2000.rds")
times                <- seq(1, 12*51)

# Run baseline ODE
out <- ode(y = init_state, times = times, func = Malaria_model_with_Array, parms = parameters)

# Event: swap 1 person Is1 -> Ir1 at time=1
events <- list(
  func = function(time, state, parameters) {
    if (time == 1) { state["Ir1"] <- state["Ir1"] + 1; state["Is1"] <- state["Is1"] - 1 }
    return(state)
  },
  time = 1
)
out_event <- ode(y = init_state, times = times, func = Malaria_model_with_Array,
                 parms = parameters, events = events)

# Find time points where rolling mean of resistance ratio hits thresholds
roll_mean <- sapply(1:(nrow(out_event)-11), function(i) mean(out_event[i:(i+11), "inc_r"] / out_event[i:(i+11), "inc"]))
thresholds <- c(0.01, 0.05, 0.1, 0.15)
pos        <- sapply(thresholds, function(th) which.min(abs(roll_mean - th)))
names(pos) <- paste0("p", gsub("\\.", "_", thresholds))

# Extract initial state at each threshold
state_vars <- c("S", "Is0","Is1","Is2","Stis1","Stis2","Fis1","Fis2",
                "GIs0","GIs1","GIs2","GAs","As","Rs",
                "Ir0","Ir1","Ir2","Stir1","Stir2","Fir1","Fir2",
                "GIr0","GIr1","GIr2","GAr","Ar","Rr")

extract_state <- function(out, row) {
  sv <- out[row, state_vars, drop = FALSE]
  setNames(unname(unlist(sv)), state_vars)
}

init_states <- lapply(pos, function(p) extract_state(out_event, p))

recheck_init_state <- function(out_event, center_pos, threshold,
                               params_run, search_range = 12) {
  
  time_check <- seq(1, 12)
  n          <- nrow(out_event)
  
  candidates <- seq(
    max(1, center_pos - search_range),
    min(n, center_pos + search_range)
  )
  
  best_idx  <- center_pos
  best_diff <- Inf
  best_mean <- NA
  
  for (t_idx in candidates) {
    
    init_try <- extract_state(out_event, t_idx)
    
    out_try <- tryCatch(
      ode(y     = init_try,
          times = time_check,
          func  = Malaria_model_with_Array,
          parms = params_run),
      error = function(e) NULL
    )
    
    if (is.null(out_try)) next
    
    # mean inc_sym_r / inc_sym ใน 12 เดือนแรก
    ratio_12  <- out_try[, "inc_sym_r"] / out_try[, "inc_sym"]
    mean_12   <- mean(ratio_12, na.rm = TRUE)
    diff_12   <- abs(mean_12 - threshold)
    
    if (diff_12 < best_diff) {
      best_diff <- diff_12
      best_idx  <- t_idx
      best_mean <- mean_12
    }
  }
  
  message(sprintf(
    "Threshold=%.2f | original_pos=%d | best_pos=%d | mean_12=%.4f | diff=%.4f",
    threshold, center_pos, best_idx, best_mean, best_diff
  ))
  
  extract_state(out_event, best_idx)
}

# =============================================
# recheck ทุก threshold
# =============================================
params_check        <- parameters
params_check$start_d <- 1000
params_check$start_b <- 1000

init_states_checked <- mapply(
  function(p, th) recheck_init_state(out_event, 
                                     center_pos  = p, 
                                     threshold   = th,
                                     params_run  = params_check),
  pos, thresholds,
  SIMPLIFY = FALSE
)
names(init_states_checked) <- names(pos)

# =============================================
# ใช้ init_states_checked แทน init_states
# =============================================
init_list <- list(
  init_states_checked$p0_01,
  init_states_checked$p0_05,
  init_states_checked$p0_1,
  init_states_checked$p0_15
)

# solid (no switch)
df_solid <- bind_rows(mapply(
  function(init, label) run_ratio_df(init, params_no_switch, time_new, label, "solid"),
  init_list, threshold_labels, SIMPLIFY = FALSE
))

# dashed (switch at t=1)
df_dash <- bind_rows(mapply(
  function(init, label) run_ratio_df(init, params_switch, time_new, label, "dashed"),
  init_list, threshold_labels, SIMPLIFY = FALSE
))

df_all <- bind_rows(df_solid, df_dash)
df_all$Threshold <- factor(df_all$Threshold, levels = threshold_labels)

time_new <- seq(1, 12*10)

run_ratio_df <- function(init_st, pars, time_run, label, linetype) {
  out <- ode(y     = init_st,
             times = time_run,
             func  = Malaria_model_with_Array,
             parms = pars)
  
  mat   <- matrix(out[, "inc_sym_r"] / out[, "inc_sym"], nrow = 12)
  ratio <- colMeans(mat)
  
  data.frame(
    Year      = 1:10,
    Ratio     = ratio,
    Threshold = label,
    Linetype  = linetype
  )
}

threshold_labels <- c("0.01", "0.05", "0.10", "0.15")
init_list        <- list(init_states$p0_01, init_states$p0_05,
                         init_states$p0_1,  init_states$p0_15)

# (no D2)
df_solid <- bind_rows(mapply(
  function(init, label) run_ratio_df(init, params_no_switch, time_new, label, "solid"),
  init_list, threshold_labels, SIMPLIFY = FALSE
))

# (start_d at t=1)
df_dash <- bind_rows(mapply(
  function(init, label) run_ratio_df(init, params_switch, time_new, label, "dashed"),
  init_list, threshold_labels, SIMPLIFY = FALSE
))

df_all <- bind_rows(df_solid, df_dash)
df_all$Threshold <- factor(df_all$Threshold, levels = threshold_labels)

# --- สี ---
colors_manual <- c(
  "0.01" = "#4169E1",
  "0.05" = "#CC0000",
  "0.10" = "#228B22",
  "0.15" = "black"
)

# --- Plot ---
ggplot(df_all, aes(x = Year, y = Ratio,
                   color    = Threshold,
                   linetype = Linetype,
                   group    = interaction(Threshold, Linetype))) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(
    values = colors_manual,
    name   = "SLDPQ Deployment Threshold\n(Symptomatic Resistant Fraction)"
  ) +
  scale_linetype_identity() +
  scale_x_continuous(
    breaks       = 1:10,
    minor_breaks = seq(1, 10, 0.5),   # minor tick แกน X
    expand       = c(0.01, 0)
  ) +
  scale_y_continuous(
    limits       = c(0, 1),
    breaks       = seq(0, 1, 0.2),
    minor_breaks = seq(0, 1, 0.05),   # minor tick ทุก 0.05 — ครอบ 0.05, 0.10, 0.15
    expand       = c(0, 0)
  ) +
  labs(
    x = "Years",
    y = "Proportion resistant\n(symptomatic incidence)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major   = element_line(color = "grey85", linewidth = 0.4),
    panel.grid.minor   = element_blank(),          # ไม่แสดง grid minor
    legend.position    = "right",
    legend.title       = element_text(face = "bold", size = 10),
    legend.text        = element_text(size = 10),
    legend.key.width   = unit(1.8, "cm"),
    axis.title         = element_text(size = 11),
    axis.text          = element_text(size = 10),
    axis.ticks         = element_line(color = "black"),
    axis.ticks.length  = unit(0.15, "cm"),          # major tick
    axis.minor.ticks.length = rel(0.5)              # minor tick สั้นกว่า major
  ) +
  guides(
    color    = guide_legend(override.aes = list(linewidth = 1.5)),
    linetype = "none",
    x        = guide_axis(minor.ticks = TRUE),      # เปิด minor tick แกน X
    y        = guide_axis(minor.ticks = TRUE)       # เปิด minor tick แกน Y
  )

# --- 2 sheet ---
df_no_int <- df_all %>%
  filter(Linetype == "solid") %>%
  mutate(Ratio = round(Ratio, 4)) %>%
  select(Threshold, Year, Ratio) %>%
  pivot_wider(
    names_from  = Year,
    values_from = Ratio,
    names_prefix = "Year_"
  ) %>%
  arrange(Threshold) %>%
  rename(`Start Resistance Ratio` = Threshold)

df_sldpq <- df_all %>%
  filter(Linetype == "dashed") %>%
  mutate(Ratio = round(Ratio, 4)) %>%
  select(Threshold, Year, Ratio) %>%
  pivot_wider(
    names_from  = Year,
    values_from = Ratio,
    names_prefix = "Year_"
  ) %>%
  arrange(Threshold) %>%
  rename(`Start Resistance Ratio` = Threshold)

# --- Export ---
wb <- createWorkbook()

addWorksheet(wb, "No_SLDPQ")
writeData(wb, "No_SLDPQ", df_no_int)

addWorksheet(wb, "With_SLDPQ")
writeData(wb, "With_SLDPQ", df_sldpq)

saveWorkbook(wb, "resistance_ratio_by_scenario.xlsx", overwrite = TRUE)
