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

# Run forward simulations from each threshold state
params_fwd        <- parameters
params_fwd$start_d <- 1000
params_fwd$start_b <- 1000
time_fwd          <- seq(0, 12*10 - 1)

run_scenario <- function(init) {
  out <- ode(y = init, times = time_fwd, func = Malaria_model_with_Array, parms = params_fwd)
  list(
    inc   = colSums(matrix(out[,"inc"],   nrow = 12)),
    inc_s = colSums(matrix(out[,"inc_s"], nrow = 12)),
    inc_r = colSums(matrix(out[,"inc_r"], nrow = 12)),
    ratio = colMeans(matrix(out[,"inc_r"] / out[,"inc"], nrow = 12))
  )
}

results <- lapply(init_states, run_scenario)

# Plot incidence resistant
colors <- c("black","purple","red","blue")
labels <- paste("Initial Proportion Resistant", thresholds)

plot(results[[1]]$inc_r, type="l", col=colors[1], lwd=2,
     ylim=c(0, 1.1*max(sapply(results, function(r) max(r$inc_r)))),
     xlab="Years", ylab="Incidence", main="Incidence Resistant")
for (i in 2:4) lines(results[[i]]$inc_r, col=colors[i], lwd=2)
legend("topleft", legend=labels, col=colors, lwd=2, cex=1)

# Plot resistance ratio
plot(results[[1]]$ratio, type="l", col=colors[1], lwd=2, ylim=c(0,1),
     xlab="Years", ylab="Ratio", main="Incidence Resistant / Total Incidence")
for (i in 2:4) lines(results[[i]]$ratio, col=colors[i], lwd=2)
for (i in 1:4) abline(h=thresholds[i], col=colors[i], lty=2)
legend("bottomright", legend=labels, col=colors, lwd=2, cex=1)

thresholds <- c(0.01, 0.05, 0.1, 0.15)
colors     <- c("0.01" = "#4472C4", "0.05" = "#FF0000", 
                "0.10" = "#00B050", "0.15" = "#000000")

# --- params 2 แบบ ---
params_no_switch        <- parameters
params_no_switch$start_d <- 1000   # ไม่ switch
params_no_switch$start_b <- 1000

params_switch           <- parameters
params_switch$start_d   <- 1      # switch ทันที (เส้นจุด)
params_switch$start_b   <- 1000

time_new <- seq(1, 12*10)

# --- รัน ทั้ง 2 แบบ ---
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

# เส้นทึบ (no switch)
df_solid <- bind_rows(mapply(
  function(init, label) run_ratio_df(init, params_no_switch, time_new, label, "solid"),
  init_list, threshold_labels, SIMPLIFY = FALSE
))

# เส้นจุด (switch at t=1)
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
  scale_linetype_identity() +   # ใช้ค่า "solid"/"dashed" ตรงๆ
  scale_x_continuous(breaks = 1:10, expand = c(0.01, 0)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0, 0)) +
  labs(
    x = "Years",
    y = "Proportion resistant\n(symptomatic incidence)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major  = element_line(color = "grey85", linewidth = 0.4),
    panel.grid.minor  = element_blank(),
    legend.position   = "right",
    legend.title      = element_text(face = "bold", size = 10),
    legend.text       = element_text(size = 10),
    legend.key.width  = unit(1.8, "cm"),
    axis.title        = element_text(size = 11),
    axis.text         = element_text(size = 10)
  ) +
  guides(
    color    = guide_legend(override.aes = list(linewidth = 1.5)),
    linetype = "none"   # ไม่แสดง linetype legend แยก
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
