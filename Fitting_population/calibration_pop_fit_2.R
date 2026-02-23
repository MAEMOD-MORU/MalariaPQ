library(deSolve)
library(ggplot2)

source("init_parameter_fitted.R")
source("model_D0_33_D1_67_baseline_change_birth_death_2betaSets_SS.R")
### population calibration #####
# data population
# 2001-2100
totalpop_data <- read.csv("data/Tanzania_pop_2000_2100.csv", header = TRUE)[1:51,2]
# 2000	2005	2010	2015	2020	2025	2030	2035	2040	2045	2050
totalpop_ssp <- read.csv("data/Tanzania_pop_SSP_2000_2100.csv", header = TRUE)[1:11,2:6]

#pop 2000
rmsd <- function(totalpop_model){
  diff <- (totalpop_model - totalpop_data)^2
  mean_sq <- sum(diff)/length(totalpop_data)
  return(sqrt(mean_sq))
}
times <- seq(0,12*51)
#pop 2000-2100
pop_run_exp <- function(pars, dat) {
  parameters$mui <- exp(pars[1])
  parameters$muo <- exp(pars[2])
  
  out_model <- tryCatch(
    ode(
      y = state, times = times,
      func = Malaria_model_with_Array,
      parms = parameters,
      rtol = 1e-8, atol = 1e-10, maxsteps = 5e5
    ),
    error = function(e) return(NULL)
  )
  
  if (is.null(out_model)) return(1e9)
  
  # penalize if solver returned early
  if (nrow(out_model) < length(times)) return(1e9)
  
  totalpop_model <- out_model[, "N"]
  
  # yearly indices (assuming monthly with initial time included)
  indx <- seq(1, 12*51, by = 12)
  
  # safety if your time grid differs
  if (max(indx) > length(totalpop_model)) return(1e9)
  
  val <- rmsd(totalpop_model[indx])
  if (!is.finite(val)) val <- 10e9
  print(val)
  val
}

# op <- optim(
#   par = c(-3.755, -3.895),
#   pop_run_exp,
#   hessian=TRUE,
#   control=list(reltol=1e-2, maxit=20)
# )

op <- optim(
  par = c(-3.5, -3.8),
  pop_run_exp,
  hessian=TRUE,
  control=list(reltol=1e-3, maxit=30)
)
op$par

saveRDS(op,"Population_fitted_optim_2.rds")
op<- readRDS("calibration_population_results.rds")
parameters_fitted <- parameters

vcov <- solve(op$hessian)

se <- sqrt(diag(vcov))

lower <- op$par - 1.96*se
upper <- op$par + 1.96*se

conf_intervals <- data.frame(
  Estimate = op$par,
  Lower_95 = lower,
  Upper_95 = upper
)


interval_max_min <- exp(conf_intervals)
interval_max_min

# --- after you already have: op, vcov_matrix, state, times, parameters, Malaria_model_with_Array,
#     totalpop_data, years, times_select, etc.

years <- 2000:2050
times_select <- seq(1, 12*51, 12)  # your yearly indices in the ODE output


# # B) Percentage parameter perturbation:
f <- 0.01  # 20% (try 0.1, 0.2, 0.3)
mui_hat <- exp(op$par[1])
muo_hat <- exp(op$par[2])
mui_low  <- mui_hat * (1 - f)
mui_high <- mui_hat * (1 + f)
muo_low  <- muo_hat * (1 - f) 
muo_high <- muo_hat * (1 + f)

## 3) Convert to natural scale (mui/muo)
mui_low  <- exp(lower_par[1])
mui_high <- exp(upper_par[1])
muo_low  <- exp(lower_par[2])
muo_high <- exp(upper_par[2])

## 4) Run ODEs for fitted / upper / lower
# fitted
parameters_fitted <- parameters
parameters_fitted$mui <- exp(op$par[1])
parameters_fitted$muo <- exp(op$par[2])
model_fitted <- ode(y = state, times = times, func = Malaria_model_with_Array, parms = parameters_fitted)
model_fitted_pop <- model_fitted[times_select, "N"]

# upper envelope: high birth (mui high) + low death (muo low)
parameters_upper <- parameters
parameters_upper$mui <- mui_high
parameters_upper$muo <- muo_low
model_upper <- ode(y = state, times = times, func = Malaria_model_with_Array, parms = parameters_upper)
model_upper_pop <- model_upper[times_select, "N"]

# lower envelope: low birth (mui low) + high death (muo high)
parameters_lower <- parameters
parameters_lower$mui <- mui_low
parameters_lower$muo <- muo_high
model_lower <- ode(y = state, times = times, func = Malaria_model_with_Array, parms = parameters_lower)
model_lower_pop <- model_lower[times_select, "N"]

## 5) Make sure upper is really above lower (avoid polygon crossing)
upper_pop <- pmax(model_upper_pop, model_lower_pop)
lower_pop <- pmin(model_upper_pop, model_lower_pop)

## 6) Plot with shaded band
ylim_use <- range(c(totalpop_data, lower_pop, upper_pop, model_fitted_pop), na.rm = TRUE)

plot(years, model_fitted_pop, type = "n",
     xlab = "years", ylab = "model_fitted_pop",
     ylim = ylim_use)

polygon(c(years, rev(years)),
        c(upper_pop, rev(lower_pop)),
        col = rgb(0.5, 0.5, 0.5, 0.25),
        border = NA)

lines(years, model_fitted_pop, lwd = 2)
lines(years, upper_pop, col = "gray50", lty = 2)
lines(years, lower_pop, col = "gray50", lty = 2)

points(years, totalpop_data, pch = 16)

legend("topleft",
       legend = c("Fitted", "Envelope (upper/lower)", "Data"),
       lty = c(1, 2, NA), lwd = c(2, 1, NA), pch = c(NA, NA, 16),
       bty = "n")

