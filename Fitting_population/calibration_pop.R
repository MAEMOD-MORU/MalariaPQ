library(deSolve)
library(ggplot2)
library(dplyr)

source("init_pop_calibration.R")
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
pop_run <- function(pars,dat){
  parameters$mui <- pars[1]
  parameters$muo <- pars[2]
  times
  out_model <- ode(y = state, times = times, func = Malaria_model_with_Array, parms = parameters)
  
  totalpop_model <- out_model[,"N"]
  # totalpop_model <- colSums(pop_model)
  #2013 - 2035 = 22 year
  indx <- seq(1,12*51,12)
  # print(totalpop_model[indx])
  print(rmsd(totalpop_model[indx]))
  rmsd(totalpop_model[indx]) 
}
#optimization L-BFGS-B
op <- optim(c(0.035,0.008),pop_run,
            hessian= T)
parameters_fitted <- parameters

parameters_fitted$mui <- op2$par[1]
parameters_fitted$muo <- op2$par[2]

model_fitted <-  ode(y = state, times = times, func = Malaria_model_with_Array, parms = parameters_fitted)

plot(model_fitted[,"N"],type = "l")
points(seq(0,12*50,12),totalpop_data)

# save op
# saveRDS(op, file = "calibration_population_results.rds")
fit <- readRDS("calibration_population_results.rds")
fisher_info<-solve(fit$hessian)
prop_sigma<-sqrt(diag(fisher_info))
prop_sigma<-diag(prop_sigma)
upper<-fit$par+prop_sigma*1.96
lower<-fit$par-prop_sigma*1.96
interval<-data.frame(value=fit$par, upper=upper, lower=lower)

interval_max_min <- data.frame(
  value = fit$par,
  upper = apply(upper, 1, max),
  lower = apply(lower, 1, min)
)

parameters_upper <- parameters
parameters_upper$mui <- interval_max_min$upper[1]
parameters_upper$muo <- interval_max_min$upper[2]

model_upper <-  ode(y = state, times = times, func = Malaria_model_with_Array, parms = parameters_upper)

parameters_lower <- parameters
parameters_lower$mui <- interval_max_min$lower[1]
parameters_lower$muo <- interval_max_min$lower[2]

model_lower <-  ode(y = state, times = times, func = Malaria_model_with_Array, parms = parameters_lower)

plot(model_fitted[,"N"],type = "l",ylim = c(0,2e8))
lines(model_upper[,"N"],col=2)
lines(model_lower[,"N"],col=3)

# 1. Get the Variance-Covariance matrix
# solve() inverts the Hessian
vcov_matrix <- solve(fit$hessian)

# 2. Extract Standard Errors (square root of the diagonal)
se <- sqrt(diag(-vcov_matrix))

# 3. Calculate 95% CI (using 1.96 for normal distribution)
lower <- fit$par - 1.96 * se
upper <- fit$par + 1.96 * se

# 4. Combine into a table
conf_intervals <- data.frame(
  Estimate = fit$par,
  Lower_95 = lower,
  Upper_95 = upper
)

# If using a standard model object (glm, nls, etc.)
confint(fit, level = 0.95)

library(boot)

boot_func <- function(data, indices) {
  # 1. Create the resampled data 'd'
  d <- data[indices] 
  
  # 2. Inside the anonymous function, we use 'd' directly
  # We do NOT pass 'd' into pop_run() if pop_run doesn't want it
  fit_boot <- optim(
    par = c(0.035, 0.008), 
    fn = function(p) {
      # Use the global/parent 'd' here instead of passing it as an argument
      # This assumes your original pop_run() uses a variable named 'totalpop_data'
      totalpop_data <<- d 
      pop_run(p)
    }
  )
  return(fit_boot$par)
}

# Run the bootstrap
results <- boot(data = totalpop_data, statistic = boot_func, R = 1000)

# Get the 'Percentile' CI (usually wider and covers data better)
boot.ci(results, type = "perc", index = 1) # For parameter 1
boot.ci(results, type = "perc", index = 2) # For parameter 2

# 1. Create a time vector (assuming 51 years based on your data length)
years <- 1:length(totalpop_data)

# 2. Extract the Lower and Upper bounds from your boot.ci results
low_p  <- c(0.021, 0.0194)
high_p <- c(0.026, 0.0242)

# 3. Create a function to project population based on your model
# (Adjust this 'project' math to match whatever is inside your pop_run)
project <- function(p) {
  # Example: Simple exponential growth or your specific model
  # Replace this with your actual model formula
  totalpop_data[1] * exp(p[1] * years) 
}

# 4. Generate the curves
fit_low  <- project(low_p)
fit_high <- project(high_p)
fit_mid  <- project(c(0.0226, 0.0203))

# 5. Plot
plot(years, totalpop_data, pch=19, col="black", xlab="Time", ylab="Population")
lines(years, fit_mid, col="blue", lwd=2) # Best fit
lines(years, fit_low, col="red", lty=2)  # Lower bound
lines(years, fit_high, col="red", lty=2) # Upper bound
polygon(c(years, rev(years)), c(fit_low, rev(fit_high)), 
        col=rgb(1, 0, 0, 0.2), border=NA) # Shaded CI area

library(numDeriv)
# Replace 'your_likelihood_function' with your actual function name
# h <- hessian(func = pop_run, x = fit$par)
h <- readRDS("hessian_pop.RDS")
se <- sqrt(diag(solve(h)))
conf_intervals <- data.frame(
  Estimate = fit$par,
  Lower_95 = fit$par - 1.96 * se,
  Upper_95 = fit$par + 1.96 * se
)

library(boot)

totalpop_data_boot <- totalpop_data

rmsd <- function(x, y) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  stopifnot(length(x) == length(y))
  sqrt(mean((x - y)^2, na.rm = TRUE))
}

pop_run_boot <- function(pars,dat) {
  
  # Do NOT mutate a global parameters object; use a local copy
  pars_local <- parameters
  pars_local$mui <- pars[1]
  pars_local$muo <- pars[2]
  
  out_model <- ode(
    y     = state,
    times = times,
    func  = Malaria_model_with_Array,
    parms = pars_local
  )
  
  totalpop_model <- out_model[, "N"]
  
  # Your yearly-ish indices (Jan of each year if monthly)
  indx <- seq(1, 12 * 51, 12)
  
  # Safety: don't exceed length
  indx <- indx[indx <= length(totalpop_model) & indx <= length(dat)]
  
  # Key fix: model vs data
  print(rmsd(totalpop_model[indx], dat[indx]))
  rmsd(totalpop_model[indx], dat[indx])
}

boot_func <- function(data, indices, state, times, parameters) {
  
  # If 'data' is a numeric vector: resample elements
  d <- data[indices]
  
  # IMPORTANT: the model is still evaluated on the full times grid.
  # If you resample individual time points, you are breaking time alignment.
  # So we rebuild a full-length vector where non-selected points are kept,
  # OR (better) we should bootstrap *residuals* or use block bootstrap.
  #
  # Minimal "drop-in" approach: bootstrap residuals around the fitted mean is better,
  # but since you asked for full code, here is a simple (time-respecting) approach:
  #
  # Resample YEARS (blocks of 12) instead of individual months:
  # (See alternative section below.)
  
  fit_boot <- optim(
    par    = c(0.022, 0.022),
    fn     = pop_run_boot,
    dat    = d,
    method = "Nelder-Mead"
  )
  
  fit_boot$par
}

# Run 3 iterations
results <- boot(data = totalpop_data_boot, statistic = boot_func, R = 3)

# saveRDS(results, file = "boot_ci_results_pop.rds")

results <- readRDS("boot_ci_results_pop.rds")

# Get the 'Percentile' CI (usually wider and covers data better)
boot.ci(results, type = "perc", index = 1) # For parameter 1

p1_sims <- results$t[, 1]
p2_sims <- results$t[, 2]

# Calculate the 95% bounds (Percentile method)
p1_bounds <- quantile(p1_sims, probs = c(0.025, 0.975))
p2_bounds <- quantile(p2_sims, probs = c(0.025, 0.975))

#out with fitted parameters
parameters_fitted_final <- parameters
parameters_fitted_final$mui <- fit$par[1]
parameters_fitted_final$muo <- fit$par[2]
model_fitted <-  ode(y = state, times = times, func = Malaria_model_with_Array, parms = parameters_fitted_final)

# Simulate upper and lower bounds using the confidence intervals
params_upper <- parameters
params_upper$mui <- 0.023
params_upper$muo <- 0.02

model_upper <-  ode(y = state, times = times, func = Malaria_model_with_Array, parms = params_upper)

params_lower <- parameters
params_lower$mui <- 0.02
params_lower$muo <- 0.019

model_lower <-  ode(y = state, times = times, func = Malaria_model_with_Array, parms = params_lower)

years <- seq(0,612/12,1/12)+2000

#plot ci95 area with fitted par is line and add points of data and all line plot have points
plot(years,model_fitted[,"N"],type = "l", xlab = "Year", ylab = "Total Population",
     main = "Total Population Calibration with 95% CI", col = "red", 
     lty=1,lwd = 2, ylim = c(0, 1.5*max(model_fitted[,"N"])))
lines(years,model_upper[,"N"], col = "gray50", lty = 2)
lines(years,model_lower[,"N"], col = "gray50", lty = 2)
# Shade the area between upper and lower bounds
x_vals <- c(years, rev((years)))
y_vals <- c(model_upper[,"N"], rev(model_lower[,"N"]))
polygon(x_vals, y_vals, col = rgb(0.5, 0.5, 0.5,alpha = 0.1), border = NA)
# add grid
grid()
lines(years,model_fitted[,"N"], col = "red", lwd = 2)
points(seq(2000,2050,1),totalpop_data, col='black', pch=16)
# add points ssp 1-3
points(seq(2000,2050,5),totalpop_ssp[,1], col='green4', pch=18, type = "b")
points(seq(2000,2050,5),totalpop_ssp[,2], col='brown4', pch=17, type = "b")
points(seq(2000,2050,5),totalpop_ssp[,3], col='purple3', pch=15, type = "b")
# abline(h = 1.3e8, col = "darkgreen", lty = 3)
legend("topleft", legend = c("Optimized Parameter", "95% CI", "UN medium-fertility","SSP1","SSP2","SSP3"),
       col = c("red", "blue", "black","green4","orange","purple3"), lwd = c(2,1,NA,1,1,1), 
       pch = c(NA,NA,16,18,17,15), ncol = 2,
       lty = c(1,2,NA,1,1,1), pt.cex = c(NA,NA,1,1,1,1))
