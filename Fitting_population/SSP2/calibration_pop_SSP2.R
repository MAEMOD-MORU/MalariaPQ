library(deSolve)
library(ggplot2)

source("init_pop_calibration.R")
source("model_D0_33_D1_67_baseline_change_birth_death_2betaSets_SS.R")
### population calibration #####
# data population
# 2000	2005	2010	2015	2020	2025	2030	2035	2040	2045	2050
totalpop_data <- read.csv("data/Tanzania_pop_SSP_2000_2100.csv", header = TRUE)[1:11,3]

#pop 2000
rmsd <- function(totalpop_model){
  diff <- (totalpop_model - totalpop_data)^2
  mean_sq <- sum(diff)/length(totalpop_data)
  return(sqrt(mean_sq))
}
times <- seq(0,12*50)
year <- seq(0,12*50/12,1/12)+2000
ssp_year <- seq(2000,2050,5)
#pop 2000-2050
pop_run <- function(pars){
  parameters$mui <- pars[1]
  parameters$muo <- pars[2]

  out_model <- ode(y = state, times = times, func = Malaria_model_with_Array, parms = parameters)
  
  totalpop_model <- out_model[,"N"]
  # totalpop_model <- colSums(pop_model)
  #2000 - 2050 = 50 year
  indx <- seq(0,12*50,12*5)+1 # index start 1
  # print(totalpop_model[indx])
  print(rmsd(totalpop_model[indx]))
  plot(year,totalpop_model,type = "l")
  points(ssp_year,totalpop_data)
  rmsd(totalpop_model[indx]) 
}
#optimization L-BFGS-B
# op <- optim(c(0.01,0.01),pop_run,
#             method = "L-BFGS-B",
#             lower = c(0.00001,0.00001),
#             upper = c(0.1,0.1),
#             hessian= T)
# op

op2 <- optim(c(0.022,0.022),pop_run,
            hessian= T)
parameters_fitted <- parameters

parameters_fitted$mui <- op2$par[1]
parameters_fitted$muo <- op2$par[2]

model_fitted <-  ode(y = state, times = times, func = Malaria_model_with_Array, parms = parameters_fitted)

plot(year,model_fitted[,"N"],type = "l")
points(ssp_year,totalpop_data)

# save op2
# saveRDS(op2, file = "calibration_population_results.rds")
fit <- readRDS("calibration_population_results.rds")
fisher_info<-solve(fit$hessian)
prop_sigma<-sqrt(diag(fisher_info))
prop_sigma<-diag(prop_sigma)
upper<-fit$par+0.01*prop_sigma
lower<-fit$par-0.01*prop_sigma
interval<-data.frame(value=fit$par, upper=upper, lower=lower)


# 1. Get the Variance-Covariance matrix
# solve() inverts the Hessian
vcov_matrix <- solve(fit$hessian)

# 2. Extract Standard Errors (square root of the diagonal)
se <- sqrt(diag(vcov_matrix))

# 3. Calculate 95% CI (using 1.96 for normal distribution)
lower <- fit$par - 1 * se
upper <- fit$par + 1 * se

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
h <- hessian(func = pop_run, x = fit$par)
se <- sqrt(diag(solve(h)))
conf_intervals <- data.frame(
  Estimate = fit$par,
  Lower_95 = fit$par - 1.96 * se,
  Upper_95 = fit$par + 1.96 * se
)
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
     lty=1,lwd = 2, ylim = c(0, 1.6*max(model_fitted[,"N"])))
lines(years,model_upper[,"N"], col = "blue", lty = 2)
lines(years,model_lower[,"N"], col = "blue", lty = 2)
# Shade the area between upper and lower bounds
x_vals <- c(years, rev((years)))
y_vals <- c(model_upper[,"N"], rev(model_lower[,"N"]))
polygon(x_vals, y_vals, col = rgb(0.5, 0.5, 0.5,alpha = 0.1), border = NA)
# add grid

lines(years,model_fitted[,"N"], col = "red", lwd = 2)
points(seq(2000,2050,1),totalpop_data, col='black', pch=16)
# add points ssp 1-3
points(seq(2000,2050,5),totalpop_ssp[,1], col='green4', pch=18, type = "b")
points(seq(2000,2050,5),totalpop_ssp[,2], col='brown4', pch=17, type = "b")
points(seq(2000,2050,5),totalpop_ssp[,3], col='purple3', pch=15, type = "b")
# abline(h = 1.3e8, col = "darkgreen", lty = 3)
legend("topleft", legend = c("Optimized Parameter", "95% CI", "Data","SSP1","SSP2","SSP3"),
       col = c("red", "blue", "black","green4","orange","purple3"), lwd = c(2,1,NA,1,1,1), 
       pch = c(NA,NA,16,18,17,15), ncol = 2,
       lty = c(1,2,NA,1,1,1), pt.cex = c(NA,NA,1,1,1,1))

# legend("topleft", legend = c("Fitted Parameter", "95% CI", "Data"),
#        col = c("red", "blue", "black"), lwd = c(2,1,NA), pch = c(NA,NA,16), 
#        lty = c(1,2,NA), pt.cex = c(NA,NA,1))