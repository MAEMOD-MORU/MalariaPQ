# plot Ci95 with Using hessian
library(deSolve)

source("init_pop_calibration.R")
source("model_D0_33_D1_67_baseline_change_birth_death_2betaSets_SS.R")

op<- readRDS("calibration_population_results_SSP1.rds")

# 2000	2005	2010	2015	2020	2025	2030	2035	2040	2045	2050
totalpop_data <- read.csv("data/Tanzania_pop_SSP_2000_2100.csv", header = TRUE)[1:11,2]

vcov <- solve(op$hessian)

# se <- sqrt(diag(vcov))
# if In sqrt(diag(vcov)) : NaNs produced change to -vcov
se <- sqrt(diag(-vcov))
lower <- op$par - 1.96*se
upper <- op$par + 1.96*se

conf_intervals <- data.frame(
  Estimate = op$par,
  Lower_95 = lower,
  Upper_95 = upper
)
times <- seq(0,12*51)
times_select <- seq(1,12*51,12)
# high birth and low death
parameters_upper <- parameters
parameters_upper$mui <- conf_intervals$Upper_95[1]
parameters_upper$muo <- conf_intervals$Lower_95[2]

model_upper <-  ode(y = state, times = times, func = Malaria_model_with_Array, parms = parameters_upper)
model_upper_pop <- model_upper[times_select, "N"]

# low birth and high death
parameters_lower <- parameters
parameters_lower$mui <- conf_intervals$Lower_95[1]
parameters_lower$muo <- conf_intervals$Upper_95[2]

model_lower <-  ode(y = state, times = times, func = Malaria_model_with_Array, parms = parameters_lower)
model_lower_pop <- model_lower[times_select, "N"]

parameters_fitted <- parameters
parameters_fitted$mui <- conf_intervals$Estimate[1]
parameters_fitted$muo <- conf_intervals$Estimate[2]

model_fitted <-  ode(y = state, times = times, func = Malaria_model_with_Array, parms = parameters_fitted)
model_fitted_pop <- model_fitted[times_select, "N"]

# plot CI95
years <- 2000:2050
years_ssp <- seq(2000,2050,5)

plot(years,model_fitted_pop,type = "l",ylim = c(0,2e8))
polygon(c(years, rev(years)), c(model_lower_pop, rev(model_upper_pop)), 
        col=rgb(0.5, 0.5, 0.5, 0.2), border=NA)
lines(years,model_upper_pop,col="gray50",lty=2)
lines(years,model_lower_pop,col="gray50",lty=2)
points(years_ssp,totalpop_data,pch=16)

