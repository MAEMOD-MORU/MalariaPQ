# plot Ci95 with Using hessian
library(deSolve)
library(dplyr)

op<- readRDS("optimized_params_equilibrium.rds")

Tanzania_data <- read.csv("data/Reported malaria cases by method of confirmation.csv")
Tanzania_Incidence <- Tanzania_data[6:14,c(1,4)]

vcov <- solve(op$hessian)

se <- sqrt(diag(vcov))

lower <- op$par - 1.96*se
upper <- op$par + 1.96*se

conf_intervals <- data.frame(
  Estimate = op$par,
  Lower_95 = lower,
  Upper_95 = upper
)
conf_intervals

interval_max_min <- exp(conf_intervals)
interval_max_min

times <- seq(1,12*101)
# high birth and low death
parameters_upper <- parameters
parameters_upper$beta_s <- rep(conf_intervals$Upper_95[1],3)
parameters_upper$beta_as <- conf_intervals$Upper_95[1]

model_upper <-  ode(y = state, times = times, func = Malaria_model_with_Array, parms = parameters_upper)
model_upper_inc <- colSums(matrix(model_upper[,"inc"],   nrow = 12))

# low birth and high death
parameters_lower <- parameters
parameters_lower$beta_s <- rep(conf_intervals$Lower_95[1],3)
parameters_lower$beta_as <- conf_intervals$Lower_95[1]

model_lower <-  ode(y = state, times = times, func = Malaria_model_with_Array, parms = parameters_lower)
model_lower_inc <- colSums(matrix(model_lower[,"inc"],   nrow = 12))

parameters_fitted <- parameters
parameters_fitted$beta_s <- rep(conf_intervals$Estimate[1],3)
parameters_fitted$beta_as <- conf_intervals$Estimate[1]

model_fitted <-  ode(y = state, times = times, func = Malaria_model_with_Array, parms = parameters_fitted)
model_fitted_inc <- colSums(matrix(model_fitted[,"inc"],   nrow = 12))

# plot CI95
years <- 0:100

plot(years,model_fitted_inc,type = "l",ylim = c(0,1e7))
polygon(c(years, rev(years)), c(model_lower_inc, rev(model_upper_inc)), 
        col=rgb(0.5, 0.5, 0.5, 0.2), border=NA)
lines(years,model_upper_inc,col="gray50",lty=2)
lines(years,model_lower_inc,col="gray50",lty=2)
points(92:100,Tanzania_Incidence[,2],pch=16)

