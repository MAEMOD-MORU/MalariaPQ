# plot Ci95 with Using hessian
library(deSolve)
library(dplyr)

#Read model
source("../model/model_D0_33_D1_67_baseline_change_birth_death_2betaSets_SS.R")

#
source("init_parameter_calibration_infection.R")
Tanzania_data <- read.csv("../data/Reported malaria cases by method of confirmation.csv")

op<- readRDS("optimized_params_equilibrium.rds")
init_state <- readRDS("init_state_2000.rds")

parameters_fitted <- parameters
parameters_fitted$beta_s <- rep(op$par[1],3)
parameters_fitted$beta_as <- op$par[1]

times <- seq(1,26*12)
year <- 2000:2025

model_fitted <-  ode(y = init_state, times = times, func = Malaria_model_with_Array, parms = parameters_fitted)
model_fitted_inc <- colSums(matrix(model_fitted[,"inc"],   nrow = 12))

plot(year,model_fitted_inc,type = "l",ylim = c(0,1e7))
points(Tanzania_Incidence[,1],Tanzania_Incidence[,2],pch=16)

