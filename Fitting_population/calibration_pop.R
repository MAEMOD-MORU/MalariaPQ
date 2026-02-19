library(deSolve)
library(dplyr)

source("init_pop_calibration.R")
source("model_D0_33_D1_67_baseline_change_birth_death_2betaSets_SS.R")
### population calibration #####
# data population
# 2001-2050
totalpop_data <- read.csv("data/Tanzania_pop_2000_2100.csv", header = TRUE)[1:51,2]
#pop 2000
rmsd <- function(totalpop_model){
  diff <- (totalpop_model - totalpop_data)^2
  mean_sq <- sum(diff)/length(totalpop_data)
  return(sqrt(mean_sq))
}
times <- seq(0,12*51)
#pop 2000-2050
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
op <- optim(c(0.035,0.00825),pop_run,
            hessian= T)
op$par
parameters_fitted <- parameters

parameters_fitted$mui <- op$par[1]
parameters_fitted$muo <- op$par[2]

model_fitted <-  ode(y = state, times = times, func = Malaria_model_with_Array, parms = parameters_fitted)

plot(model_fitted[,"N"],type = "l")
points(seq(0,12*50,12),totalpop_data)

# save op
saveRDS(op, file = "calibration_population_results.rds")
