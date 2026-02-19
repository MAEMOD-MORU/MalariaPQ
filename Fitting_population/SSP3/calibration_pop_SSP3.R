library(deSolve)
library(ggplot2)

source("init_pop_calibration.R")
source("model_D0_33_D1_67_baseline_change_birth_death_2betaSets_SS.R")
### population calibration #####
# data population
# 2000	2005	2010	2015	2020	2025	2030	2035	2040	2045	2050
totalpop_data <- read.csv("data/Tanzania_pop_SSP_2000_2100.csv", header = TRUE)[1:11,4]

#pop 2000
rmsd <- function(totalpop_model){
  diff <- (totalpop_model - totalpop_data)^2
  mean_sq <- sum(diff)/length(totalpop_data)
  return(sqrt(mean_sq))
}
times <- seq(0,12*51)
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
  points(ssp_year,totalpop_data)
  rmsd(totalpop_model[indx]) 
}

op <- optim(c(0.035,0.00825),pop_run,
            hessian= T)
op$par
parameters_fitted <- parameters

parameters_fitted$mui <- op$par[1]
parameters_fitted$muo <- op$par[2]
times <- seq(0,12*50)
model_fitted <-  ode(y = state, times = times, func = Malaria_model_with_Array, parms = parameters_fitted)

plot(year,model_fitted[,"N"],type = "l")
points(ssp_year,totalpop_data)

# save op
saveRDS(op, file = "calibration_population_results_SSP3.rds")
