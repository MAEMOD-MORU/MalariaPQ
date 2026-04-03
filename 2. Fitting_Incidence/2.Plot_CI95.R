# Run Simulation
library(deSolve)
library(BayesianTools)
library(scales)
#Read model
source("model_D0_33_D1_67_baseline_change_birth_death_2betaSets_SS.R")
source("init_parameter_optim_fitted.R")

Tanzania_data <- read.csv("data/Reported malaria cases by method of confirmation.csv")
Tanzania_Incidence <- Tanzania_data[6:14,c(1,4)]
Tanzania_k13_Allele_frequency <- read.csv("data/Tanzania_K13_Allele_frequency.csv")[,2]

init_state <- readRDS("init_state_2000.rds") # form 1.2 Find_init_state_equilibrium_incidence

times <- seq(0,12*51)

# Poputation
totalpop_2000_2100 <- read.csv("data/Tanzania_pop_2000_2100.csv", header = TRUE)


MCMC_out <- readRDS("MCMC_out_example.rds")

######Posterior Predictive Plot ######

parameters_changing <- function(parameters,pars){
  
  parameters$beta_s   <- c(pars[1], pars[1], pars[1])
  parameters$beta_as  <- pars[1]
  parameters$beta_s_2 <- c(pars[2], pars[2], pars[2])
  parameters$beta_as_2<- pars[2]
  parameters$beta_r   <- c(pars[1], pars[1], pars[1]) * pars[3]
  parameters$beta_ar  <- pars[1] * pars[3]
  parameters$beta_r_2   <- c(pars[2], pars[2], pars[2]) * pars[3]
  parameters$beta_ar_2  <- pars[2] * pars[3]
  parameters$c_beta_r <- pars[4]
  
  return(parameters)
}

# Pull a sample of 1000 parameter sets from our Markov chain
paramsample_test <- getSample(MCMC_out, parametersOnly = T,
                              numSamples = 300,start = 1000 )

parameters_list <- list()

times_fit <- seq(1, 12*51, 1) # 12 months per year, 9 years-

for (i in 1:length(paramsample_test[,1])) {
  parameters_list[[i]] <- parameters_changing(parameters_fitting,paramsample_test[i,])
}

pb <- txtProgressBar(min = 0, max = length(paramsample_test[,1]), style = 3)

ysample <- NULL
ysample2 <- NULL

for(i in 1:length(paramsample_test[,1])){
  parameters_run <- parameters_list[[i]]
  
  out_run <- ode(y = init_state, times = times_fit,
                 func = Malaria_model_with_Array, parms = parameters_run,
                 events = events)
  
  mat_inc <- matrix(out_run[,"inc"], nrow = 12)
  inc_model_2000_2050 <- colSums(mat_inc)
  
  mat_ratio <- matrix(out_run[,"inc_r"] / out_run[,"inc"], nrow = 12)
  ratio_2000_2050 <- colMeans(mat_ratio)
  
  if(i==1){
    ysample <- list(inc_model_2000_2050)
    ysample2 <- list(ratio_2000_2050) 
  }else{
    ysample <-append(ysample,list(inc_model_2000_2050))
    ysample2 <- append(ysample2,list(ratio_2000_2050))
  }
  setTxtProgressBar(pb, i)
}

length(ysample)

# incidence plot
plot(2000:2050,inc_model_2000_2050,xlab= "month",ylab ="Incidence",
     type = "l" ,col=2,ylim = c(0,3e7))
for (i in 1:length(ysample)) {
  x <- ysample[[i]]
  lines(2000:2050,x,col=alpha('gray', 1))
}
lines(2000:2050,inc_model_2000_2050,col="blue",lwd=2)
points(2015:2023,Tanzania_Incidence[, 2], col = "red",pch=19)

# K13 resistance ratio
plot(2000:2050,ratio_2000_2050,xlab= "month",ylab ="Resistance Ratio",
     type = "l" ,col=2,ylim = c(0,1))
for (i in 1:length(ysample2)) {
  x <- ysample2[[i]]
  lines(2000:2050,x,col=alpha('gray', 1))
}
lines(2000:2050,ratio_2000_2050,col="blue",lwd=2)
points(2016:2022,Tanzania_k13_Allele_frequency,col = "red",pch=19)
points(2025,0.4,col = "red",pch=19)

##### CI plot #####
# Incidence
params.rand <- matrix(data = NA,ncol = length(2000:2050),
                      nrow = nrow(paramsample_test))
params.rand2 <- matrix(data = NA,ncol = length(2000:2050),
                       nrow = nrow(paramsample_test))

for(i in 1:nrow(paramsample_test)){
  x <- ysample[[i]]
  
  params.rand[i,] <- rpois(n=2000:2050,
                           lambda = x
  )
  
  # 1. Get the predicted mean proportion from the model
  p_vector <- ysample2[[i]]
  
  # 2. Get the phi (precision) for this MCMC iteration
  # Replace 5 with the correct column index for phi in your MCMC chain
  phi_val <- paramsample_test[i, 5] 
  
  # 3. Constraint: p must be strictly between 0 and 1 for rbeta to work
  p_vector <- p_vector
  
  # 4. Sample from Beta distribution
  # shape1 = p * phi, shape2 = (1 - p) * phi
  params.rand2[i,] <- rbeta(n = length(2000:2050), 
                            shape1 = p_vector * phi_val, 
                            shape2 = (1 - p_vector) * phi_val)
}

# Remove Nan,NA,NULL
params.rand_complete <- complete.cases(params.rand)
params.rand_remove_nan <- params.rand[params.rand_complete,]
params.rand_complete2 <- complete.cases(params.rand2)
params.rand_remove_nan2 <- params.rand2[params.rand_complete2,]

CI95 <- data.frame(getCredibleIntervals(params.rand_remove_nan,c(0.025, 0.5,0.975)))

times_year <- 2000:2050

up <- unlist( CI95[3,])
low <- unlist( CI95[1,])
plot(times_year ,unlist( CI95[2,]),type = "l",xlab= "month",ylab ="Incidence",ylim = c(0, 3e+7))
#make polygon where coordinates start with lower limit and 
# then upper limit in reverse order
polygon(c(times_year,rev(times_year)),c(low,rev(up)),col = "grey75", border = FALSE)

lines(times_year ,unlist( CI95[2,]), lwd = 2,col=4)
#add red lines on borders of polygon
lines(times_year ,up, col="black",lty=2)
lines(times_year ,low, col="black",lty=2)
points(2015:2023 ,Tanzania_Incidence[, 2],col=2,type="p",pch=20,cex=2)
title("Incidence")

# Resistance Ratio
CI95_2 <- data.frame(getCredibleIntervals(params.rand_remove_nan2,c(0.025, 0.5,0.975)))
times_year <- 2000:2050                                                                    
up2 <- unlist( CI95_2[3,])
low2 <- unlist( CI95_2[1,])
plot(times_year ,unlist( CI95_2[2,]),type = "l",xlab= "month",ylab ="Resistance Ratio",ylim = c(0, 1))
#make polygon where coordinates start with lower limit and
# then upper limit in reverse order
polygon(c(times_year,rev(times_year)),c(low2,rev(up2)),col ="grey75", border = FALSE)
lines(times_year ,unlist( CI95_2[2,]), lwd = 2,col=4)
#add red lines on borders of polygon
lines(times_year ,up2, col="black",lty=2)
lines(times_year ,low2, col="black",lty=2)
points(2016:2022 ,Tanzania_k13_Allele_frequency,col= 2,type="p",pch=20,cex=2)
points(2025,0.4,col= 2,type="p",pch=20,cex=2)

# save csv CI95
incidence_ci95 <- data.frame(
  Year = times_year,
  Median = unlist( CI95[2,]),
  `2.5%` = low,
  `97.5%` = up
)
write.csv(incidence_ci95, "incidence_ci95.csv", row.names = FALSE)
resistance_ratio_ci95 <- data.frame(
  Year = times_year,
  Median = unlist( CI95_2[2,]),
  `2.5%` = low2,
  `97.5%` = up2
)
write.csv(resistance_ratio_ci95, "K13_resistance_ratio.csv", row.names = FALSE)

# If your sample has column names, this keeps them
param_names <- c("Beta 1",	"Beta 2",	"Multiply Beta R",	
                 "Resistant Fitness Cost",	"Phi For Beta Distribution")


param_ci95_base <- t(apply(paramsample_test, 2, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE))
param_ci95_base <- data.frame(
  Parameter = param_names,
  `2.5%` = param_ci95_base[,1],
  Median = param_ci95_base[,2],
  `97.5%`= param_ci95_base[,3],
  row.names = NULL
)

param_ci95_base
write.csv(param_ci95_base, "parameter_ci95.csv", row.names = FALSE)
