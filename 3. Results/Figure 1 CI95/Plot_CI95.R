# Run Simulation
library(deSolve)
library(BayesianTools)
library(scales)
#Read model
source("../../model/model_D0_33_D1_67_baseline_change_birth_death_2betaSets_SS.R")
source("../init_parameter_MCMC.R")

Tanzania_data <- read.csv("../../data/Reported malaria cases by method of confirmation.csv")
Tanzania_Incidence <- Tanzania_data[6:14,c(1,4)]
Tanzania_k13_Allele_frequency <- read.csv("../../data/Tanzania_K13_Allele_frequency.csv")[,2]

init_state <- readRDS("../init_state_2000.rds") # form 1.2 Find_init_state_equilibrium_incidence

times <- seq(0,12*51)

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
  parameters$start_b <- 12*21
  
  return(parameters)
}

events <- list(
  func = function(time, state, parameters) {
    if (time == 1) {
      state["Ir0"] <- state["Ir0"] + 1
      state["Ir1"] <- state["Ir1"] + 1
      # state["Ir2"] <- state["Ir2"] + 1
      state["Ar"] <- state["Ar"] + 1
      
      # reduce Is
      state["Is0"] <- state["Is0"] -1
      state["Is1"] <- state["Is1"] -1
      # state["Is2"] <- state["Is2"] -1
      state["As"] <- state["As"] -1
      
    }
    return(state)
  },
  time =1  # only one time point
)

# Pull a sample of 1000 parameter sets from our Markov chain
paramsample_test <- getSample(MCMC_out, parametersOnly = T,
                              numSamples = 90,start = 1000 )

parameters_list <- list()

times_fit <- seq(1, 12*51, 1) # 12 months per year, 9 years-

for (i in 1:length(paramsample_test[,1])) {
  parameters_list[[i]] <- parameters_changing(parameters,paramsample_test[i,])
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

CI95 <- data.frame(getCredibleIntervals(params.rand_remove_nan, c(0.025, 0.5, 0.975)))
times_year <- 2000:2050
up  <- unlist(CI95[3,])
low <- unlist(CI95[1,])
med <- unlist(CI95[2,])

par(mfrow = c(1,1), mar = c(3, 3, 3, 1))  # leave right margin for legend

plot(times_year, med,
     type  = "l", lwd = 0,          # draw nothing yet; polygon goes first
     log   = "y",                   # <-- LOG SCALE
     xlab  = "Year",
     ylab  = "Malaria Incidence",
     ylim  = c(100, 3e7),
     xlim  = c(2000, 2050),
     panel.first = grid(col = "grey90", lty = 1))
polygon(c(times_year, rev(times_year)),
        c(low, rev(up)),
        col    = "grey80",
        border = NA)               # no border on polygon
lines(times_year, med, lwd = 2, col = "blue")
points(2015:2023, Tanzania_Incidence[, 2],
       col = "black", pch = 21, bg = "white",
       cex = 1.4, lwd = 1.5)
legend("bottomleft",
       legend = c("Confirmed Data", "Model Predicted Median"),
       pch    = c(1, NA),
       lty    = c(NA, 1),
       lwd    = c(1.5, 2),
       cex = 1,
       col    = c("black", "blue")
)
title("Malaria Incidence")

## --- Panel 2: Resistance Ratio ---
CI95_2 <- data.frame(getCredibleIntervals(params.rand_remove_nan2, c(0.025, 0.5, 0.975)))
up2  <- unlist(CI95_2[3,])
low2 <- unlist(CI95_2[1,])
med2 <- unlist(CI95_2[2,])

plot(times_year, med2,
     type = "l", lwd = 0,
     xlab = "Year",
     ylab = "Resistance Ratio",
     ylim = c(0, 1),
     xlim = c(2000, 2050),
     panel.first = grid(col = "grey90", lty = 1))
polygon(c(times_year, rev(times_year)),
        c(low2, rev(up2)),
        col    = "grey80",
        border = NA)
lines(times_year, med2, lwd = 2, col = "red") 
points(2016:2022, Tanzania_k13_Allele_frequency,
       col = "black", pch = 21, bg = "white",
       cex = 1.4, lwd = 1.5)
points(2025, 0.4,
       col = "black", pch = 21, bg = "white",
       cex = 1.4, lwd = 1.5)
legend("topleft",
       legend = c("Observed Data", "Model Predicted Median"),
       pch    = c(1, NA),
       lty    = c(NA, 1),
       lwd    = c(1.5, 2),
       col    = c("black", "red")
)
title("Resistance Ratio")

