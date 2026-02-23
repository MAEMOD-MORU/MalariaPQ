# Run Simulation
library(deSolve)
library(BayesianTools)
library(scales)
#Read model
source("model_D0_33_D1_67_baseline_change_birth_death_2betaSets_SS.R")
source("init_parameter_fitted.R")

Tanzania_data <- read.csv("data/Reported malaria cases by method of confirmation.csv")
Tanzania_Incidence <- Tanzania_data[6:14,c(1,4)]
Tanzania_k13_Allele_frequency <- read.csv("data/Tanzania_K13_Allele_frequency.csv")[,2]

init_state <- readRDS("init_state_2000.rds") # form 1.2 Find_init_state_equilibrium_incidence

times <- seq(0,12*51)

# Poputation
totalpop_2000_2100 <- read.csv("data/Tanzania_pop_2000_2100.csv", header = TRUE)

#test model ode
out <- ode(y = init_state, times = times, func = Malaria_model_with_Array, parms = parameters)

#
plot(out[,"time"],out[,"N"], type = "l", xlab = "Months", ylab = "Incidence Sensitive",
     main = "Incidence Sensitive", col = "purple", lwd = 2, ylim = c(0, 1.1*max(out[,"N"])))

# Plot the results inc_s and inc_r 0 : 6000
plot(out[,"time"],out[,"inc"], type = "l", xlab = "Months", ylab = "Incidence Sensitive",
     main = "Incidence Sensitive", col = "purple", lwd = 2, ylim = c(0, 1.1*max(out[,"inc_s"])))
# lines(out[,"inc_r"], col = "red", lwd = 2)

plot(out[,"inc_sym"], type = "l", xlab = "Months", ylab = "Incidence Sensitive",
     main = "Incidence Sensitive", col = "purple", lwd = 2, ylim = c(0, 1.1*max(out[,"inc_s"])))

plot(seq(2000,2051,1/12), out[,"N"], type = "l",
     xlab = "Time (Years)", ylab = "Population",
     main = "Total Population from Model with Optimized Parameters",
     col = "blue", lwd = 2, xlim = c(2015,2050),
     ylim = c(20000000, max(out[,"N"])))
points(2000:2100, totalpop_2000_2100[,2], col = "red", pch = 19)

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

out_event <- ode(y = init_state, times = times, 
                 func = Malaria_model_with_Array, 
                 parms = parameters, events = events)


# Plot the results inc_s and inc_r 0 : 6000
plot(out_event[,"inc"], type = "l", xlab = "Months", ylab = "Incidence Sensitive",
     main = "Incidence Sensitive", col = "purple", lwd = 2, ylim = c(0, 1.1*max(out[,"inc_s"])))
# lines(out[,"inc_r"], col = "red", lwd = 2)

plot(out_event[,"inc_s"], type = "l", xlab = "Months", ylab = "Incidence Sensitive",
     main = "Incidence Sensitive", col = "purple", lwd = 2, ylim = c(0, 1.1*max(out[,"inc_s"])))
lines(out_event[,"inc_r"], col = "red", lwd = 2)
legend("topleft", legend = c("Incidence Sensitive", "Incidence Resistant"),
       col = c("purple", "red"), lwd = 2, cex = 0.8)

parameters_fitting <- parameters
times_fit <- seq(1, 12*26, 1) # 12 months per year, 18 years



initial_guess <- c(1.5,1.3,1.05, 0.75) 

parameters_fitting$beta_s   <- c(initial_guess[1], initial_guess[1], initial_guess[1])
parameters_fitting$beta_as  <- initial_guess[1]
parameters_fitting$beta_s_2 <- c(initial_guess[2], initial_guess[2], initial_guess[2])
parameters_fitting$beta_as_2<- initial_guess[2]
parameters_fitting$beta_r   <- c(initial_guess[1], initial_guess[1], initial_guess[1]) * initial_guess[3]
parameters_fitting$beta_ar  <- initial_guess[1] * initial_guess[3]
parameters_fitting$beta_r_2   <- c(initial_guess[2], initial_guess[2], initial_guess[2]) * initial_guess[3]
parameters_fitting$beta_ar_2  <- initial_guess[2] * initial_guess[3]
parameters_fitting$c_beta_r <- initial_guess[4]

parameters_fitting$start_b <- 12*21

out <- ode(y = init_state, times = times_fit,
           func = Malaria_model_with_Array, parms = parameters_fitting,
           events = events)

mat_inc <- matrix(out[,"inc"], nrow = 12)
inc_model_2015_2023 <- colSums(mat_inc)

plot(2000:2025, inc_model_2015_2023, type = "l",
     xlab = "Years", ylab = "Incidence",
     main = "Total Incidence by Year",
     # xlim = c(2014,2035),
     lwd = 2,
     col = 1, ylim = c(0, 1.1 * max(inc_model_2015_2023 *1.1)))
points(2015:2023,Tanzania_Incidence[, 2], col = "red",pch=19)

obs_inc <- Tanzania_Incidence[,2]
obs_k13 <- c(Tanzania_k13_Allele_frequency,0.4)

loglik_function <- function(pars) {
  parameters_fitting$beta_s   <- c(pars[1], pars[1], pars[1])
  parameters_fitting$beta_as  <- pars[1]
  parameters_fitting$beta_s_2 <- c(pars[2], pars[2], pars[2])
  parameters_fitting$beta_as_2<- pars[2]
  parameters_fitting$beta_r   <- c(pars[1], pars[1], pars[1]) * pars[3]
  parameters_fitting$beta_ar  <- pars[1] * pars[3]
  parameters_fitting$beta_r_2   <- c(pars[2], pars[2], pars[2]) * pars[3]
  parameters_fitting$beta_ar_2  <- pars[2] * pars[3]
  parameters_fitting$c_beta_r <- pars[4]
  par_phi <- pars[5]
  
  
  out <- ode(y = init_state, times = times_fit,
        func = Malaria_model_with_Array, parms = parameters_fitting,
        events = events)
  
  

  if (is.null(out) |  any(is.na(out)) |  any(is.nan(out))| any(is.infinite(out)) | length(out[,"inc"]) < 108) return(-1e8)
  
  mat_inc <- matrix(out[,"inc"], nrow = 12)
  inc_model_2000_2025 <- colSums(mat_inc)
  # x <- data.frame(year = 2000:2025,
  #            incidence = inc_model_2000_2025)
  
  mat_ratio <- matrix(out[,"inc_r"] / out[,"inc"], nrow = 12)
  ratio_2000_2025 <- colMeans(mat_ratio)
  
  par(mfrow = c(2, 1))
  
  years <- 2000:2025
  
  # Top: total incidence
  plot(
    years, inc_model_2000_2025,
    xlab = "Year", ylab = "Incidence",
    lwd  = 2, type = "l",
    main = "Total Incidence",
    xlim = c(2000, 2025),
    ylim = c(0, 2e7)
    # ylim = c(0, max(inc_model_2000_2025) * 1.1)
  )
  # lines(years, inc_model_2000_2025_r, col = "red")
  # lines(years, inc_model_2000_2025_n, col = "green3")
  points(2015:2023, Tanzania_Incidence[, 2], col = "blue", lwd = 2)
  
  # Bottom: resistance / total ratio
  plot(
    years, ratio_2000_2025,
    type = "l", col = "red",
    xlab = "Year", ylab = "Ratio",
    xlim = c(2000, 2025),
    ylim = c(0, 1),
    main = "Ratio Resistant Incidence / Total Incidence"
  )
  # abline(h = ref_val, lty = 2)  # 75% or 50% line
  points(
    c(2016:2022, 2025),
    c(Tanzania_k13_Allele_frequency, 0.4),
    col = 2, lwd = 2, pch = 3
  )
  abline(v = 2025, lty = 2)
  text(2023,0.4,round(ratio_2000_2025[26],2))
  
  # Poisson log-likelihood
  ll_inc <- sum(dpois(x = obs_inc, lambda = inc_model_2000_2025[16:24], log = TRUE))
  
  model_ratio <- c(ratio_2000_2025[17:23],ratio_2000_2025[26])
  
  # คำนวณ alpha และ beta สำหรับ dbeta โดยอิงจากค่าเฉลี่ย (p) ที่โมเดลทำนาย
  # shape1 = p * phi, shape2 = (1 - p) * phi
  shape1 <- model_ratio * par_phi
  shape2 <- (1 - model_ratio) * par_phi
  
  # คำนวณ Log-likelihood
  ll_k13 <- sum(dbeta(x = obs_k13, shape1 = shape1, shape2 = shape2, log = TRUE))
  
  ll <- ll_inc + ll_k13
  if(is.nan(ll) | is.infinite(ll)|is.na(ll)){
    ll <- -1e9
  }
  return(ll)
}

options(digits = 3)

density = function(par){
  d.beta1 <-  dunif(par[1], min = 0.01, max = 5, log =TRUE)
  d.beta2 <-  dunif(par[2], min = 0.01, max = 5, log =TRUE)
  d.beta_r_multiply <-  dunif(par[3], min = 1, max = 1.5, log =TRUE)
  d.c_beta_r <- dnorm(par[4], mean= 0.75, sd = 0.1, log =TRUE)
  d.phi <- dunif(par[5], min= 1, max = 25, log =TRUE)
  return(d.beta1+d.beta2+d.beta_r_multiply + d.c_beta_r + d.phi)
}

sampler = function(n=1){
  s.beta1 <- runif(n, min = 0.01, max = 5)
  s.beta2 <- runif(n, min = 0.01, max = 5)
  s.beta_r_multiply <- runif(n, min = 1, max = 1.5)
  s.c_beta_r <- rnorm(n, mean= 0.75, sd = 0.1)
  s.phi <- runif(n, min= 1, max = 25)
   return(cbind(s.beta1,
                 s.beta2,
                 s.beta_r_multiply,
                 s.c_beta_r,
                 s.phi))
}

prior <- createPrior(density = density, sampler = sampler,
                     lower = c(rep(0.01, 2),
                               1,0.01,
                               1
                               ),
                     upper = c(rep(5, 2),
                               1.5,1,
                               50
                               ),
                     best = c(initial_guess)
                     )

# par3 exponential distribution
# c_beta_r normal mean= 0.5 sd = 0.1 c_beta_r > 0

bayes_setup <- createBayesianSetup(
  likelihood = loglik_function,
  prior = prior
)

burnin <- 1500
iterations <- 7500

settings <- list(burnin=burnin,iterations = iterations+burnin, nrChains   = 3)
MCMC_out <- runMCMC(bayes_setup, settings = settings, sampler = "DEzs")

summary(MCMC_out)
plot(MCMC_out)
gelmanDiagnostics(MCMC_out)
correlationPlot(MCMC_out,start=1000,scaleCorText =F)
marginalPlot(MCMC_out,start=1000)
# plotDiagnostic(MCMC_out)


# saveRDS(MCMC_out5, "MCMC_out_08_02_2026_1_03.rds")

MCMC_out <- readRDS("MCMC_out_08_02_2026_1_03.rds")

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

# par(mfrow = c(1, 1))
# plot(2015:2023,inc_model_2015_2023,xlab= "month",ylab ="Incidence",
#      type = "l" ,col=2,ylim = c(0,1e7))
# for (i in 1:length(ysample)) {
#   x <- ysample[[i]]
#   lines(2015:2023,x,col=alpha('gray', 1))
# }
# lines(2015:2023,inc_model_2015_2023,col="blue",lwd=2)
# points(2015:2023,Tanzania_Incidence[, 2], col = "red",pch=19)

plot(2000:2050,inc_model_2000_2050,xlab= "month",ylab ="Incidence",
     type = "l" ,col=2,ylim = c(0,3e7))
for (i in 1:length(ysample)) {
  x <- ysample[[i]]
  lines(2000:2050,x,col=alpha('gray', 1))
}
lines(2000:2050,inc_model_2000_2050,col="blue",lwd=2)
points(2015:2023,Tanzania_Incidence[, 2], col = "red",pch=19)

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
write.csv(incidence_ci95, "incidence_ci95_09_feb.csv", row.names = FALSE)
resistance_ratio_ci95 <- data.frame(
  Year = times_year,
  Median = unlist( CI95_2[2,]),
  `2.5%` = low2,
  `97.5%` = up2
)
write.csv(resistance_ratio_ci95, "K13_resistance_ratio_09_feb.csv", row.names = FALSE)

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
write.csv(param_ci95_base, "parameter_ci95_09_feb.csv", row.names = FALSE)
