# Run Simulation
library(deSolve)
library(BayesianTools)
library(scales)
#Read model
source("../model/model_D0_33_D1_67_baseline_change_birth_death_2betaSets_SS.R")
source("init_parameter_optim_fitted.R")

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



initial_guess <- c(1.4,1.3,1.05, 0.73)

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
plotDiagnostic(MCMC_out)

saveRDS(MCMC_out, "MCMC_out.rds")

# if Gelman Rubin multivariate psrf:  > 1.1 
# it indicates that the chains have not converged well 
# and may require more iterations or adjustments to the sampling process.
# In such cases, you might consider increasing the number of iterations,
# adjusting the proposal distribution, or checking for any issues 
# in the model specification that could be affecting convergence.

# we will run more iterations until the multivariate psrf is less than 1.1 to 
# ensure better convergence of the MCMC chains.
# MCMC_out_2nd <- runMCMC(MCMC_out)

# if you want new settings, you can specify them in the runMCMC function, for example:
# MCMC_out_2nd <- runMCMC(MCMC_out, settings = list(iterations = 15000, burnin = 5000))


# summary(MCMC_out_2nd)
# plot(MCMC_out_2nd)
# gelmanDiagnostics(MCMC_out_2nd)
# correlationPlot(MCMC_out_2nd,start=1000,scaleCorText =F)
# marginalPlot(MCMC_out_2nd,start=1000)
# plotDiagnostic(MCMC_out_2nd)
# 
# saveRDS(MCMC_out_2nd, "MCMC_out_2nd.rds")

# if it don't work or psrf don't go down
# you can try to adjust the proposal distribution,
# for example by changing the scale of the proposal distribution 
# or using a different sampling method.
# exmaple:
# density = function(par){
# d.beta1 <-  dunif(par[1], min = 0.01, max = 3, log =TRUE) 
# ...
# sampler = function(n=1){
# s.beta1 <- runif(n, min = 0.01, max = 3)
# ...
# change 5 to 3 in the prior distribution to make it narrower, 
# which can help with convergence if the original range was too wide.
# or change prior lower, upper, or best values to better reflect 
# the expected parameter space based on prior knowledge or preliminary runs.
# prior <- createPrior(density = density, sampler = sampler,
# lower = c(rep(0.1, 2),1,0.1,1),
# upper = c(rep(3, 2),1.5,1,50),
# best = c(c(1.5,1.2,1.1, 0.75))