# Run Simulation
library(deSolve)
setwd("D:/Work/2025/malaria/16 Jun/Baseline2_new_model/lambda_sym_asym")
#Read model
source("model_D0_33_D1_67_baseline_change_birth_death_2betaSets.r")

#
# source("function.r")

#
source("init_parameter_equilibrium_birth_2_rate_2beta_res_new.R")

totalpop <- read.csv("data/uganda_totalpop_2000_2100.csv", header = TRUE)
uganda_Incidence <- read.csv("data/uganda_Incidence.csv")

init_state <- readRDS("init_state.rds")

# parameters$beta_s <- rep(parameters_fit[1],3)
# parameters$beta_as <- parameters_fit[2]
parameters$beta_r <- parameters$beta_s_2 *0.625
parameters$beta_ar <-  parameters$beta_as_2 *0.625


target_sym_asym <- rep((25/75),276) # Target incidence for symptomatic asymptomatic cases
target_s_r <- rep((1.5),12) #

times <- seq(0,12*51)
# parameters$beta_s_2 <- rep(5.7,3)
# parameters$beta_s_3 <- rep(6.2,3)

# Poputation
totalpop_2000_2100 <- read.csv("data/uganda_totalpop_2000_2100.csv", header = TRUE)

#test model ode
out <- ode(y = init_state, times = times, func = Malaria_model_with_Array, parms = parameters)

# Plot the results inc_s and inc_r 0 : 6000
plot(out[,"inc"], type = "l", xlab = "Months", ylab = "Incidence Sensitive",
     main = "Incidence Sensitive", col = "purple", lwd = 2, ylim = c(0, 1.1*max(out[,"inc_s"])))
# lines(out[,"inc_r"], col = "red", lwd = 2)

plot(out[,"inc_sym"], type = "l", xlab = "Months", ylab = "Incidence Sensitive",
     main = "Incidence Sensitive", col = "purple", lwd = 2, ylim = c(0, 1.1*max(out[,"inc_s"])))

plot(seq(2000,2051,1/12), out[,"N"], type = "l",
     xlab = "Time (Years)", ylab = "Population",
     main = "Total Population from Model with Optimized Parameters",
     col = "blue", lwd = 2, xlim = c(2000,2050),
     ylim = c(20000000, max(out[,"N"])))
points(2000:2100, totalpop_2000_2100[,2], col = "red", pch = 19)

events <- list(
  func = function(time, state, parameters) {
    if (time == 12 * 14) {
      print(paste("Event triggered at time", time))
      state["Ir0"] <- state["Ir0"] + state["Is0"] * 0.05
      state["Ir1"] <- state["Ir1"] + state["Is1"] * 0.05
      state["Ir2"] <- state["Ir2"] + state["Is2"] * 0.05
      state["Ar"] <- state["Ar"] + state["As"] * 0.05
      
      # reduce Is
      state["Is0"] <- state["Is0"] *0.95
      state["Is1"] <- state["Is1"] *0.95
      state["Is2"] <- state["Is2"] *0.95
      state["As"] <- state["As"] *0.95
      
      # state["S"] <- state["S"] - (
      #   state["Is0"] * 0.05 +
      #     state["Is1"] * 0.05 +
      #     state["Is2"] * 0.05 +
      #     state["As"] * 0.05
      # )
    }
    return(state)
  },
  time = 12 * 14  # only one time point
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

#calculate the root mean square deviation
rmsd <- function(model,data){
  diff <- (model - data)^2
  mean_sq <- sum(diff)/length(data)
  return(sqrt(mean_sq))
}

parameters_fitting <- parameters
times_fit <- seq(1, 12*23, 1) # 12 months per year, 18 years

# Set the initial values for the parameters to be fitted
parameters_fitting$beta_s <- c(4.7,4.7,4.7)  # assuming same for all Is subgroups
parameters_fitting$beta_as <- 1.8  # assuming same as beta_s0
parameters_fitting$beta_s_2 <- c(3.9,3.9,3.9)  # assuming same for all Is subgroups
parameters_fitting$beta_as_2 <- 1.5  # assuming same as beta_s0
parameters_fitting$beta_r <- 2.4
parameters_fitting$beta_ar <- 1

#the objective function for using with the 'optim' command
loss_function <- function(pars) {
  print(pars)
  parameters_fitting$beta_s <- c(pars[1],pars[1],pars[1])
  parameters_fitting$beta_as <- pars[2]  # assuming same as beta_s0
  parameters_fitting$beta_s_2 <- c(pars[3],pars[3],pars[3])
  parameters_fitting$beta_as_2 <- pars[4]  # assuming same as beta_s0
  parameters_fitting$beta_r <- c(pars[5],pars[5],pars[5])
  parameters_fitting$beta_ar <- pars[6]  # assuming same as beta_s0
  if(pars[5]<pars[6]){
    return(10e8)
  }
  #init state is last_values
  init_state <- init_state
  
  out <- ode(y = init_state, times = times_fit,
             func = Malaria_model_with_Array, parms = parameters_fitting,
             events = events)
  # plot(seq(2000+1/12,2023,1/12),(out[,"inc_sym"] / out[,"inc_asym"]),
  #      type = "l", xlab = "Months", ylab = "Incidence", col="green",
  #      main = "Model Incidence vs Target Incidence",ylim=c(0,3))
  # lines(seq(2000+1/12,2023,1/12),target_sym_asym, col = "blue")
  # plot(out[,"inc_r"])
  if(length(out[,"inc"]) < 276){
    return(10e8)
  }
  #
  mat_inc <- matrix(out[,"inc"], nrow = 12)
  inc_model_2000_2022 <- (colSums(mat_inc))
  
  
  
  plot(2000:2022,inc_model_2000_2022,lwd=2,type = "l",ylim = c(0,1.1*max(uganda_Incidence[,2])))
  lines(2000:2022,uganda_Incidence[,2],col=2,lwd=2)
  
  #ratio inc_sym and inc_asym
  ratio_inc_s_r_last_12 <- out[(276-11):276,"inc_s"] / out[(276-11):276,"inc_r"]
  
  #rmsd for each component
  rmsd_inc_s <- rmsd(inc_model_2000_2022, uganda_Incidence[,2])
  rmsd_inc_s_r <- rmsd(ratio_inc_s_r_last_12,target_s_r)*10e4
  
  if(is.na(rmsd_inc_s) || is.nan(rmsd_inc_s) || is.infinite(rmsd_inc_s) ||
     is.na(rmsd_inc_s_r) || is.nan(rmsd_inc_s_r) || is.infinite(rmsd_inc_s_r)) {
    return(10e8) # If the sum of rmsd is too high, return a large value to penalize
  }
  
  
  print(paste0("Means Symptomatic/Asymptomatic Incidence Ratio: ", mean(out[,"inc_sym"] / out[,"inc_asym"])))
  #ratio of incidence
  print(sum(rmsd_inc_s,rmsd_inc_s_r))
  
  # root mean square error
  # sum of rmsd
  sum(rmsd_inc_s,rmsd_inc_s_r) # sum of rmsd
}
# Initial guess for the parameters
# initial_guess <- c(4.5,2,3, 1.25,2.5 , 1  )  # Initial guess for βs and βr
initial_guess <- c(5,2,4, 2,3 , 2 )  # Initial guess for βs and βr

result <- optim(
  par = initial_guess,
  fn = loss_function,
  method = "L-BFGS-B",
  lower = rep(0.001,length(initial_guess)),  # Lower bounds for βs and βr
  upper = rep(10,length(initial_guess)),
  hessian = T
)

result

# Extract the optimized parameters
optimized_params <- result$par
cat("Optimized Parameters:\n")
cat("βs:", optimized_params[1], "\n") #4.43311
cat("βas:", optimized_params[2], "\n") #1.694944
cat("βs_2:", optimized_params[3], "\n") #6
cat("βas_2:", optimized_params[4], "\n") #2.277
cat("βr:", optimized_params[5], "\n") #6
cat("βar:", optimized_params[6], "\n") #2.277

# result <-readRDS("optimized_params_uganda_2beta_betares.rds")
optimized_params <- result$par
parameters_optim <- parameters

# Run the model with optimized parameters
parameters_optim$beta_s <- rep(optimized_params[1], 3)  # assuming same for all Is subgroups
parameters_optim$beta_as <- optimized_params[2]  # assuming same as beta_s0
parameters_optim$beta_s_2 <- rep(optimized_params[3], 3)  # assuming same for all Is subgroups
parameters_optim$beta_as_2 <- optimized_params[4]  # assuming same as beta_s0
parameters_optim$beta_r <- rep(optimized_params[5], 3)  # assuming same for all Is subgroups
parameters_optim$beta_ar <- optimized_params[6]  # assuming same as beta_s0

#save the optimized parameters
saveRDS(result,"optimized_params_uganda_2beta_betares.rds")


time_optimized <- seq(0, 12*30, by = 1)  # Define the time points for optimized parameters

out_optimized <- ode(
  y = init_state,
  times = time_optimized,
  func = Malaria_model_with_Array,
  parms = parameters_optim,
  events=events
)

# divide each into 12
monthly <- rep(uganda_Incidence[,2] / 12, each = 1)  # 12 points per value
pop_monthly <- as.vector(sapply(monthly, function(x) rep(x, 12)))


# Plot the results with optimized parameters
plot(out_optimized[,"time"], out_optimized[,"inc"], type = "l",
     xlab = "Time (Months)", ylab = "Incidence", main = "Malaria Incidence with Optimized Parameters",
     col = "blue", lwd = 2)
points(pop_monthly)
legend("topleft", legend = c("Optimized Incidence", "Target Incidence"),
       col = c("blue", "red"), lwd = 2, cex = 0.8)

#inc_s and inc_r
plot(out_optimized[,"time"], out_optimized[,"inc_s"], type = "l",
     xlab = "Time (Months)", ylab = "Incidence Sensitive", main = "Malaria Incidence Sensitive",
     col = "purple", lwd = 2, ylim = c(0, 1.1*max(out_optimized[,"inc_s"])))
lines(out_optimized[,"time"], out_optimized[,"inc_r"], col = "red", lwd = 2)

#plot inc_sym and inc_asym
plot(out_optimized[,"time"], out_optimized[,"inc_asym"], type = "l",
     xlab = "Time (Months)", ylab = "Incidence Symptomatic", main = "Malaria Incidence Symptomatic",
     col = "purple", lwd = 2, ylim = c(0, 1.1*max(out_optimized[,"inc_asym"])))
lines(out_optimized[,"time"], out_optimized[,"inc_sym"], col = "red", lwd = 2)


parameters_optim <- parameters
parameters_optim$beta_r <- rep(optimized_params[1], 3)  # assuming same for all Is subgroups
parameters_optim$beta_ar <- optimized_params[2]  # assuming same as beta_s0
time_optimized <- seq(0, 12*51, by = 1)  # Define the time points for optimized parameters
library(manipulate)
manipulate(
  {
    parameters_optim$beta_s <- rep(beta_s, 3)  # assuming same for all Is subgroups
    parameters_optim$beta_as <- beta_as  # assuming same as beta_s0
    parameters_optim$beta_s_2 <- rep(beta_s_2, 3)  # assuming same for all Is subgroups
    parameters_optim$beta_as_2 <- beta_as_2  # assuming same as beta_s0
    parameters_optim$beta_r <- rep(beta_r, 3)  # assuming same for all Is subgroups
    parameters_optim$beta_ar <- beta_ar  # assuming same as beta_s0
    
    out_optimized <- ode(
      y = init_state,
      times = time_optimized,
      func = Malaria_model_with_Array,
      parms = parameters_optim,
      events = events
    )
    
    mat_inc <- matrix(out_optimized[1:612,"inc"], nrow = 12)
    inc_model_2000_2022 <- (colSums(mat_inc))
    
    plot(out_optimized[,"inc_s"], type = "l", xlab = "Months", ylab = "Incidence Sensitive",
         main = "Incidence Sensitive", col = "purple", lwd = 2, ylim = c(0, 1.1*max(out[,"inc_s"])))
    lines(out_optimized[,"inc_r"], col = "red", lwd = 2)
    
    plot(out_optimized[,"inc_asym"], type = "l", xlab = "Months", ylab = "Incidence ",
         main = "Incidence asymptomatic", col = "purple", lwd = 2, ylim = c(0, 1.1*max(out[,"inc_asym"])))
    lines(out_optimized[,"inc_sym"], col = "red", lwd = 2)

    plot(seq(2000,2050,1),inc_model_2000_2022,type = "l", xlab = "Months", ylab = "Incidence",
         main = "Malaria Incidence with Optimized Parameters",
         col = "blue", lwd = 2, ylim = c(0, 1.1*max(inc_model_2000_2022)))
    lines(2000:2022,uganda_Incidence[,2],col=2,lwd=2)
    # print(length(out_optimized[,"inc"]))
    ratio_inc_sym_asym <- out_optimized[,"inc_sym"] / out_optimized[,"inc_asym"]
    ratio_inc_s_r <- out_optimized[(276-12):276,"inc_s"] / out_optimized[(276-12):276,"inc_r"]
    print(mean(ratio_inc_sym_asym[1:(12*14)]))
    # print(mean(ratio_inc_sym_asym[(1+(12*9)):(12*14)]))
    print(mean(ratio_inc_sym_asym[(1+(12*14)):276]))
    print(mean(ratio_inc_s_r))
    
  },
  beta_s = slider(0.001,10,step=0.001,initial =4.5),
  beta_as = slider(0.001,10,step=0.001,initial =2),
  beta_s_2 = slider(0.001,10,step=0.001,initial =3),
  beta_as_2 = slider(0.001,10,step=0.001,initial =1.25),
  beta_r = slider(0.001,10,step=0.001,initial =2.5),
  beta_ar = slider(0.001,10,step=0.001,initial =1)
)

library(matrixStats)
info<-solve(result$hessian) # original solve(-fit$hessian)
prop_sigma<-sqrt(diag(info))
prop_sigma<-diag(prop_sigma)
upper<-result$par+1.96*prop_sigma
lower<-result$par-1.96*prop_sigma
interval<-data.frame(value=result$par, upper=rowMaxs(upper), lower=rowMins(lower))

# #Extract the Hessian matrix
# hessian_matrix <- result$hessian
# #Compute the variance-covariance matrix
# variance_covariance_matrix <- solve(hessian_matrix)
# #Compute standard errors (SEs):
# standard_errors <- sqrt(diag(variance_covariance_matrix))
# #Calculate the 95% confidence intervals
# estimates <- result$par
# lower <- estimates - 1.96 * standard_errors
# upper <- estimates + 1.96 * standard_errors
# 
# ci_95 <- data.frame(Estimate = estimates,
#                     SE = standard_errors,
#                     Lower95 = lower,
#                     Upper95 = upper)


#plot ci95
parameters_lower <- parameters_fitting
parameters_lower$beta_s <- c(interval$lower[1],interval$lower[1],interval$lower[1])
parameters_lower$beta_as <- interval$lower[2]  # assuming same as beta_s0
parameters_lower$beta_s_2 <- c(interval$lower[3],interval$lower[3],interval$lower[3])
parameters_lower$beta_as_2 <- lower[4]  # assuming same as beta_s0
parameters_lower$beta_r <- c(interval$lower[5],interval$lower[5],interval$lower[5])
parameters_lower$beta_ar <- interval$lower[6]  # assuming same as beta_s0

out_lower <- ode(y = init_state, times = times_fit,
                 func = Malaria_model_with_Array, parms = parameters_lower)

parameters_upper <- parameters_fitting
parameters_upper$beta_s <- c(interval$upper[1],interval$upper[1],interval$upper[1])
parameters_upper$beta_as <- interval$upper[2]  # assuming same as beta_s0
parameters_upper$beta_s_2 <- c(interval$upper[3],interval$upper[3],interval$upper[3])
parameters_upper$beta_as_2 <- interval$upper[4]  # assuming same as beta_s0
parameters_upper$beta_r <- c(interval$upper[5],interval$upper[5],interval$upper[5])
parameters_upper$beta_ar <- interval$upper[6]  # assuming same as beta_s0
out_upper <- ode(y = init_state, times = times_fit,
                 func = Malaria_model_with_Array, parms = parameters_upper)

# Plot the results with confidence intervals
plot(out_optimized[,"time"]/12, out_optimized[,"inc"], type = "l",
     xlab = "Time (Years)", ylab = "Population",
     main = "Total Population from Model with Confidence Intervals",
     col = "blue", lwd = 2, ylim = c(0, 1.1*max(out_optimized[,"inc"])))
lines(out_lower[,"time"]/12, out_lower[,"inc"], col = "green", lty = 2)
lines(out_upper[,"time"]/12, out_upper[,"inc"], col = "red", lty = 2)
#polygon
polygon(c(out_lower[,"time"]/12, rev(out_upper[,"time"]/12)), 
        c(out_lower[,"N"], rev(out_upper[,"N"])), 
        col = rgb(0, 1, 0, 0.2), border = NA)
points(0:24, totalpop[,2], col = "red", pch = 19)
# Add legend
legend("topright", legend = c("Optimized Population", "Lower CI", "Upper CI", "Observed Population"),
       col = c("blue", "green", "red", "red"), pch = c(NA, NA, NA, 19), lwd = 2, lty = c(1, 2, 2, NA), cex = 0.8)

