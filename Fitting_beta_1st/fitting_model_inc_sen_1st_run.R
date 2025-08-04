# Run Simulation
library(deSolve)

#Read model
source("model_D0_33_D1_67_baseline_change_birth_death_2betaSets.R")

#
source("init_parameter_sen_only.R")
uganda_Incidence <- read.csv("data/uganda_Incidence.csv")

target_sym_asym <- rep((25/75),276) # Target incidence for symptomatic asymptomatic cases

#test model ode
out <- ode(y = state, times = times, func = Malaria_model_with_Array, parms = parameters)

# Plot the results inc_s and inc_r 0 : 6000
plot(out[,"inc_s"], type = "l", xlab = "Months", ylab = "Incidence Sensitive",
     main = "Incidence Sensitive", col = "purple", lwd = 2, ylim = c(0, 1.1*max(out[,"inc_s"])))
lines(out[,"inc_r"], col = "red", lwd = 2)

# Plot the results inc_s and inc_r 0 : 6000
plot(out[,"Is0"], type = "l", xlab = "Months", ylab = "Incidence Sensitive",
     main = "Incidence Sensitive", col = "purple", lwd = 2, ylim = c(0, 1.1*max(out[,"inc_s"])))
lines(out[,"As"], col = "red", lwd = 2)


# Plot the results N
plot(out[,"N"], type = "l", xlab = "Months", ylab = "Population",
     main = "Total Population", col = "blue", lwd = 2, ylim = c(24000100, 24000200))

#calculate the root mean square deviation
rmsd <- function(model,target){
  diff <- (model - target)^2
  mean_sq <- mean(diff)
  return(sqrt(mean_sq))
}

# Fitting function
# Define the loss function

parameters_fitting <- parameters # Copy the parameters for fitting
times <- seq(0,6000)

loss_function <- function(pars) {
  # print(beta_vec)
  parameters_fitting$beta_s <- c(pars[1],pars[1],pars[1])
  parameters_fitting$beta_as <- pars[2]  # assuming same as beta_s0
  
  #init state is last_values
  init_state <- state
  
  out <- ode(y = init_state, times = times,
             func = Malaria_model_with_Array, parms = parameters_fitting)
  plot(5725:6000,(out[5725:6000,"inc_sym"] / out[5725:6000,"inc_asym"]),
       type = "l", xlab = "Months", ylab = "Incidence", col="green",
       main = "Model Incidence vs Target Incidence",ylim=c(0,3))
  lines(5725:6000,target_sym_asym, col = "blue")
  
  #
  mat_inc_s <- matrix(out[5725:6000,"inc_s"], nrow = 12)
  inc_s_model_2000_2022 <- (colSums(mat_inc_s))
  
  plot(2000:2022,inc_s_model_2000_2022,lwd=2,type = "l",ylim = c(0,1.1*max(uganda_Incidence[,2])))
  lines(2000:2022,uganda_Incidence[,2],col=2,lwd=2)
  
  #ratio inc_sym and inc_asym
  ratio_inc_sym_asym <- out[5725:6000,"inc_sym"] / out[5725:6000,"inc_asym"]
  
  
  #rmsd for each component
  rmsd_inc_sym_asym <- rmsd(ratio_inc_sym_asym, target_sym_asym)*10e7
  rmsd_inc_s <- rmsd(inc_s_model_2000_2022, uganda_Incidence[,2])
  
  if(is.na(rmsd_inc_s) || is.nan(rmsd_inc_s) || is.infinite(rmsd_inc_s) ||
     is.na(rmsd_inc_sym_asym) || is.nan(rmsd_inc_sym_asym) || is.infinite(rmsd_inc_sym_asym)) {
    return(10e8) # If the sum of rmsd is too high, return a large value to penalize
  }
  
  print(paste0("Means Symptomatic/Asymptomatic Incidence Ratio: ", mean(out[5725:6000,"inc_sym"] / out[5725:6000,"inc_asym"])))
  #ratio of incidence
  print(sum(rmsd_inc_s,rmsd_inc_sym_asym))
  
  # root mean square error
  # sum of rmsd
  sum(rmsd_inc_s,rmsd_inc_sym_asym) # sum of rmsd
}
# Initial guess for the parameters
initial_guess <- c(6 , 2.277  )  # Initial guess for βs and βr
# Perform optimization
result <- optim(
  par = initial_guess,
  fn = loss_function,
  method = "L-BFGS-B",
  lower = rep(0.001,length(initial_guess)),  # Lower bounds for βs and βr
  upper = rep(100,length(initial_guess))
)
# Extract the optimized parameters
optimized_params <- result$par
cat("Optimized Parameters:\n")
cat("βs:", optimized_params[1], "\n") #2.82
cat("βas:", optimized_params[2], "\n") #1.07
parameters_optim <- parameters

# Run the model with optimized parameters
parameters_optim$beta_s <- rep(optimized_params[1], 3)  # assuming same for all Is subgroups
# parameters_optim$beta_r <- rep(optimized_params[2], 3)  # assuming same for all Ir subgroups
parameters_optim$beta_as <- optimized_params[2]  # assuming same as beta_s0
# parameters_optim$beta_ar <- optimized_params[4]  # assuming same as beta_r0

# Save the optimized parameters for back up
saveRDS(result,"optimized_params_result_uganda.rds")

time_optimized <- seq(0, 6000, by = 1)  # Define the time points for optimized parameters

out_optimized <- ode(
  y = state,
  times = time_optimized,
  func = Malaria_model_with_Array,
  parms = parameters_optim
)

last_values <- tail(out_optimized, 1)

# Initialize the last state for the next simulation
init_last_state <- c(
  S = last_values[,"S"],
  Is0 = last_values[,"Is0"],
  Is1 = last_values[,"Is1"],
  Is2 = last_values[,"Is2"],
  Stis1 = last_values[,"Stis1"],
  Stis2 = last_values[,"Stis2"],
  Fis1 = last_values[,"Fis1"],
  Fis2 = last_values[,"Fis2"],
  GIs0 = last_values[,"GIs0"],
  GIs1 = last_values[,"GIs1"],
  GIs2 = last_values[,"GIs2"],
  GAs = last_values[,"GAs"],
  As = last_values[,"As"],
  Rs = last_values[,"Rs"],
  Ir0 = last_values[,"Ir0"],
  Ir1 = last_values[,"Ir1"],
  Ir2 = last_values[,"Ir2"],
  Stir1 = last_values[,"Stir1"],
  Stir2 = last_values[,"Stir2"],
  Fir1 = last_values[,"Fir1"],
  Fir2 = last_values[,"Fir2"],
  GIr0 = last_values[,"GIr0"],
  GIr1 = last_values[,"GIr1"],
  GIr2 = last_values[,"GIr2"],
  GAr = last_values[,"GAr"],
  Ar = last_values[,"Ar"],
  Rr = last_values[,"Rr"]
)

# Save the last state for future use
saveRDS(init_last_state,"init_state.rds")

# Plot the results with optimized parameters
plot(out_optimized[5725:6000,"time"], out_optimized[5725:6000,"inc"], type = "l",
     xlab = "Time (Months)", ylab = "Incidence", main = "Malaria Incidence with Optimized Parameters",
     col = "blue", lwd = 2)
legend("topright", legend = c("Optimized Incidence", "Target Incidence"),
       col = c("blue", "red"), lwd = 2, cex = 0.8)



#plot ratio inc_sym and inc_asym
plot(out_optimized[5725:6000,"inc_sym"] / out_optimized[5725:6000,"inc_asym"], type = "l", xlab = "Months", ylab = "Incidence Ratio (Symptomatic/Asymptomatic)",
     main = "Incidence Ratio (Symptomatic/Asymptomatic)\n with Optimized Parameters", col = "green", lwd = 2, ylim = c(0, 0.5))
lines(target_sym_asym, col = "red", lwd = 2)
legend("bottomright", legend = c("Model", "Ratio 25/75(target)"),
       col = c("green", "red"), lwd = 2, cex = 0.8)

plot(out_optimized[5725:6000,"inc_asym"], type = "l", xlab = "Months", ylab = "Incidence Sensitive",
     main = "Incidence Sensitive (Sym/asym)", col = "green", lwd = 2, ylim = c(0, max(out_optimized[5725:6000,"inc_asym"])))
lines(out_optimized[5725:6000,"inc_sym"], col = "blue", lwd = 2)


paramter_optimized_drug <- parameters_optim
paramter_optimized_drug$start_d <- 120
times_optimized <- seq(0, 240, by = 1)  # Define the time points for optimized parameters with new start_d

#run model with optimized parameters and new start_d
out_optimized_drug <- ode(
  y = init_last_state,
  times = times_optimized,
  func = Malaria_model_with_Array,
  parms = paramter_optimized_drug
)
# Plot the results with optimized parameters and new start_d
plot(out_optimized_drug[,"time"], out_optimized_drug[,"inc"], type = "l",
     xlab = "Time (Months)", ylab = "Incidence", main = "Malaria Incidence with Optimized Parameters and New Start Day",
     col = "blue", lwd = 2)