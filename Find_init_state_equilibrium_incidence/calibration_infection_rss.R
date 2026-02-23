# Find init state 
library(deSolve)

#Read model
source("model_D0_33_D1_67_baseline_change_birth_death_2betaSets_SS.r")

#
source("init_parameter_calibration_infection.R")
Tanzania_data <- read.csv("data/Reported malaria cases by method of confirmation.csv")
Tanzania_Incidence <- Tanzania_data[6:14,c(1,4)]

target_sym_asym <- rep((25/75),276) # Target incidence for symptomatic asymptomatic cases

#test model ode
out <- ode(y = state, times = times, func = Malaria_model_with_Array, parms = parameters)

# Plot the results inc_s and inc_r 0 : 6000
plot(out[,"inc_s"], type = "l", xlab = "Months", ylab = "Incidence Sensitive",
     main = "Incidence Sensitive", col = "purple", lwd = 2, ylim = c(0, 1.1*max(out[,"inc_s"])))
lines(out[,"inc_r"], col = "red", lwd = 2)

# Plot the results N 
plot(out[,"N"], type = "l", xlab = "Months", ylab = "Population",
     main = "Total Population", col = "blue", lwd = 2, ylim = c(34260138, 34260140))

#calculate the root mean square deviation
rss <- function(model, target) {
  sum((model - target)^2)
}

# Fitting function
# Define the loss function

parameters_fitting <- parameters # Copy the parameters for fitting
times <- seq(0,1200)

parameters_fitting$mui <- 0.01
parameters_fitting$muo <- 0.01

parameters_fitting$beta_s <- c(1,1,1)
parameters_fitting$beta_as <- 1

#test model ode
# out <- ode(y = state, times = times, func = Malaria_model_with_Array, parms = parameters)


loss_function <- function(pars) {
  print(pars)
  
  # Fit one beta shared across beta_s and beta_as
  parameters_fitting$beta_s  <- rep(pars[1], 3)
  parameters_fitting$beta_as <- pars[1]
  
  init_state <- state
  
  out <- tryCatch(
    ode(
      y = init_state, times = times,
      func = Malaria_model_with_Array,
      parms = parameters_fitting
    ),
    error = function(e) NULL
  )
  
  if (is.null(out)) return(1e99)
  
  # 2015–2023 monthly window (your indices)
  mat_inc_s <- matrix(out[1093:1200, "inc_s"], nrow = 12)
  inc_s_model_2015_2023 <- colSums(mat_inc_s)
  
  # plot (optional; can slow optim a lot)
  plot(2015:2023, inc_s_model_2015_2023, lwd = 2, type = "l",
       ylim = c(0, 1.1 * max(Tanzania_Incidence[, 2], na.rm = TRUE)))
  points(2015:2023, Tanzania_Incidence[, 2], col = 2, lwd = 2)
  
  RSS_inc_s <- rss(inc_s_model_2015_2023, Tanzania_Incidence[, 2])
  
  if (!is.finite(RSS_inc_s)) return(1e99)
  
  print(RSS_inc_s)
  RSS_inc_s
}
# Initial guess for the parameters
initial_guess <- c(2)  # Initial guess for βs and βr
# Perform optimization
# result <- optim(
#   par = initial_guess,
#   fn = loss_function,
#   method = "L-BFGS-B",
#   lower = rep(0.001,length(initial_guess)),  # Lower bounds for βs and βr
#   upper = rep(10,length(initial_guess))
# )

result <- optim(
  par = initial_guess,
  fn = loss_function,
  hessian =T
)
# Extract the optimized parameters
optimized_params <- result$par
cat("Optimized Parameters:\n")
cat("βs:", optimized_params[1], "\n")
parameters_optim <- parameters_fitting

# Run the model with optimized parameters
parameters_optim$beta_s <- rep(optimized_params[1], 3)  # assuming same for all Is subgroups
parameters_optim$beta_as <- optimized_params[1]  # assuming same as beta_s0

saveRDS(result,"optimized_params_equilibrium.rds")

time_optimized <- seq(0, 6000, by = 1)  # Define the time points for optimized parameters

out_optimized <- ode(
  y = state,
  times = time_optimized,
  func = Malaria_model_with_Array,
  parms = parameters_optim
)

last_values <- tail(out_optimized, 1)

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

saveRDS(init_last_state,"init_state_2000.rds")

# Plot the results with optimized parameters
plot(out_optimized[5725:6000,"time"], out_optimized[5725:6000,"inc"], type = "l",
     xlab = "Time (Months)", ylab = "Incidence", main = "Malaria Incidence with Optimized Parameters",
     col = "blue", lwd = 2)

# in year
mat_inc_s_optimized <- matrix(out_optimized[5725:6000,"inc_s"], nrow = 12)
inc_s_model_2015_2023_optimized <- (colSums(mat_inc_s_optimized))
target_incidence <- Tanzania_Incidence[,2]

plot(2001:2023,inc_s_model_2015_2023_optimized,lwd=2,type = "l",ylim = c(0,1.1*max(Tanzania_Incidence[,2])),
     xlab = "Years", ylab = "Incidence", main = "Model Incidence vs Target Incidence with Optimized Parameters")
points(2015:2023,Tanzania_Incidence[,2],col=2,lwd=2)


plot(out_optimized[5725:6000,"inc_asym"], type = "l", xlab = "Months", ylab = "Incidence Sensitive",
     main = "Incidence Sensitive (Sym/asym)", col = "green", lwd = 2, ylim = c(0, max(out_optimized[5725:6000,"inc_asym"])))
lines(out_optimized[5725:6000,"inc_sym"], col = "blue", lwd = 2)

legend("bottomleft", legend = c( "Asymptomatic Incidence", "Symptomatic Incidence"),
       col = c( "green", "blue"), lty = 1, cex = 0.8)





