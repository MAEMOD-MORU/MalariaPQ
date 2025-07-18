# Run Simulation
library(deSolve)

#Read model
source("model_D0_33_D1_67_baseline_change_birth_death.r")

#
source("function.r")

#
source("init_parameter_sen_only.r")

# Poputation
totalpop_2000_2100 <- read.csv("data/uganda_totalpop_2000_2100.csv", header = TRUE)

#calculate the root mean square deviation
rmsd <- function(model,data){
  diff <- (model - data)^2
  mean_sq <- sum(diff)/length(data)
  return(sqrt(mean_sq))
}

parameters_fitting <- parameters
times_fit <- seq(0, 12*51, 1) # 12 months per year, 18 years

#the objective function for using with the 'optim' command
pop_run <- function(pars){
  # update the parameters
  parameters_fitting$mui <- pars[1]
  parameters_fitting$muo <- pars[2]
  parameters_fitting$mui_2025 <- pars[3]
  parameters_fitting$muo_2025 <- pars[4]
  # Event: change beta at t = 120
  
  
  out <- ode(y = state, times = times_fit,
             func = Malaria_model_with_Array, 
             parms = parameters_fitting
             )
  # get the total population from the model output
  pop_model <- out[,"N"]
  # total population from the model
  # 2000 - 2050
  totalpop_model <- pop_model[seq(0, 12*51, 12)]
  #line + points
  plot(2000:2050, totalpop_model, type = "o", col = "blue", 
       xlab = "Time (Years)", ylab = "Population",
       main = "Total Population from Model", ylim = c(20000000, 86000000))
  points(2000:2050,totalpop_2000_2100[1:51,2],col="red", pch = 19)
  
  
  # calculate the root mean square deviation
  rmsd_value <- rmsd(totalpop_model, totalpop_2000_2100[1:51,2])
  print(paste("RMSD:", rmsd_value))
  if (is.na(rmsd_value) || rmsd_value == Inf || rmsd_value == -Inf) {
    return(1e10)  # Return a large value if RMSD is not valid
  }
  return(rmsd_value)
}

# Initial guess for the parameters
initial_guess <- c(0.0125, 0.01,0.012,0.01)  # Initial guess for mui and muo
# Perform optimization
result <- optim(par = initial_guess, fn = pop_run, method = "L-BFGS-B",
                 lower = c(0.0001, 0.0001,0.0001,0.0001), 
                upper = c(0.5, 0.5,0.5,0.5),hessian =T)
result
# Extract the optimized parameters
optimized_params <- result$par
cat("Optimized Parameters:\n")
cat("mui:", optimized_params[1], "\n")
cat("muo:", optimized_params[2], "\n")
cat("mui_2025:", optimized_params[3], "\n")
cat("muo_2025:", optimized_params[4], "\n")

parameters_optimized <- parameters_fitting
# Run the model with optimized parameters
parameters_optimized$mui <- optimized_params[1]
parameters_optimized$muo <- optimized_params[2]
parameters_optimized$mui_2025 <- optimized_params[3]
parameters_optimized$muo_2025 <- optimized_params[4]

out_optimized <- ode(y = state, times = times_fit,
                      func = Malaria_model_with_Array, parms = parameters_optimized)

# Plot the results with optimized parameters
plot(seq(2000,2051,1/12), out_optimized[,"N"], type = "l",
     xlab = "Time (Years)", ylab = "Population",
     main = "Total Population from Model with Optimized Parameters",
     col = "blue", lwd = 2, xlim = c(2000,2050),
     ylim = c(20000000, max(out_optimized[1:600,"N"])))
points(2000:2100, totalpop_2000_2100[,2], col = "red", pch = 19)
# Add legend
legend("topleft", legend = c("Optimized Population", "Observed Population"),
       col = c("blue", "red"), pch = c(NA, 19), lwd = 2, cex = 0.8)

#Extract the Hessian matrix
hessian_matrix <- result$hessian
#Compute the variance-covariance matrix
variance_covariance_matrix <- solve(hessian_matrix)
#Compute standard errors (SEs):
standard_errors <- sqrt(diag(variance_covariance_matrix))
#Calculate the 95% confidence intervals
estimates <- result$par
lower <- estimates - 1.96 * se
upper <- estimates + 1.96 * se

ci_95 <- data.frame(Estimate = estimates,
                    SE = se,
                    Lower95 = lower,
                    Upper95 = upper)

#Use the generalized inverse instead
library(MASS)
vcov <- ginv(result$hessian)
se <- sqrt(diag(vcov))  # may still have imaginary values
# Calculate 95% confidence intervals
lower <- estimates - 1.96 * se
upper <- estimates + 1.96 * se
ci_95 <- data.frame(Estimate = estimates,
                    SE = se,
                    Lower95 = lower,
                    Upper95 = upper)

#plot ci95
parametert_lower <- parameters_fitting
parametert_lower$mui <- ci_95$Lower95[1]
parametert_lower$muo <- ci_95$Lower95[2]
out_lower <- ode(y = state, times = times_fit,
                      func = Malaria_model_with_Array, parms = parametert_lower)
parametert_upper <- parameters_fitting
parametert_upper$mui <- ci_95$Upper95[1]
parametert_upper$muo <- ci_95$Upper95[2]
out_upper <- ode(y = state, times = times_fit,
                      func = Malaria_model_with_Array, parms = parametert_upper)
# Plot the results with confidence intervals
plot(out_optimized[,"time"]/12, out_optimized[,"N"], type = "l",
     xlab = "Time (Years)", ylab = "Population",
     main = "Total Population from Model with Confidence Intervals",
     col = "blue", lwd = 2, ylim = c(20000000, 55000000))
lines(out_lower[,"time"]/12, out_lower[,"N"], col = "green", lty = 2)
lines(out_upper[,"time"]/12, out_upper[,"N"], col = "red", lty = 2)
#polygon
polygon(c(out_lower[,"time"]/12, rev(out_upper[,"time"]/12)), 
         c(out_lower[,"N"], rev(out_upper[,"N"])), 
         col = rgb(0, 1, 0, 0.2), border = NA)
points(0:24, totalpop[,2], col = "red", pch = 19)
# Add legend
legend("topright", legend = c("Optimized Population", "Lower CI", "Upper CI", "Observed Population"),
       col = c("blue", "green", "red", "red"), pch = c(NA, NA, NA, 19), lwd = 2, lty = c(1, 2, 2, NA), cex = 0.8)
