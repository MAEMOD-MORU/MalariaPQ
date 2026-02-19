#detach all library
rm(list = ls())

source("init_parameter_fitted.R")
source("init_parameter_SSP1.R")
source("init_parameter_SSP2.R")
source("init_parameter_SSP3.R")
source("model_D0_33_D1_67_baseline_change_birth_death_2betaSets_SS.R")

library(deSolve)

init_state <- readRDS("init_state_d1_33.rds")
totalpop_data <- read.csv("data/Tanzania_pop_2000_2100.csv", header = TRUE)[1:51,2]
data_ssp <- read.csv("data/Tanzania_pop_SSP_2000_2100.csv", header = TRUE)[1:11,2:6]
ssp_year <- seq(2000, 2050, by = 5)
inc_data <- read.csv("data/Reported malaria cases by method of confirmation.csv", header = TRUE)[6:14,4]
k13_data <- read.csv("data/Tanzania_K13_Allele_frequency.csv", header = TRUE)[2:8,2]

events <- list(
  func = function(time, state, parameters) {
    if (time == 12 * 1) {
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
  time = 12 * 1  # only one time point
)

times <- seq(1, 50*12, by = 1)
parameters_fitted <- parameters 
out_model <- ode(y = init_state, times = times, events = events,
                 func = Malaria_model_with_Array, parms = parameters_fitted)
out_model_pop <- out_model[,"N"]
out_model_inc <- colSums(matrix(out_model[, "inc"],   nrow = 12))
out_model_ratio_inc_r <- colMeans(matrix(out_model[, "inc_r"] / out_model[, "inc"], nrow = 12))

# Run the model for SSP1
out_ssp1 <- ode(y = init_state, times = times, events = events,
                func = Malaria_model_with_Array, parms = parameters_ssp1)
out_ssp1_pop <- out_ssp1[,"N"] 
out_ssp1_inc <- colSums(matrix(out_ssp1[,"inc"],   nrow = 12))
out_ssp1_ratio_inc_r <- colMeans(matrix(out_ssp1[,"inc_r"]/ out_ssp1[,"inc"], nrow = 12))

# Run the model for SSP2
out_ssp2 <- ode(y = init_state, times = times, events = events,
                func = Malaria_model_with_Array, parms = parameters_ssp2)
out_ssp2_pop <- out_ssp2[,"N"] 
out_ssp2_inc <- colSums(matrix(out_ssp2[,"inc"],   nrow = 12))
out_ssp2_ratio_inc_r <- colMeans(matrix(out_ssp2[,"inc_r"]/ out_ssp2[,"inc"], nrow = 12))

# Run the model for SSP3
out_ssp3 <- ode(y = init_state, times = times, events = events,
                func = Malaria_model_with_Array, parms = parameters_ssp3)
out_ssp3_pop <- out_ssp3[,"N"] 
out_ssp3_inc <- colSums(matrix(out_ssp3[,"inc"],   nrow = 12))
out_ssp3_ratio_inc_r <- colMeans(matrix(out_ssp3[,"inc_r"]/ out_ssp3[,"inc"], nrow = 12))

year_m <- seq(0, 599/12, by = 1/12)+2000
year <- 2000:2050

# Plot the results
plot(year_m, out_ssp1_pop, type = "l", col = "blue",
     xlab = "Time (Year)", ylab = "Population", 
     main = "Total Population for Different SSPs",
     ylim = range(c(0,out_ssp1_pop, out_ssp2_pop, out_ssp3_pop)))
lines(year_m, out_model_pop, col = "black",lwd=1)
lines(year_m, out_ssp2_pop, col = "red")
lines(year_m, out_ssp3_pop, col = "green3")
points(2000:2050,totalpop_data, col = "black", pch = 16)
points(ssp_year, data_ssp$SSP1, col = "blue", pch = 18)
points(ssp_year, data_ssp$SSP2, col = "red", pch = 15)
points(ssp_year, data_ssp$SSP3, col =" green", pch = 17)
grid()
legend("topleft", legend = c("SSP1 Model", "SSP2 Model", "SSP3 Model", "Baseline Model",
                            "SSP1 Data", "SSP2 Data", "SSP3 Data","UN medium-fertility"),ncol = 2,
       col = c("blue", "red", "green3", "black", "blue", "red", "green3","black"),
       pch = c(NA, NA, NA, NA, 18, 15, 17,16), lty = c(1, 1, 1, 1, NA, NA, NA,NA),lwd=3)

# plot inc
plot(year, out_ssp1_inc, type = "l", col = "blue",
     xlab = "Time (months)", ylab = "Incidence", 
     main = "Incidence over Time for Different SSPs",
     ylim = range(c(0,out_ssp1_inc, out_ssp2_inc, out_ssp3_inc)))
lines(year, out_model_inc, col = "black",lwd=1)
lines(year, out_ssp2_inc, col = "red")
lines(year, out_ssp3_inc, col = "green")
points(2015:2023,inc_data, col = "black", pch = 16)
legend("topleft", legend = c("SSP1 Model", "SSP2 Model ", "SSP3 Model", "Baseline Model"),
       col = c("blue", "red", "green", "black"),
       lty = c(1, 1, 1, 1))

#plot ratio inc r
plot(year, out_ssp1_ratio_inc_r, type = "l", col = "blue",
     xlab = "Time (months)", ylab = "Ratio Incidence Resistance", 
     main = "Ratio Incidence Resistance over Time for Different SSPs",
     ylim = range(c(0,out_ssp1_ratio_inc_r, out_ssp2_ratio_inc_r, out_ssp3_ratio_inc_r)))
lines(year, out_model_ratio_inc_r, col = "black",lwd=1)
lines(year, out_ssp2_ratio_inc_r, col = "red")
lines(year, out_ssp3_ratio_inc_r, col = "green")
points(2016:2022,k13_data, col = "red", pch = 16)
points(2025,0.4, col = "red", pch = 16)
#line on 2025
abline(v = 2025, col = "darkgreen", lty = 3)
legend("topleft", legend = c("SSP1 Model", "SSP2 Model ", "SSP3 Model", "Baseline Model"),
       col = c("blue", "red", "green", "black"),
       lty = c(1, 1, 1, 1))

