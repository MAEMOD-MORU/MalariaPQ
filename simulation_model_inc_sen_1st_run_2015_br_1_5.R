# Run Simulation
library(deSolve)
# setwd("D:/Work/2025/malaria/16 Jun/Baseline2_new_model/lambda_sym_asym")
#Read model
source("model_D0_33_D1_67_baseline_change_birth_death_2betaSets.r")

#
# source("function.r")

#
source("init_parameter_after_1st_run_2015_2.R")

# totalpop <- read.csv("data/Tanzania_pop_2010_2100.csv", header = TRUE)
Tanzania_data <- read.csv("data/Reported malaria cases by method of confirmation.csv")
Tanzania_Incidence <- Tanzania_data[6:14,c(1,4)]
Tanzania_k13_Allele_frequency <- read.csv("data/Tanzania_K13_Allele_frequency.csv")[,2]


init_state <- readRDS("init_state_2015_2_31Oct.rds")

# parameters$beta_s <- rep(parameters_fit[1],3)
# parameters$beta_as <- parameters_fit[2]
# parameters$beta_r <- parameters$beta_s_2 *0.625
# parameters$beta_ar <-  parameters$beta_as_2 *0.625


target_sym_asym <- rep((25/75),276) # Target incidence for symptomatic asymptomatic cases
target_s_r <- rep((1.5),12) #

times <- seq(0,12*81)
# parameters$beta_s_2 <- rep(5.7,3)
# parameters$beta_s_3 <- rep(6.2,3)
# parameters$beta_s_2 <- parameters$beta_s_2*2
parameters$start_b <- 1000

# parameters$mui <- parameters$muo

# Poputation
totalpop_2010_2100 <- read.csv("data/Tanzania_pop_2010_2100.csv", header = TRUE)

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

plot(seq(2015,2096,1/12), out[,"N"], type = "l",
     xlab = "Time (Years)", ylab = "Population",
     main = "Total Population from Model with Optimized Parameters",
     col = "blue", lwd = 2, xlim = c(2015,2050),
     ylim = c(20000000, max(out[,"N"])))
points(2010:2100, totalpop_2010_2100[,2], col = "red", pch = 19)

# 0.0000001 m = 1.1 s=1.6 as =1.4 or s=1.8 as =1.3
# 0.0000000001 m = 1.3 s=2 as=1.4
# init_mply <- 0.0000001 

events <- list(
  func = function(time, state, parameters) {
    if (time == 12 * 1) {
      # state["Ir0"] <- state["Ir0"] + state["Is0"] * init_mply
      # state["Ir1"] <- state["Ir1"] + state["Is1"] * init_mply
      # state["Ir2"] <- state["Ir2"] + state["Is2"] * init_mply
      # state["Ar"] <- state["Ar"] + state["As"] * init_mply
      state["Ir0"] <- state["Ir0"] + 1
      state["Ir1"] <- state["Ir1"] + 1
      state["Ir2"] <- state["Ir2"] + 1
      state["Ar"] <- state["Ar"] + 1
      
      # reduce Is
      # state["Is0"] <- state["Is0"] *(1-init_mply)
      # state["Is1"] <- state["Is1"] *(1-init_mply)
      # state["Is2"] <- state["Is2"] *(1-init_mply)
      # state["As"] <- state["As"] *(1-init_mply)
      state["Is0"] <- state["Is0"] -1
      state["Is1"] <- state["Is1"] -1
      state["Is2"] <- state["Is2"] -1
      state["As"] <- state["As"] -1
      
      # state["S"] <- state["S"] - (
      #   state["Is0"] * 0.05 +
      #     state["Is1"] * 0.05 +
      #     state["Is2"] * 0.05 +
      #     state["As"] * 0.05
      # )
    }
    return(state)
  },
  time = 12 * 1  # only one time point
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

parameters_optim <- parameters
# parameters_optim$beta_r <- rep(optimized_params[1], 3)  # assuming same for all Is subgroups
# parameters_optim$beta_ar <- optimized_params[2]  # assuming same as beta_s0
time_optimized <- seq(1, 12*51, by = 1)  # Define the time points for optimized parameters
library(manipulate)
manipulate(
  {
    parameters_optim$beta_s <- rep(beta_s, 3)  # assuming same for all Is subgroups
    parameters_optim$beta_as <- beta_as  # assuming same as beta_s0
    parameters_optim$beta_r <- m*rep(beta_s, 3)  # assuming same for all Is subgroups
    parameters_optim$beta_ar <- m*beta_as  # assuming same as beta_s0
    
    init_state_test <- init_state*i

    events <- list(
      func = function(time, state, parameters) {
        if (time == time_start_res_after*12) {
          # state["Ir0"] <- state["Ir0"] + state["Is0"] * init_mply
          # state["Ir1"] <- state["Ir1"] + state["Is1"] * init_mply
          # state["Ir2"] <- state["Ir2"] + state["Is2"] * init_mply
          # state["Ar"] <- state["Ar"] + state["As"] * init_mply
          state["Ir0"] <- state["Ir0"] + 1
          state["Ir1"] <- state["Ir1"] + 1
          state["Ir2"] <- state["Ir2"] + 1
          state["Ar"] <- state["Ar"] + 1
          
          # reduce Is
          # state["Is0"] <- state["Is0"] *(1-init_mply)
          # state["Is1"] <- state["Is1"] *(1-init_mply)
          # state["Is2"] <- state["Is2"] *(1-init_mply)
          # state["As"] <- state["As"] *(1-init_mply)
          state["Is0"] <- state["Is0"] -1
          state["Is1"] <- state["Is1"] -1
          state["Is2"] <- state["Is2"] -1
          state["As"] <- state["As"] -1
          
          # state["S"] <- state["S"] - (
          #   state["Is0"] * 0.05 +
          #     state["Is1"] * 0.05 +
          #     state["Is2"] * 0.05 +
          #     state["As"] * 0.05
          # )
        }
        return(state)
      },
      time = 12 * time_start_res_after  # only one time point
    )
    
    
    out <- ode(y = init_state_test, times = time_optimized,
               func = Malaria_model_with_Array, parms = parameters_optim
               ,events = events
               )
    
    mat_inc <- matrix(out[,"inc"], nrow = 12)
    inc_model_2015_2023 <- (colSums(mat_inc))
    
    mat_inc <- matrix(out[,"inc_r"], nrow = 12)
    inc_model_2015_2023_r <- (colSums(mat_inc))
    
    mat_inc <- matrix(out[,"N"], nrow = 12)
    inc_model_2015_2023_n <- (colMeans(mat_inc))
    
    mat_inc <- matrix(out[,"Is0"], nrow = 12)
    inc_model_2015_2023_GIs0 <- (colMeans(mat_inc))
    
    mat_inc <- matrix(out[,"Is1"], nrow = 12)
    inc_model_2015_2023_GIs1 <- (colMeans(mat_inc))
    
    mat_inc <- matrix(out[,"Is2"], nrow = 12)
    inc_model_2015_2023_GIs2 <- (colMeans(mat_inc))
    
    mat_ratio <- matrix(out[,"inc_r"]/(out[,"inc_r"]+out[,"inc_s"]), nrow = 12)
    ratio_2015_2023 <- (colMeans(mat_ratio))
    
    par(mfrow =c(2,1))
    plot(2000:2050,xlab="year",
         inc_model_2015_2023,lwd=2,type = "l",main="total inc",
         xlim=c(year_plot,2050)
         )
    lines(2000:2050,inc_model_2015_2023_r,col="red")
    lines(2000:2050,inc_model_2015_2023_n,col="green3")
    points(2015:2023,Tanzania_Incidence[,2],col="blue",lwd=2)
    # Add secondary axis
    # par(new = TRUE)
    plot(2000:2050,ratio_2015_2023, type = "l", col = "red", xlab = "ratio inc_r/total_inc", ylab = "",
         xlim=c(year_plot,2050)
         )
    points(Tanzania_k13_Allele_frequency,col=2,lwd=2,pch=3)
    points(2025,0.4,col=2,lwd=2,pch=3)
    # axis(side = 4)
    # plot(inc_model_2015_2023_GIs0,type = "l",col="red")
    # lines(inc_model_2015_2023_GIs1,col="green")
    # lines(inc_model_2015_2023_GIs2,col="blue")
    
  },
  beta_s = slider(0,2,step=0.00001,initial =0.68),
  beta_as = slider(0,10,step=0.00001,initial =6.8),
  m = slider(0.1,3,step=0.01,initial =1.05),
  i = slider(0.00001,2,step=0.00001,initial =0.08),
  year_plot = slider(2000,2035,step=1,initial =2000),
  time_start_res_after = slider(0,25,step=1,initial =16)
)
