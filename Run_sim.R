# Run Simulation
library(deSolve)

#Read model
source("model.r")

#
source("function.r")

# Initial Population
initP <- 70000000

# Initial States from UI Inputs
init_total_IAs <- 10
init_total_IAr <- 10

prob_s = 0.1 # Ratio of infection on sym/asym sensitive group
prob_r = 0.1 # Ratio of infection on sym/asym resistance group

Is_total <- init_total_IAs * prob_s
As_total <- init_total_IAs * (1-prob_s)

Ir_total <- init_total_IAr * prob_r
Ar_total <- init_total_IAr * (1-prob_r)

init_D1s <- 0.33
init_D2s <- 0.33

init_D1r <- 0.33
init_D2r <- 0.33

propIs <- c(1 - (init_D1s + init_D2s), init_D1s, init_D2s) # ratio of
propIr <- c(1 - (init_D1r + init_D2r), init_D1r, init_D2r)

initIs <- Is_total * propIs  # Sensitive subgroups
initAs <- As_total
initGIs <- c(0, 0, 0)
initGAs <- 0
initRs <- 0

initIr <- Ir_total * propIr  # Resistant subgroups
initAr <- Ar_total
initGIr <- c(0, 0, 0)
initGAr <- 0
initRr <- 0

# Calculate initial susceptible population
initS <- initP - sum(initIs) - sum(initGIs) - initGAs - initAs - initRs -
  sum(initIr) - sum(initGIr) - initGAr - initAr - initRr

# Dynamic tau from UI
tau <- c(
  0.01, 0.01, 0.01, 0.01,
  0.01, 0.01, 0.01, 0.01
)* 30

beta_s0 = 0.05* 30
beta_s1 = 0.05* 30
beta_s2 = 0.05* 30
beta_r0 = 0.05* 30
beta_r1 = 0.05* 30
beta_r2 = 0.05* 30

parameters <- list(
  mui = (1 / 50) / 12,
  muo = (1 / 50) / 12,
  flo = 1,
  amp = 0.5,
  period = 12,
  phase = 0,
  start_d = 3620,# 10 years after 3500 months

  beta_s = c(beta_s0, beta_s1, beta_s2),
  beta_r = c(beta_r0, beta_r1,beta_r2),
  beta_as = beta_s0,
  beta_ar = beta_r0,

  alpha = 1.07 / 12, #  doi:10.1371/journal.pone.0001767
  tau = tau,

  d0 = (1 / 44) * 30, # Gametocyte recovery rate from Palang paper, table 2 at row Placebo (DP), column microscopy
  d1 = (1 / 20) * 30, # Can't find the ACT.
  d2 = (1 / 7.88) * 30, # Gametocyte recovery rate from Palang paper, table 2 at row 0.0625 mg/kg, column microscopy

  gamma_i = 0.02* 30,
  gamma_a = 0.02* 30,
  prob_s = prob_s,
  propIs = propIs,
  prob_r = prob_r,
  propIr = propIr
)

state <- c(
  S = initS,
  Is0 = initIs[1], Is1 = initIs[2], Is2 = initIs[3], # Flatten Is array into individual components
  GIs0 = initGIs[1], GIs1 = initGIs[2], GIs2 = initGIs[3],  # Gametocyte from Is
  GAs = initGAs,  # Gametocyte from As
  As = initAs,
  Rs = initRs,
  Ir0 = initIr[1], Ir1 = initIr[2], Ir2 = initIr[3], # Flatten Ir array into individual components
  GIr0 = initGIr[1], GIr1 = initGIr[2], GIr2 = initGIr[3],  # Gametocyte from Ir
  GAr = initGAr,  # Gametocyte from Ar
  Ar = initAr,
  Rr = initRr
)

times <- seq(0, 400*12, by = 1)
out<- ode(y = state, times = times, func = Malaria_model_with_Array, parms = parameters)

plot(out[1:2000,"inc"],type="l",xlab="Months",ylab="Incidence")

# incidence
plot(x=0:240,out[3500:3740,"inc"],type="l",xlab="months",ylab="incidence")

# gametocyte sensitive
plot(x=0:240, out[3500:3740,c("GIs0")],type="l",xlab="months",ylab="gametocyte sensitive")
lines(out[3500:3740,c("GIs1")],col=2)
lines(out[3500:3740,c("GIs2")],col=3)
title("Sensitive")

# gametocyte res
plot(x=0:240, out[3500:3740,c("GIr0")],type="l",ylim=c(3e5,6e5),xlab="months",ylab="gametocyte resistance")
lines(out[3500:3740,c("GIr1")],col=2)
lines(out[3500:3740,c("GIr2")],col=3)
title("Resistance")




#matrix
n_run <- 10

matrix_output_inc <- matrix(nrow = length(times),ncol = n_run)

matrix_output_GIs0 <- matrix(nrow = length(times),ncol = n_run)
matrix_output_GIs1 <- matrix(nrow = length(times),ncol = n_run)
matrix_output_GIs2 <- matrix(nrow = length(times),ncol = n_run)

matrix_output_GIr0 <- matrix(nrow = length(times),ncol = n_run)
matrix_output_GIr1 <- matrix(nrow = length(times),ncol = n_run)
matrix_output_GIr2 <- matrix(nrow = length(times),ncol = n_run)

matrix_output_total_inf <- matrix(nrow = length(times),ncol = n_run)
matrix_output_total_G <- matrix(nrow = length(times),ncol = n_run)

matrix_output_Gs <- matrix(nrow = length(times),ncol = n_run)
matrix_output_Gr <- matrix(nrow = length(times),ncol = n_run)

matrix_output_Gsym <- matrix(nrow = length(times),ncol = n_run)
matrix_output_Gasym <- matrix(nrow = length(times),ncol = n_run)

pb <- txtProgressBar(max =n_run)
for (i in 1:n_run) {
  # d0 ~ Uniform(10.9–75.6)
  # d1 ~ Uniform(2.21–60.0)
  # d2 ~ Uniform(1.37–11.0)
  rd0 <- runif(n=1,min = 10.9,max = 44) # 44.1
  rd1 <- runif(n=1,min = 2.21,max = 20) #7.88
  rd2 <- runif(n=1,min = 1.37,max = 11.0) #4.88
  parameters$d0 <- 1/rd0*30
  parameters$d1 <- 1/rd1*30
  parameters$d2 <- 1/rd2*30

  out<- ode(y = state, times = times, func = Malaria_model_with_Array, parms = parameters)
  matrix_output_inc[,i] <- out[,"inc"]

  matrix_output_GIs0[,i] <- out[,"GIs0"]
  matrix_output_GIs1[,i] <- out[,"GIs1"]
  matrix_output_GIs2[,i] <- out[,"GIs2"]

  matrix_output_GIr0[,i] <- out[,"GIr0"]
  matrix_output_GIr1[,i] <- out[,"GIr1"]
  matrix_output_GIr2[,i] <- out[,"GIr2"]
  
  matrix_output_total_inf[,i] <- out[,"total_inf"]
  matrix_output_total_G[,i] <- out[,"total_G"]
  
  matrix_output_Gs[,i] <- out[,"Gs"]
  matrix_output_Gr[,i] <- out[,"Gr"]
  
  matrix_output_Gsym[,i] <- out[,"Gsym"]
  matrix_output_Gasym[,i] <- out[,"Gasym"]
  setTxtProgressBar(pb, i)
  }

matrix_output <- run_matrix_output(100,360)


plot_out_line(matrix_output_inc[3500:3740,])
plot_out_line(matrix_output_GIs0[3500:3740,],ylim=c(0,1e6))
plot_out_line(matrix_output_GIs1[3500:3740,],ylim=c(0,1e6))
plot_out_line(matrix_output_GIs2[3500:3740,],ylim=c(0,1e6))

plot_out_CI95(matrix_output_inc[3500:3740,])
plot_out_CI95(matrix_output_GIs0[3500:3740,],c(ylim=c(0,1e6),xlim=c(1,10)))
plot_out_CI95(matrix_output_GIs1[3500:3740,],ylim=c(0,1e6))
plot_out_CI95(matrix_output_GIs2[3500:3740,],ylim=c(0,1e6))


#### inc ######
# Compute median, lower (2.5th), and upper (97.5th) percentiles
median_inc <- apply(matrix_output_inc[3500:3740,], 1, median, na.rm = TRUE)
lower_inc <- apply(matrix_output_inc[3500:3740,], 1, quantile, probs = 0.025, na.rm = TRUE)
upper_inc <- apply(matrix_output_inc[3500:3740,], 1, quantile, probs = 0.975, na.rm = TRUE)

# Plot with confidence interval
plot(1:nrow(matrix_output_inc[3500:3740,]), median_inc, type = "l", col = "blue", lwd = 1.5, ylim = range(lower_inc, upper_inc),
     ylab = "Incidence", xlab = "Months", main = "Total Incidence with 95% CI")

# Add shaded region for confidence interval
polygon(c(1:nrow(matrix_output_inc[3500:3740,]), rev(1:nrow(matrix_output_inc[3500:3740,]))), c(lower_inc, rev(upper_inc)), col = rgb(0.5, 0.5, 0.5, 0.75), border = NA)
lines(1:nrow(matrix_output_inc[3500:3740,]),lower_inc, col="grey50",lty=2)
lines(1:nrow(matrix_output_inc[3500:3740,]),upper_inc, col="grey50",lty=2)

# Add median line
lines(1:nrow(matrix_output_inc[3500:3740,]), median_inc, col = "blue", lwd = 1.5)


#### gametocyte sensitive ######
median_GIs0 <- apply(matrix_output_GIs0[3500:3740,], 1, median, na.rm = TRUE)
lower_GIs0 <- apply(matrix_output_GIs0[3500:3740,], 1, quantile, probs = 0.025, na.rm = TRUE)
upper_GIs0 <- apply(matrix_output_GIs0[3500:3740,], 1, quantile, probs = 0.975, na.rm = TRUE)

median_GIs1 <- apply(matrix_output_GIs1[3500:3740,], 1, median, na.rm = TRUE)
lower_GIs1 <- apply(matrix_output_GIs1[3500:3740,], 1, quantile, probs = 0.025, na.rm = TRUE)
upper_GIs1 <- apply(matrix_output_GIs1[3500:3740,], 1, quantile, probs = 0.975, na.rm = TRUE)

median_GIs2 <- apply(matrix_output_GIs2[3500:3740,], 1, median, na.rm = TRUE)
lower_GIs2 <- apply(matrix_output_GIs2[3500:3740,], 1, quantile, probs = 0.025, na.rm = TRUE)
upper_GIs2 <- apply(matrix_output_GIs2[3500:3740,], 1, quantile, probs = 0.975, na.rm = TRUE)


# Plot with confidence interval
plot(1:nrow(matrix_output_GIs0[3500:3740,]), median_GIs0, type = "l", col = "green", lwd = 1.5, ylim = range(matrix_output_GIs0, matrix_output_GIs0),
     ylab = "Number of People Carrying Gametocytes", xlab = "Months", main = "Sensitive with 95% CI")

# Add shaded region for confidence interval
polygon(c(1:nrow(matrix_output_GIs0[3500:3740,]), rev(1:nrow(matrix_output_GIs0[3500:3740,]))), c(lower_GIs0, rev(upper_GIs0)), col = rgb(0.2, 0.5, 0.2, 0.25), border = NA)
lines(1:nrow(matrix_output_GIs0[3500:3740,]),lower_GIs0, col="green",lty=2)
lines(1:nrow(matrix_output_GIs0[3500:3740,]),upper_GIs0, col="green",lty=2)

polygon(c(1:nrow(matrix_output_GIs1[3500:3740,]), rev(1:nrow(matrix_output_GIs1[3500:3740,]))), c(lower_GIs1, rev(upper_GIs1)), col = rgb(0.5, 0.2, 0.2, 0.25), border = NA)
lines(1:nrow(matrix_output_GIs1[3500:3740,]),lower_GIs1, col="red",lty=2)
lines(1:nrow(matrix_output_GIs1[3500:3740,]),upper_GIs1, col="red",lty=2)

polygon(c(1:nrow(matrix_output_GIs2[3500:3740,]), rev(1:nrow(matrix_output_inc[3500:3740,]))), c(lower_GIs2, rev(upper_GIs2)), col = rgb(0.2, 0.2, 0.5, 0.25), border = NA)
lines(1:nrow(matrix_output_GIs2[3500:3740,]),lower_GIs2, col="blue",lty=2)
lines(1:nrow(matrix_output_GIs2[3500:3740,]),upper_GIs2, col="blue",lty=2)

# Add median line
lines(1:nrow(matrix_output_GIs0[3500:3740,]), median_GIs0, col = rgb(0.1, 0.7, 0.1, 1), lwd = 2)
lines(1:nrow(matrix_output_GIs1[3500:3740,]), median_GIs1, col = rgb(0.7, 0.1, 0.1, 1), lwd = 2)
lines(1:nrow(matrix_output_GIs2[3500:3740,]), median_GIs2, col = rgb(0.1, 0.1, 0.7, 1), lwd = 2)
legend("bottomleft", legend = c("Placebo", "0.0625 mg/kg", "0.25 mg/kg"),
       col = c(rgb(0.1, 0.7, 0.1, 1), rgb(0.7, 0.1, 0.1, 1), rgb(0.1, 0.1, 0.7, 1)),
       lwd = 2, lty = 1, bg = "white",cex = 1.5)


#### gametocyte resistance ######
median_GIr0 <- apply(matrix_output_GIr0[3500:3740,], 1, median, na.rm = TRUE)
lower_GIr0 <- apply(matrix_output_GIr0[3500:3740,], 1, quantile, probs = 0.025, na.rm = TRUE)
upper_GIr0 <- apply(matrix_output_GIr0[3500:3740,], 1, quantile, probs = 0.975, na.rm = TRUE)

median_GIr1 <- apply(matrix_output_GIr1[3500:3740,], 1, median, na.rm = TRUE)
lower_GIr1 <- apply(matrix_output_GIr1[3500:3740,], 1, quantile, probs = 0.025, na.rm = TRUE)
upper_GIr1 <- apply(matrix_output_GIr1[3500:3740,], 1, quantile, probs = 0.975, na.rm = TRUE)

median_GIr2 <- apply(matrix_output_GIr2[3500:3740,], 1, median, na.rm = TRUE)
lower_GIr2 <- apply(matrix_output_GIr2[3500:3740,], 1, quantile, probs = 0.025, na.rm = TRUE)
upper_GIr2 <- apply(matrix_output_GIr2[3500:3740,], 1, quantile, probs = 0.975, na.rm = TRUE)

# Plot with confidence interval
plot(1:nrow(matrix_output_GIr0[3500:3740,]), median_GIr0, type = "l", col = "green", lwd = 1.5, ylim = range(matrix_output_GIr0, matrix_output_GIr0),
     ylab = "Number of People Carrying Gametocytes", xlab = "Months", main = "Resistance with 95% CI")


# Add shaded region for confidence interval
polygon(c(1:nrow(matrix_output_GIr0[3500:3740,]), rev(1:nrow(matrix_output_GIr0[3500:3740,]))), c(lower_GIr0, rev(upper_GIr0)), col = rgb(0.2, 0.5, 0.2, 0.25), border = NA)
lines(1:nrow(matrix_output_GIr0[3500:3740,]),lower_GIr0, col="green",lty=2)
lines(1:nrow(matrix_output_GIr0[3500:3740,]),upper_GIr0, col="green",lty=2)

polygon(c(1:nrow(matrix_output_GIr1[3500:3740,]), rev(1:nrow(matrix_output_GIr1[3500:3740,]))), c(lower_GIr1, rev(upper_GIr1)), col = rgb(0.5, 0.2, 0.2, 0.25), border = NA)
lines(1:nrow(matrix_output_GIr1[3500:3740,]),lower_GIr1, col="red",lty=2)
lines(1:nrow(matrix_output_GIr1[3500:3740,]),upper_GIr1, col="red",lty=2)

polygon(c(1:nrow(matrix_output_GIr2[3500:3740,]), rev(1:nrow(matrix_output_inc[3500:3740,]))), c(lower_GIr2, rev(upper_GIr2)), col = rgb(0.2, 0.2, 0.5, 0.25), border = NA)
lines(1:nrow(matrix_output_GIr2[3500:3740,]),lower_GIr2, col="blue",lty=2)
lines(1:nrow(matrix_output_GIr2[3500:3740,]),upper_GIr2, col="blue",lty=2)

# Add median line
lines(1:nrow(matrix_output_GIr0[3500:3740,]), median_GIr0, col =rgb(0.1, 0.7, 0.1, 1), lwd = 2)
lines(1:nrow(matrix_output_GIr1[3500:3740,]), median_GIr1, col = rgb(0.7, 0.1, 0.1, 1), lwd = 2)
lines(1:nrow(matrix_output_GIr2[3500:3740,]), median_GIr2, col = rgb(0.1, 0.1, 0.7, 1), lwd = 2)
legend("bottomleft", legend = c("Placebo", "0.0625 mg/kg", "0.25 mg/kg"),
       col = c(rgb(0.1, 0.7, 0.1, 1), rgb(0.7, 0.1, 0.1, 1), rgb(0.1, 0.1, 0.7, 1)),
       lwd = 2, lty = 1, bg = "white",cex = 1.5)

#### total g/total inf(IS+AS+IR+AR) * 100% ######
median_G_per_inf <- apply((matrix_output_total_G[3500:3740,]/matrix_output_total_inf[3500:3740,]), 1, median, na.rm = TRUE)
lower_G_per_inf  <- apply((matrix_output_total_G[3500:3740,]/matrix_output_total_inf[3500:3740,]), 1, quantile, probs = 0.025, na.rm = TRUE)
upper_G_per_inf  <- apply((matrix_output_total_G[3500:3740,]/matrix_output_total_inf[3500:3740,]), 1, quantile, probs = 0.975, na.rm = TRUE)


# Plot with confidence interval
plot(1:nrow(matrix_output_inc[3500:3740,]), median_G_per_inf, type = "l", col = "yellow", lwd = 2, ylim = range(lower_G_per_inf, upper_G_per_inf),
     ylab = "Ratio", xlab = "Months", main = "Total Carrying Gametocytes/total Infected")

# Add shaded region for confidence interval
polygon(c(1:nrow(matrix_output_inc[3500:3740,]), rev(1:nrow(matrix_output_inc[3500:3740,]))), c(lower_G_per_inf, rev(upper_G_per_inf)), col = rgb(0.75, 0.75, 0.75, 0.5), border = NA)
lines(1:nrow(matrix_output_inc[3500:3740,]),lower_G_per_inf, col="grey75",lty=2)
lines(1:nrow(matrix_output_inc[3500:3740,]),upper_G_per_inf, col="grey75",lty=2)

# Add median line
lines(1:nrow(matrix_output_inc[3500:3740,]), median_G_per_inf, col = "yellow", lwd = 2)
#### total gr/total gs ######
median_Gr_per_Gs <- apply((matrix_output_Gs[3500:3740,]/matrix_output_Gr[3500:3740,]), 1, median, na.rm = TRUE)
lower_Gr_per_Gs  <- apply((matrix_output_Gs[3500:3740,]/matrix_output_Gr[3500:3740,]), 1, quantile, probs = 0.025, na.rm = TRUE)
upper_Gr_per_Gs  <- apply((matrix_output_Gs[3500:3740,]/matrix_output_Gr[3500:3740,]), 1, quantile, probs = 0.975, na.rm = TRUE)


# Plot with confidence interval
plot(1:nrow(matrix_output_inc[3500:3740,]), median_Gr_per_Gs, type = "l", col = "orange", lwd = 2, ylim = range(lower_Gr_per_Gs, upper_Gr_per_Gs),
     ylab = "Ratio", xlab = "Months", main = "Total Carrying Gametocytes/total Infected")

# Add shaded region for confidence interval
polygon(c(1:nrow(matrix_output_inc[3500:3740,]), rev(1:nrow(matrix_output_inc[3500:3740,]))), c(lower_Gr_per_Gs, rev(upper_Gr_per_Gs)), col = rgb(0.75, 0.75, 0.75, 0.5), border = NA)
lines(1:nrow(matrix_output_inc[3500:3740,]),lower_Gr_per_Gs, col="grey75",lty=2)
lines(1:nrow(matrix_output_inc[3500:3740,]),upper_Gr_per_Gs, col="grey75",lty=2)

# Add median line
lines(1:nrow(matrix_output_inc[3500:3740,]), median_Gr_per_Gs, col = "orange", lwd = 2)

#### total gsym/total gasym ######

median_Gsym_per_Gasym <- apply((matrix_output_Gsym[3500:3740,]/matrix_output_Gasym[3500:3740,]), 1, median, na.rm = TRUE)
lower_Gsym_per_Gasym  <- apply((matrix_output_Gsym[3500:3740,]/matrix_output_Gasym[3500:3740,]), 1, quantile, probs = 0.025, na.rm = TRUE)
upper_Gsym_per_Gasym  <- apply((matrix_output_Gsym[3500:3740,]/matrix_output_Gasym[3500:3740,]), 1, quantile, probs = 0.975, na.rm = TRUE)


# Plot with confidence interval
plot(1:nrow(matrix_output_inc[3500:3740,]), median_Gsym_per_Gasym, type = "l", col = "purple2", lwd = 2, ylim = range(lower_Gsym_per_Gasym, upper_Gsym_per_Gasym),
     ylab = "Ratio", xlab = "Months", main = "Total Symptomatic Carrying Gametocytes/Total Asymptomatic Carrying Gametocytes")

# Add shaded region for confidence interval
polygon(c(1:nrow(matrix_output_inc[3500:3740,]), rev(1:nrow(matrix_output_inc[3500:3740,]))), c(lower_Gsym_per_Gasym, rev(upper_Gsym_per_Gasym)), col = rgb(0.75, 0.75, 0.75, 0.5), border = NA)
lines(1:nrow(matrix_output_inc[3500:3740,]),lower_Gsym_per_Gasym, col="grey75",lty=2)
lines(1:nrow(matrix_output_inc[3500:3740,]),upper_Gsym_per_Gasym, col="grey75",lty=2)

# Add median line
lines(1:nrow(matrix_output_inc[3500:3740,]), median_Gsym_per_Gasym, col = "purple2", lwd = 2)

