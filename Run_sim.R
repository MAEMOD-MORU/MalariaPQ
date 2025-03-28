# Run Simulation
library(deSolve)

#Read model
source("model.r")

# Initial Population
initP <- 70000000

# Initial States from UI Inputs
init_total_IAs <- 100
init_total_IAr <- 10

prob_s = 0.5 # Ratio of infection on sensitive group
prob_r = 0.5 # Ratio of infection on resistance group

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
)

beta_s0 = 0.01
beta_s1 = 0.04
beta_s2 = 0.04
beta_r0 = 0.04
beta_r1 = 0.04
beta_r2 = 0.04

parameters <- list(
  mui = (1 / 50) / 12,
  muo = (1 / 50) / 12,
  flo = 1,
  amp = 1,
  period = 12,
  phase = 0,

  beta_s = c(beta_s0, beta_s1, beta_s2),
  beta_r = c(beta_r0, beta_r1,beta_r2),
  beta_as = beta_s0,
  beta_ar = beta_r0,

  alpha = 1 / 12,
  tau1 = 1 / 4 / 12,
  tau2 = 1 / 26 / 12,
  tau = tau,

  d0 = (1 / 44) / 30, # Gametocyte recovery rate from Palang paper, table 2 at row Placebo (DP), column microscopy
  d1 = (1 / 20) / 30, # Can't find the ACT
  d2 = (1 / 7.88) / 30, # Gametocyte recovery rate from Palang paper, table 2 at row 0.0625 mg/kg, column microscopy

  gamma_i = 0.02,
  gamma_a = 0.02,
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

times <- seq(0, 160, by = 1)
out<- ode(y = state, times = times, func = Malaria_model_with_Array, parms = parameters)

plot(out,select=c("inc","inc_s"),type="l")

#matrix
n_run <- 100

matrix_output <- matrix(nrow = length(times),ncol = n_run)

for (i in 1:n_run) {
  # d0 <- (10.9–75.6)
  # d1 <- (2.21–60.0)
  # d2 <- (1.37–11.0)
  rd0 <- runif(n=1,min = 10.9,max = 75.9)
  rd1 <- runif(n=1,min = 2.21,max = 60.0)
  rd2 <- runif(n=1,min = 1.37,max = 11.0)
  parameters$d0 <- 1/rd0/30
  parameters$d1 <- 1/rd1/30
  parameters$d2 <- 1/rd2/30

  out<- ode(y = state, times = times, func = Malaria_model_with_Array, parms = parameters)
  matrix_output[,i] <- out[,"inc"]
}
plot(matrix_output[,1],type = "l")
for (y in 2:n_run) {
  lines(matrix_output[,y],col=y)
}

