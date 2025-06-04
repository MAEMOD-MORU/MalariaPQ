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

prob_s = 0.25 # Ratio of infection on sym/asym sensitive group
prob_r = 0.25 # Ratio of infection on sym/asym resistance group

Is_total <- init_total_IAs * prob_s
As_total <- init_total_IAs * (1-prob_s)

Ir_total <- init_total_IAr * prob_r
Ar_total <- init_total_IAr * (1-prob_r)

init_D1s <- 1/3
init_D2s <- 1/3

init_D1r <- 1/3
init_D2r <- 1/3

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

initStis <- c(0,0) # Succeed to Treatment (Sensitive) (D1,D2)
initFis <- c(0,0) # Fail to Treatment (Sensitive) (D1,D2)

initStir <- c(0,0) # Succeed to Treatment (Resistant) (D1,D2)
initFir <- c(0,0) # Fail to Treatment (Resistant) (D1,D2)

# Calculate initial susceptible population
initS <- initP - sum(initIs) - sum(initGIs) - initGAs - initAs - initRs -
  sum(initIr) - sum(initGIr) - initGAr - initAr - initRr

state <- c(
  S = initS,
  Is0 = initIs[1], Is1 = initIs[2], Is2 = initIs[3], # Flatten Is array into individual components
  Stis1 = initStis[1],Stis2=initStis[2],
  Fis1 = initFis[1],Fis2=initFis[2],
  GIs0 = initGIs[1], GIs1 = initGIs[2], GIs2 = initGIs[3],  # Gametocyte from Is
  GAs = initGAs,  # Gametocyte from As
  As = initAs,
  Rs = initRs,
  Ir0 = initIr[1], Ir1 = initIr[2], Ir2 = initIr[3], # Flatten Ir array into individual components
  Stir1 = initStir[1],Stir2 = initStir[2],
  Fir1 = initFir[1],Fir2 = initFir[2],
  GIr0 = initGIr[1], GIr1 = initGIr[2], GIr2 = initGIr[3],  # Gametocyte from Ir
  GAr = initGAr,  # Gametocyte from Ar
  Ar = initAr,
  Rr = initRr
)


# Dynamic tau from UI
τnfs <- rep(4.3,2) # Succeed to Treatment group to Recovered (Sensitive)
τigs <- c(0.33,0.75,3.7) # Symptomatic Gametocytes to Recovered (Sensitive)
τfs <-rep(0.33,2) # Fail to Treatment group to Recovered (Sensitive)
τntsd0 <- 0.33 # D0 to Recovered (Sensitive)
Tas <- 0.33 # Asymptomatic to Recovered (Sensitive)
Tags <- 0.33 # Asymptomatic Gametocytes to Recovered (Sensitive)

τnfr <- rep(4.3,2)  # Succeed to Treatment group to Recovered (Resistant)
τigr <- c(0.33,0.75,3.7) # Symptomatic Gametocytes to Recovered (Resistant)
τfr <- rep(0.33,2) # Fail to Treatment group to Recovered (Resistant)
τntrd0 <- 0.33 # D0 to Recovered (Resistant)
Tar <- 0.33 # Asymptomatic to Recovered (Resistant)
Tagr <- 0.33 # Asymptomatic Gametocytes to Recovered (Resistant)

beta_s0 = 1.6
beta_s1 = 1.6
beta_s2 = 1.6
beta_r0 = 1.6
beta_r1 = 1.6
beta_r2 = 1.6

parameters_Fail_rate_1 <- list(
  mui = (1 / 50) / 12,
  muo = (1 / 50) / 12,
  flo = 1,
  amp = 0.25,
  period = 12,
  phase = 0,
  prob_s = prob_s,
  prob_r = prob_r,
  start_d = 3620,# 10 years after 3500 months
  t_long = 12*10, # 2,5,10 years 
  Fail_rate_s = 0.1, # Fail Treatment of ACT
  Fail_rate_r = 0.1,
  
  beta_s = c(beta_s0, beta_s1, beta_s2),
  beta_r = c(beta_r0, beta_r1,beta_r2),
  beta_as = beta_s0,
  beta_ar = beta_r0,
  
  alpha = 1.07 / 12, #  doi:10.1371/journal.pone.0001767
  τigs=τigs,
  tau_as=Tas,
  tau_ags=Tagr,
  
  τigr=τigr,
  tau_ar=Tar,
  tau_agr=Tags,
  
  tau_ntsd0= τntsd0,
  tau_ntrd0= τntrd0,
  
  τnfs = τnfs ,
  τfs = τfs ,
  
  τnfr=τnfr,
  τfr=τfr,
  
  d0 = (1 / 90) * 30, # no drugs
  d1 = (1 / 44.1) * 30, # Gametocyte recovery rate from Palang paper, table 2 at row Placebo (DP), column microscopy
  d2 = (1 / 7.88) * 30, # Gametocyte recovery rate from Palang paper, table 2 at row 0.0625 mg/kg, column microscopy
  
  gamma_is = c(0.1,0.1,0.1), #0.1
  gamma_ir = c(0.1,0.1,0.1),
  gamma_ar = 0.1,
  gamma_as = 0.1,
  gamma_infs = c(0.01,0.01),
  gamma_ifs = c(1,1) *10, #1
  gamma_infr = c(0.01,0.01),
  gamma_ifr = c(1,1) *10,
  prob_lam = 0.7
)

parameters_Fail_rate_2 <- list(
  mui = (1 / 50) / 12,
  muo = (1 / 50) / 12,
  flo = 1,
  amp = 0.25,
  period = 12,
  phase = 0,
  prob_s = prob_s,
  prob_r = prob_r,
  start_d = 3620,# 10 years after 3500 months
  t_long = 12*10, # 2,5,10 years 
  Fail_rate_s = 0.1, # Fail Treatment of ACT
  Fail_rate_r = 0.3,
  
  beta_s = c(beta_s0, beta_s1, beta_s2),
  beta_r = c(beta_r0, beta_r1,beta_r2),
  beta_as = beta_s0,
  beta_ar = beta_r0,
  
  alpha = 1.07 / 12, #  doi:10.1371/journal.pone.0001767
  τigs=τigs,
  tau_as=Tas,
  tau_ags=Tagr,
  
  τigr=τigr,
  tau_ar=Tar,
  tau_agr=Tags,
  
  tau_ntsd0= τntsd0,
  tau_ntrd0= τntrd0,
  
  τnfs = τnfs ,
  τfs = τfs ,
  
  τnfr=τnfr,
  τfr=τfr,
  
  d0 = (1 / 90) * 30, # no drugs
  d1 = (1 / 44.1) * 30, # Gametocyte recovery rate from Palang paper, table 2 at row Placebo (DP), column microscopy
  d2 = (1 / 7.88) * 30, # Gametocyte recovery rate from Palang paper, table 2 at row 0.0625 mg/kg, column microscopy
  
  gamma_is = c(0.1,0.1,0.1), #0.1
  gamma_ir = c(0.1,0.1,0.1),
  gamma_ar = 0.1,
  gamma_as = 0.1,
  gamma_infs = c(0.01,0.01),
  gamma_ifs = c(1,1) *10, #1
  gamma_infr = c(0.01,0.01),
  gamma_ifr = c(1,1) *10,
  prob_lam = 0.7
)

parameters_Fail_rate_3 <- list(
  mui = (1 / 50) / 12,
  muo = (1 / 50) / 12,
  flo = 1,
  amp = 0.25,
  period = 12,
  phase = 0,
  prob_s = prob_s,
  prob_r = prob_r,
  start_d = 3620,# 10 years after 3500 months
  t_long = 12*10, # 2,5,10 years 
  Fail_rate_s = 0.1, # Fail Treatment of ACT
  Fail_rate_r = 0.5,
  
  beta_s = c(beta_s0, beta_s1, beta_s2),
  beta_r = c(beta_r0, beta_r1,beta_r2),
  beta_as = beta_s0,
  beta_ar = beta_r0,
  
  alpha = 1.07 / 12, #  doi:10.1371/journal.pone.0001767
  τigs=τigs,
  tau_as=Tas,
  tau_ags=Tagr,
  
  τigr=τigr,
  tau_ar=Tar,
  tau_agr=Tags,
  
  tau_ntsd0= τntsd0,
  tau_ntrd0= τntrd0,
  
  τnfs = τnfs ,
  τfs = τfs ,
  
  τnfr=τnfr,
  τfr=τfr,
  
  d0 = (1 / 90) * 30, # no drugs
  d1 = (1 / 44.1) * 30, # Gametocyte recovery rate from Palang paper, table 2 at row Placebo (DP), column microscopy
  d2 = (1 / 7.88) * 30, # Gametocyte recovery rate from Palang paper, table 2 at row 0.0625 mg/kg, column microscopy
  
  gamma_is = c(0.1,0.1,0.1), #0.1
  gamma_ir = c(0.1,0.1,0.1),
  gamma_ar = 0.1,
  gamma_as = 0.1,
  gamma_infs = c(0.01,0.01),
  gamma_ifs = c(1,1) *10, #1
  gamma_infr = c(0.01,0.01),
  gamma_ifr = c(1,1) *10,
  prob_lam = 0.7
)

parameters_Fail_rate_4 <- list(
  mui = (1 / 50) / 12,
  muo = (1 / 50) / 12,
  flo = 1,
  amp = 0.25,
  period = 12,
  phase = 0,
  prob_s = prob_s,
  prob_r = prob_r,
  start_d = 3620,# 10 years after 3500 months
  t_long = 12*10, # 2,5,10 years 
  Fail_rate_s = 0.1, # Fail Treatment of ACT
  Fail_rate_r = 0.7,
  
  beta_s = c(beta_s0, beta_s1, beta_s2),
  beta_r = c(beta_r0, beta_r1,beta_r2),
  beta_as = beta_s0,
  beta_ar = beta_r0,
  
  alpha = 1.07 / 12, #  doi:10.1371/journal.pone.0001767
  τigs=τigs,
  tau_as=Tas,
  tau_ags=Tagr,
  
  τigr=τigr,
  tau_ar=Tar,
  tau_agr=Tags,
  
  tau_ntsd0= τntsd0,
  tau_ntrd0= τntrd0,
  
  τnfs = τnfs ,
  τfs = τfs ,
  
  τnfr=τnfr,
  τfr=τfr,
  
  d0 = (1 / 90) * 30, # no drugs
  d1 = (1 / 44.1) * 30, # Gametocyte recovery rate from Palang paper, table 2 at row Placebo (DP), column microscopy
  d2 = (1 / 7.88) * 30, # Gametocyte recovery rate from Palang paper, table 2 at row 0.0625 mg/kg, column microscopy
  
  gamma_is = c(0.1,0.1,0.1), #0.1
  gamma_ir = c(0.1,0.1,0.1),
  gamma_ar = 0.1,
  gamma_as = 0.1,
  gamma_infs = c(0.01,0.01),
  gamma_ifs = c(1,1) *10, #1
  gamma_infr = c(0.01,0.01),
  gamma_ifr = c(1,1) *10,
  prob_lam = 0.7
)

parameters_Fail_rate_5 <- list(
  mui = (1 / 50) / 12,
  muo = (1 / 50) / 12,
  flo = 1,
  amp = 0.25,
  period = 12,
  phase = 0,
  prob_s = prob_s,
  prob_r = prob_r,
  start_d = 3620,# 10 years after 3500 months
  t_long = 12*10, # 2,5,10 years 
  Fail_rate_s = 0.1, # Fail Treatment of ACT
  Fail_rate_r = 0.9,
  
  beta_s = c(beta_s0, beta_s1, beta_s2),
  beta_r = c(beta_r0, beta_r1,beta_r2),
  beta_as = beta_s0,
  beta_ar = beta_r0,
  
  alpha = 1.07 / 12, #  doi:10.1371/journal.pone.0001767
  τigs=τigs,
  tau_as=Tas,
  tau_ags=Tagr,
  
  τigr=τigr,
  tau_ar=Tar,
  tau_agr=Tags,
  
  tau_ntsd0= τntsd0,
  tau_ntrd0= τntrd0,
  
  τnfs = τnfs ,
  τfs = τfs ,
  
  τnfr=τnfr,
  τfr=τfr,
  
  d0 = (1 / 90) * 30, # no drugs
  d1 = (1 / 44.1) * 30, # Gametocyte recovery rate from Palang paper, table 2 at row Placebo (DP), column microscopy
  d2 = (1 / 7.88) * 30, # Gametocyte recovery rate from Palang paper, table 2 at row 0.0625 mg/kg, column microscopy
  
  gamma_is = c(0.1,0.1,0.1), #0.1
  gamma_ir = c(0.1,0.1,0.1),
  gamma_ar = 0.1,
  gamma_as = 0.1,
  gamma_infs = c(0.01,0.01),
  gamma_ifs = c(1,1) *10, #1
  gamma_infr = c(0.01,0.01),
  gamma_ifr = c(1,1) *10,
  prob_lam = 0.7
)

times <- seq(0, 5000, by = 1)

out_30_Fail_rate <- ode(y = state, times = times, func = Malaria_model_with_Array_2, parms = parameters_Fail_rate_1)
out_40_Fail_rate <- ode(y = state, times = times, func = Malaria_model_with_Array_2, parms = parameters_Fail_rate_2)
out_50_Fail_rate <- ode(y = state, times = times, func = Malaria_model_with_Array_2, parms = parameters_Fail_rate_3)
out_60_Fail_rate <- ode(y = state, times = times, func = Malaria_model_with_Array_2, parms = parameters_Fail_rate_4)
out_70_Fail_rate <- ode(y = state, times = times, func = Malaria_model_with_Array_2, parms = parameters_Fail_rate_5)

# 70% Fail rate
plot(rowSums(out_70_Fail_rate[(3500-240):(3500),c("Fis1","Fis2","Fir1","Fir2")]),type="l",xlab="Months",ylab="Incidence",main="Fail to Treatment (2,5,10 year treatment)",col=1,ylim = c(0,1e3))
# 60% Fail rate
lines(rowSums(out_60_Fail_rate[(3500-240):(3500),c("Fis1","Fis2")]),type="l",col=2)
# 50% Fail rate
lines(rowSums(out_50_Fail_rate[(3500-240):(3500),c("Fis1","Fis2")]),type="l",col=3)
# 40% Fail rate
lines(rowSums(out_40_Fail_rate[(3500-240):(3500),c("Fis1","Fis2")]),type="l",col=4)
# 30% Fail rate
lines(rowSums(out_30_Fail_rate[(3500-240):(3500),c("Fis1","Fis2")]),type="l",col=5)

legend("topright", legend = c("90% Fail rate", "70% Fail rate", "50% Fail rate","30% Fail rate","10% Fail rate"),
       col = 1:5,
       lwd = 2, lty = 1, bg = "white",cex = 1)


# 70% Fail rate
plot(rowSums(out_70_Fail_rate[(3500-240):(3500),c("GIr0","GIr1","GIr2")]),type="l",xlab="Months",ylab="Number of People Carrying Gametocytes",main="People Carrying Gametocytes (Resistance)",col=1,ylim = c(0,1.5e5))
# 60% Fail rate
lines(rowSums(out_60_Fail_rate[(3500-240):(3500),c("GIr0","GIr1","GIr2")]),type="l",col=2)
# 50% Fail rate
lines(rowSums(out_50_Fail_rate[(3500-240):(3500),c("GIr0","GIr1","GIr2")]),type="l",col=3)
# 40% Fail rate
lines(rowSums(out_40_Fail_rate[(3500-240):(3500),c("GIr0","GIr1","GIr2")]),type="l",col=4)
# 30% Fail rate
lines(rowSums(out_30_Fail_rate[(3500-240):(3500),c("GIr0","GIr1","GIr2")]),type="l",col=5)

legend("topleft", legend = c("90% Fail rate", "70% Fail rate", "50% Fail rate","30% Fail rate","10% Fail rate"),
       col = 1:5,
       lwd = 2, lty = 1, bg = "white",cex = 1)

# 70% Fail rate
plot(rowSums(out_70_Fail_rate[(3500):(3500+360),c("GIr0","GIr1","GIr2")]),type="l",xlab="Months",ylab="Number of People Carrying Gametocytes",main="People Carrying Gametocytes (Resistance)",col=1,ylim = c(0,1.5e5))
# 60% Fail rate
lines(rowSums(out_60_Fail_rate[(3500):(3500+360),c("GIr0","GIr1","GIr2")]),type="l",col=2)
# 50% Fail rate
lines(rowSums(out_50_Fail_rate[(3500):(3500+360),c("GIr0","GIr1","GIr2")]),type="l",col=3)
# 40% Fail rate
lines(rowSums(out_40_Fail_rate[(3500):(3500+360),c("GIr0","GIr1","GIr2")]),type="l",col=4)
# 30% Fail rate
lines(rowSums(out_30_Fail_rate[(3500):(3500+360),c("GIr0","GIr1","GIr2")]),type="l",col=5)

abline(v=120,col="blue")
text(85,110000,"start\n treatment")
abline(v=240,col="black")
text(210,50000,"stop\n treatment")

legend("topright", legend = c("70% Fail rate", "60% Fail rate", "50% Fail rate","40% Fail rate","30% Fail rate"),
       col = 1:5,
       lwd = 2, lty = 1, bg = "white",cex = 1)

# 70% Fail rate
plot(out_70_Fail_rate[(3500-240):(3500),"GIr2"],type="l",xlab="Months",ylab="Number of People Carrying Gametocytes",main="People Carrying Gametocytes (Resistance)",col=1,ylim = c(0,1e5))
lines(out_70_Fail_rate[(3500-240):(3500),"GIr0"],type="l",col=2)

lines(out_30_Fail_rate[(3500-240):(3500),"GIr0"],type="l",col=3,lty=2)
lines(out_30_Fail_rate[(3500-240):(3500),"GIr2"],type="l",col=3)
lines(out_40_Fail_rate[(3500-240):(3500),"GIr0"],type="l",col=2,lty=2)
lines(out_40_Fail_rate[(3500-240):(3500),"GIr2"],type="l",col=2)

