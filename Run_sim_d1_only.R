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

init_D1s <- 2/3
init_D2s <- 0/3

init_D1r <- 2/3
init_D2r <- 0/3

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

beta_s0 = 1.45
beta_s1 = 1.45
beta_s2 = 1.45
beta_r0 = 1.45
beta_r1 = 1.45
beta_r2 = 1.45

parameters <- list(
  mui = (1 / 50) / 12,
  muo = (1 / 50) / 12,
  flo = 1,
  amp = 0.25,
  period = 12,
  phase = 0,
  prob_s = prob_s,
  prob_r = prob_r,
  start_d = 4620,# 10 years after 3500 months
  t_long = 12*10, # 2,5,10 years 
  Fail_rate_s = 0.1, # Fail Treatment of ACT
  Fail_rate_r = 0.5,
  
  beta_s = c(beta_s0, beta_s1, beta_s2),
  beta_r = c(beta_r0, beta_r1,beta_r2),
  beta_as = beta_s0,
  beta_ar = beta_r0,
  
  alpha = 1.07 / 12, #  doi:10.1371/journal.pone.0001767
  τigs=τigs,
  Tas=Tas,
  Tags=Tags,
  
  τigr=τigr,
  Tar=Tar,
  Tagr=Tagr,
  
  τntsd0= τntsd0,
  τntrd0= τntrd0,
  
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
  gamma_ifs = c(1,1) , #1
  gamma_infr = c(0.01,0.01),
  gamma_ifr = c(1,1) ,
  prob_lam = 0.7
)

parameters_5 <- list(
  mui = (1 / 50) / 12,
  muo = (1 / 50) / 12,
  flo = 1,
  amp = 0.25,
  period = 12,
  phase = 0,
  prob_s = prob_s,
  prob_r = prob_r,
  start_d = 4620,# 10 years after 3500 months
  t_long = 12*5, # 2,5,10 years 
  Fail_rate_s = 0.1, # Fail Treatment of ACT
  Fail_rate_r = 0.5,
  
  beta_s = c(beta_s0, beta_s1, beta_s2),
  beta_r = c(beta_r0, beta_r1,beta_r2),
  beta_as = beta_s0,
  beta_ar = beta_r0,
  
  alpha = 1.07 / 12, #  doi:10.1371/journal.pone.0001767
  τigs=τigs,
  Tas=Tas,
  Tags=Tags,
  
  τigr=τigr,
  Tar=Tar,
  Tagr=Tagr,
  
  τntsd0= τntsd0,
  τntrd0= τntrd0,
  
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
  gamma_ifs = c(1,1) , #1
  gamma_infr = c(0.01,0.01),
  gamma_ifr = c(1,1) ,
  prob_lam = 0.7
)

parameters_2 <- list(
  mui = (1 / 50) / 12,
  muo = (1 / 50) / 12,
  flo = 1,
  amp = 0.25,
  period = 12,
  phase = 0,
  prob_s = prob_s,
  prob_r = prob_r,
  start_d = 4620,# 10 years after 3500 months
  t_long = 12*2, # 2,5,10 years 
  Fail_rate_s = 0.1, # Fail Treatment of ACT
  Fail_rate_r = 0.5,
  
  beta_s = c(beta_s0, beta_s1, beta_s2),
  beta_r = c(beta_r0, beta_r1,beta_r2),
  beta_as = beta_s0,
  beta_ar = beta_r0,
  
  alpha = 1.07 / 12, #  doi:10.1371/journal.pone.0001767
  τigs=τigs,
  Tas=Tas,
  Tags=Tags,
  
  τigr=τigr,
  Tar=Tar,
  Tagr=Tagr,
  
  τntsd0= τntsd0,
  τntrd0= τntrd0,
  
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
  gamma_ifs = c(1,1) , #1
  gamma_infr = c(0.01,0.01),
  gamma_ifr = c(1,1) ,
  prob_lam = 0.7
)

times <- seq(0, 5000, by = 1)
out<- ode(y = state, times = times, func = Malaria_model_with_Array_2, parms = parameters)
out_5<- ode(y = state, times = times, func = Malaria_model_with_Array_2, parms = parameters_5)
out_2<- ode(y = state, times = times, func = Malaria_model_with_Array_2, parms = parameters_2)



out

plot(out[,"inc_s"],type="l",xlab="Months",ylab="Incidence",main="Incidence Sensitive/Resistant",col=1)
lines(out[,"inc_r"],type="l",col=2)
legend("topleft",c("Sensitive","Resistant"),col = 1:2,lty =c(1,1))

plot(out[,"inc"],type="l",xlab="Months",ylab="Incidence",main="Total Incidence")

plot(out[,"malaria_case"],type="l",xlab="Months",ylab="Incidence",main="Total Malaria Case")
lines(out[,"total_G"],col=2)
lines(out[,"total_inf"],col=3)

plot(out[,"Ir"],type="l",xlab="Months",ylab="Incidence",main="Total Malaria Case")

plot(out[960:1200,"malaria_case"],type="l",xlab="Months",ylab="Incidence",main="Total Malaria Case")


plot(out[0:4500,"inc"],type="l",xlab="Months",ylab="Incidence",main="Total Incidence (No treatment)")
plot(out[4000:4012,"inc"],type="l",xlab="Months",ylab="Incidence",main="Total Incidence")
plot(0:240,out[(4500-240):4500,"inc"],type="l",xlab="Months",ylab="Incidence",main="Total Incidence (No treatment)")
plot(out[,"N"],type="l",xlab="Months",ylab="Total")
plot(out[,"N"],type="l",xlab="Months",ylab="Population",main="Total Population", ylim = c(70000050,70000050))

plot(cumsum(out[0:3500,"inc"]))

plot(out[0:3500,"Ir0"],type = "l")

plot(out[0:3500,"Fis1"],type = "l")
lines(out[0:3500,"Stis1"],type = "l",col=2)

plot(out[(3500-240):3500,"inc_s"],type="l",xlab="Months",ylab="Incidence",main="Incidence Sensitive/Resistant",col=1,ylim = c(0,1.2e5))
lines(out[(3500-240):3500,"inc_r"],type="l",col=2)
legend("bottomright",c("Sensitive","Resistant"),col = 1:2,lty =c(1,1),cex = 0.9)

# 10 year
plot(out[(4500):(4500+360),"inc"],type="l",xlab="Months",ylab="Incidence",main="Total Incidence (10 year treatment)",col=1,ylim = c(0,1.5e5))
# lines(out[(3500):(3500+360),"inc_r"],type="l",col=2)
# legend("topright",c("Sensitive","Resistant"),col = 1:2,lty =c(1,1),cex = 0.9)
abline(v=120,col="green")
text(170,140000,"start treatment")
abline(v=240,col="blue")
text(300,140000,"Stop treatment")

# 5 year
plot(out_5[(4500):(4500+300),"inc"],type="l",xlab="Months",ylab="Incidence",main="Total Incidence (5 year treatment)",col=1,ylim = c(0,1.5e5))
# lines(out[(3500):(3500+360),"inc_r"],type="l",col=2)
# legend("topright",c("Sensitive","Resistant"),col = 1:2,lty =c(1,1),cex = 0.9)
abline(v=120,col="green")
text(150,140000,"start\n treatment")
abline(v=180,col="blue")
text(230,140000,"Stop treatment")

# 2 year
plot(out[(4500):(4500+300-36),"inc"],type="l",xlab="Months",ylab="Incidence",main="Total Incidence (2 year treatment)",col=1,ylim = c(0,1.5e5))
# lines(out[(3500):(3500+360),"inc_r"],type="l",col=2)
# legend("topright",c("Sensitive","Resistant"),col = 1:2,lty =c(1,1),cex = 0.9)
abline(v=120,col="green")
text(90,50000,"start\n treatment")
abline(v=180-36,col="blue")
text(180,140000,"Stop treatment")

plot(out[960:1200,"inc_s"],type="l",xlab="Months",ylab="Incidence",main="Incidence Sensitive",col=1,ylim = c(4e5,7e5))
plot(out[960:1200,"inc_r"],type="l",xlab="Months",ylab="Incidence",main="Incidence Resistant",col=2,ylim = c(4e5,7e5))

plot(out[(3500-240):3500,"inc_asym"],type="l",xlab="Months",ylab="Incidence",main="Incidence Asymptomatic/Symptomatic",col=3,ylim = c(0,1.2e5))
lines(out[(3500-240):3500,"inc_sym"],type="l",col=4)
legend("right",c("Asymptomatic","Symptomatic"),col = 3:4,lty =c(1,1),cex = 0.9)

plot(out[960:1200,"Gs"],type="l",xlab="Months",ylab="Component",main="Gametocyte Infections Sensitive/Resistant",col=1,ylim = c(3e5,8e5))
lines(out[960:1200,"Gr"],type="l",col=2)
legend("topleft",c("Sensitive","Resistant"),col = 1:2,lty =c(1,1))

plot(out[960:1200,"Gs"],type="l",xlab="Months",ylab="Incidence",main="Gametocyte Infections Sensitive",col=1)
plot(out[960:1200,"Gr"],type="l",xlab="Months",ylab="Incidence",main="Gametocyte Infections Resistant",col=2)

plot(out[(4500-240):4500,"GS_inc"],type="l",xlab="Months",ylab="Incidence",main="Gametocyte Infections Sensitive/Resistant",col=1,ylim = c(0e5,2.5e4))
lines(out[(4500-240):4500,"GR_inc"],type="l",col=2)
legend("bottomright",c("Sensitive","Resistant"),col = 1:2,lty =c(1,1))

plot(out[(4500-240):4500,"Gasym_inc"],type="l",xlab="Months",ylab="Incidence",main="Gametocyte Infections Sensitive/Resistant",col=3,ylim = c(0e5,3e4))
lines(out[(4500-240):4500,"Gsym_inc"],type="l",col=4)
legend("bottomright",c("Asymptomatic","Symptomatic"),col = 3:4,lty =c(1,1))

plot(rowSums(out[(4500-240):4500,c("Fir1","Fir2")]),type="l",xlab="Months",ylab="Patient",main="Fail to Treatment Patient",col=2,ylim = c(0,2.5e3))
lines(rowSums(out[(4500-240):4500,c("Fis1","Fis2")]),type="l",col=1)
legend("bottomright",c("Sensitive","Resistant"),col = 1:2,lty =c(1,1))

plot(rowSums(out[960:1200,c("Fir1","Fir2")]),type="l",xlab="Months",ylab="Incidence",main="Fail to Treatment Sensitive",col=1)
plot(rowSums(out[960:1200,c("Fis1","Fis2")]),type="l",xlab="Months",ylab="Incidence",main="Fail to Treatment Resistant",col=2)

plot(out[0:4500,c("GIs0")],type="l",xlab="Months",ylab="Patient",main="Fail to Treatment Patient",col="green",ylim = c(0,6e3))
lines(out[0:4500,c("GIs1")],col="red")
lines(out[0:4500,c("GIs2")],col="blue")
lines(out[0:4500,c("GAs")],col=1)

plot(out[0:4500,c("GIr0")],type="l",xlab="Months",ylab="Patient",main="Fail to Treatment Patient",col="green",ylim = c(0,6e3))
lines(out[0:4500,c("GIr1")],col="red")
lines(out[0:4500,c("GIr2")],col="blue",lty=3,lwd=0.5)
lines(out[0:4500,c("GAr")],col=1)

plot(out)

# incidence
plot(x=0:240,out[4500:4740,"inc"],type="l",xlab="months",ylab="incidence")

plot(x=0:240,out[(1000-240):1000,"inc"],type="l",xlab="months",ylab="incidence")



# gametocyte sensitive
plot(x=0:240, out[4500:4740,c("GIs0")],type="l",xlab="months",ylab="gametocyte sensitive")
lines(out[4500:4740,c("GIs1")],col=2)
lines(out[4500:4740,c("GIs2")],col=3)
title("Sensitive")

# gametocyte res
plot(x=0:240, out[4500:4740,c("GIr0")],type="l",ylim=c(0e5,1e4),xlab="months",ylab="Gametocytes Resistance")
lines(out[4500:4740,c("GIr1")],col=2)
lines(out[4500:4740,c("GIr2")],col=3)

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

matrix_output_inc_s <- matrix(nrow = length(times),ncol = n_run)
matrix_output_inc_r <- matrix(nrow = length(times),ncol = n_run)

matrix_output_inc_sym <- matrix(nrow = length(times),ncol = n_run)
matrix_output_inc_asym <- matrix(nrow = length(times),ncol = n_run)

matrix_output_Gs <- matrix(nrow = length(times),ncol = n_run)
matrix_output_Gr <- matrix(nrow = length(times),ncol = n_run)

matrix_output_Gs_sym <- matrix(nrow = length(times),ncol = n_run)
matrix_output_Gs_asym <- matrix(nrow = length(times),ncol = n_run)

matrix_output_Gr_sym <- matrix(nrow = length(times),ncol = n_run)
matrix_output_Gr_asym <- matrix(nrow = length(times),ncol = n_run)

matrix_output_Gsym <- matrix(nrow = length(times),ncol = n_run)
matrix_output_Gasym <- matrix(nrow = length(times),ncol = n_run)

matrix_output_Gsym_inc_s <- matrix(nrow = length(times),ncol = n_run)
matrix_output_Gsym_inc_r <- matrix(nrow = length(times),ncol = n_run)

matrix_output_Gasym_inc_s <- matrix(nrow = length(times),ncol = n_run)
matrix_output_Gasym_inc_r <- matrix(nrow = length(times),ncol = n_run)

pb <- txtProgressBar(max =n_run)
set.seed(283)
for (i in 1:n_run) {
  # 0.33 (90 days)
  # d1 ~ Uniform(10.9–75.6)
  # d2 ~ Uniform(2.21–60.0)
  # rd0 <- runif(n=1,min = 10.9,max = 44) # 44.1
  rd1 <- runif(n=1,min = 10.9,max = 75.6) #7.88
  rd2 <- runif(n=1,min = 2.21,max = 60.0) #4.88
  # parameters$d0 <- 1/rd0*30
  parameters$d1 <- 1/rd1*30
  parameters$d2 <- 1/rd2*30
  
  out<- ode(y = state, times = times, func = Malaria_model_with_Array_2, parms = parameters)
  matrix_output_inc[,i] <- out[,"inc"]
  
  matrix_output_GIs0[,i] <- out[,"GIs0"]
  matrix_output_GIs1[,i] <- out[,"GIs1"]
  matrix_output_GIs2[,i] <- out[,"GIs2"]
  
  matrix_output_GIr0[,i] <- out[,"GIr0"]
  matrix_output_GIr1[,i] <- out[,"GIr1"]
  matrix_output_GIr2[,i] <- out[,"GIr2"]
  
  matrix_output_total_inf[,i] <- out[,"total_inf"]
  matrix_output_total_G[,i] <- out[,"total_G"]
  
  matrix_output_inc_s[,i] <- out[,"inc_s"]
  matrix_output_inc_r[,i] <- out[,"inc_r"]
  
  matrix_output_inc_sym[,i] <- out[,"inc_sym"]
  matrix_output_inc_asym[,i] <- out[,"inc_asym"]
  
  matrix_output_Gs[,i] <- out[,"Gs"]
  matrix_output_Gr[,i] <- out[,"Gr"]
  
  matrix_output_Gs_sym[,i] <- out[,"Gs_sym"]
  matrix_output_Gs_asym[,i] <- out[,"Gs_asym"]
  
  matrix_output_Gr_sym[,i] <- out[,"Gr_sym"]
  matrix_output_Gr_asym[,i] <- out[,"Gr_asym"]
  
  matrix_output_Gsym[,i] <- out[,"Gsym"]
  matrix_output_Gasym[,i] <- out[,"Gasym"]
  
  matrix_output_Gsym_inc_s[,i] <- out[,"Gsym_inc_s"]
  matrix_output_Gsym_inc_r[,i] <- out[,"Gsym_inc_r"]
  
  matrix_output_Gasym_inc_s[,i] <- out[,"Gasym_inc_s"]
  matrix_output_Gasym_inc_r[,i] <- out[,"Gasym_inc_r"]
  
  setTxtProgressBar(pb, i)
}

# matrix_output <- run_matrix_output(100,360)


plot_out_line(matrix_output_inc[4500:4740,])
plot_out_line(matrix_output_GIs0[4500:4740,])
plot_out_line(matrix_output_GIs1[4500:4740,])
plot_out_line(matrix_output_GIs2[4500:4740,])

plot_out_CI95(matrix_output_inc[4500:4740,])
plot_out_CI95(matrix_output_GIs0[4500:4740,])
plot_out_CI95(matrix_output_GIs1[4500:4740,])
plot_out_CI95(matrix_output_GIs2[4500:4740,])


#### inc ######
# Compute median, lower (2.5th), and upper (97.5th) percentiles
median_inc <- apply(matrix_output_inc[4500:4740,], 1, median, na.rm = TRUE)
lower_inc <- apply(matrix_output_inc[4500:4740,], 1, quantile, probs = 0.025, na.rm = TRUE)
upper_inc <- apply(matrix_output_inc[4500:4740,], 1, quantile, probs = 0.975, na.rm = TRUE)

# Plot with confidence interval
plot(1:nrow(matrix_output_inc[4500:4740,]), median_inc, type = "l", col = "blue", lwd = 1.5, ylim = range(lower_inc, upper_inc),
     ylab = "Incidence", xlab = "Months", main = "Total Incidence with 95% CI")

# Add shaded region for confidence interval
polygon(c(1:nrow(matrix_output_inc[4500:4740,]), rev(1:nrow(matrix_output_inc[4500:4740,]))), c(lower_inc, rev(upper_inc)), col = rgb(0.5, 0.5, 0.5, 0.75), border = NA)
lines(1:nrow(matrix_output_inc[4500:4740,]),lower_inc, col="grey50",lty=2)
lines(1:nrow(matrix_output_inc[4500:4740,]),upper_inc, col="grey50",lty=2)

# Add median line
lines(1:nrow(matrix_output_inc[4500:4740,]), median_inc, col = "blue", lwd = 1.5)
abline(v=120,col="red",lwd=2,lty = 2)
text(85,50000,"Start Treatment")
abline(v=180,col="purple3",lwd=2,lty = 2)
text(160,150000,"5 year")

#### gametocyte sensitive ######
median_GIs0 <- apply(matrix_output_GIs0[4500:4740,], 1, median, na.rm = TRUE)
lower_GIs0 <- apply(matrix_output_GIs0[4500:4740,], 1, quantile, probs = 0.025, na.rm = TRUE)
upper_GIs0 <- apply(matrix_output_GIs0[4500:4740,], 1, quantile, probs = 0.975, na.rm = TRUE)

median_GIs1 <- apply(matrix_output_GIs1[4500:4740,], 1, median, na.rm = TRUE)
lower_GIs1 <- apply(matrix_output_GIs1[4500:4740,], 1, quantile, probs = 0.025, na.rm = TRUE)
upper_GIs1 <- apply(matrix_output_GIs1[4500:4740,], 1, quantile, probs = 0.975, na.rm = TRUE)

median_GIs2 <- apply(matrix_output_GIs2[4500:4740,], 1, median, na.rm = TRUE)
lower_GIs2 <- apply(matrix_output_GIs2[4500:4740,], 1, quantile, probs = 0.025, na.rm = TRUE)
upper_GIs2 <- apply(matrix_output_GIs2[4500:4740,], 1, quantile, probs = 0.975, na.rm = TRUE)


# Plot with confidence interval
plot(1:nrow(matrix_output_GIs0[4500:4740,]), median_GIs0, type = "l", col = "green", lwd = 1.5,ylim = c(0,10000),
     ylab = "Number of People Carrying Gametocytes", xlab = "Months", main = "Gametocyte Sensitive with 95% CI")

# Add shaded region for confidence interval
polygon(c(1:nrow(matrix_output_GIs0[4500:4740,]), rev(1:nrow(matrix_output_GIs0[4500:4740,]))), c(lower_GIs0, rev(upper_GIs0)), col = rgb(0.2, 0.5, 0.2, 0.25), border = NA)
lines(1:nrow(matrix_output_GIs0[4500:4740,]),lower_GIs0, col="green",lty=2)
lines(1:nrow(matrix_output_GIs0[4500:4740,]),upper_GIs0, col="green",lty=2)

polygon(c(1:nrow(matrix_output_GIs1[4500:4740,]), rev(1:nrow(matrix_output_GIs1[4500:4740,]))), c(lower_GIs1, rev(upper_GIs1)), col = rgb(0.5, 0.2, 0.2, 0.25), border = NA)
lines(1:nrow(matrix_output_GIs1[4500:4740,]),lower_GIs1, col="red",lty=2)
lines(1:nrow(matrix_output_GIs1[4500:4740,]),upper_GIs1, col="red",lty=2)

# polygon(c(1:nrow(matrix_output_GIs2[4500:4740,]), rev(1:nrow(matrix_output_inc[4500:4740,]))), c(lower_GIs2, rev(upper_GIs2)), col = rgb(0.2, 0.2, 0.5, 0.25), border = NA)
# lines(1:nrow(matrix_output_GIs2[4500:4740,]),lower_GIs2, col="blue",lty=2)
# lines(1:nrow(matrix_output_GIs2[4500:4740,]),upper_GIs2, col="blue",lty=2)

# Add median line
lines(1:nrow(matrix_output_GIs0[4500:4740,]), median_GIs0, col = rgb(0.1, 0.7, 0.1, 1), lwd = 2)
lines(1:nrow(matrix_output_GIs1[4500:4740,]), median_GIs1, col = rgb(0.7, 0.1, 0.1, 1), lwd = 2)
# lines(1:nrow(matrix_output_GIs2[4500:4740,]), median_GIs2, col = rgb(0.1, 0.1, 0.7, 1), lwd = 2)

abline(v=120,col="red",lwd=2,lty = 2)
text(85,2000,"Start Treatment")
abline(v=180,col="purple3",lwd=2,lty = 2)
text(165,7000,"5 year")

legend("topright", legend = c("No treatment", "ACT"),
       col = c(rgb(0.1, 0.7, 0.1, 1), rgb(0.7, 0.1, 0.1, 1)),
       lwd = 2, lty = 1, bg = "white",cex = 1)



#### gametocyte sensitive (Asym/Sym) ######
median_Gs_sym <- apply(matrix_output_Gs_sym[4500:4740,], 1, median, na.rm = TRUE)
lower_Gs_sym <- apply(matrix_output_Gs_sym[4500:4740,], 1, quantile, probs = 0.025, na.rm = TRUE)
upper_Gs_sym <- apply(matrix_output_Gs_sym[4500:4740,], 1, quantile, probs = 0.975, na.rm = TRUE)

median_Gs_asym <- apply(matrix_output_Gs_asym[4500:4740,], 1, median, na.rm = TRUE)
lower_Gs_asym <- apply(matrix_output_Gs_asym[4500:4740,], 1, quantile, probs = 0.025, na.rm = TRUE)
upper_Gs_asym <- apply(matrix_output_Gs_asym[4500:4740,], 1, quantile, probs = 0.975, na.rm = TRUE)

# Plot with confidence interval
plot(1:nrow(matrix_output_Gs_asym[4500:4740,]), median_Gs_asym, type = "l", col = "green", lwd = 1.5,ylim = c(0,55000),
     ylab = "Number of People Carrying Gametocytes", xlab = "Months", main = "Gametocyte Sensitive(Asym/Sym) with 95% CI")

# Add shaded region for confidence interval
polygon(c(1:nrow(matrix_output_Gs_asym[4500:4740,]), rev(1:nrow(matrix_output_Gs_asym[4500:4740,]))), c(lower_Gs_asym, rev(upper_Gs_asym)), col = rgb(0.2, 0.5, 0.2, 0.25), border = NA)
lines(1:nrow(matrix_output_Gs_asym[4500:4740,]),lower_Gs_asym, col="green4",lty=2)
lines(1:nrow(matrix_output_Gs_asym[4500:4740,]),upper_Gs_asym, col="green4",lty=2)

polygon(c(1:nrow(matrix_output_Gs_sym[4500:4740,]), rev(1:nrow(matrix_output_Gs_sym[4500:4740,]))), c(lower_Gs_sym, rev(upper_Gs_sym)), col = rgb(0.2, 0.2, 0.5, 0.25), border = NA)
lines(1:nrow(matrix_output_Gs_sym[4500:4740,]),lower_Gs_sym, col="blue",lty=2)
lines(1:nrow(matrix_output_Gs_sym[4500:4740,]),upper_Gs_sym, col="blue",lty=2)

lines(1:nrow(matrix_output_Gs_asym[4500:4740,]), median_Gs_asym, col = "green4", lwd = 2)
lines(1:nrow(matrix_output_Gs_sym[4500:4740,]), median_Gs_sym, col = "blue", lwd = 2)

abline(v=120,col="red",lwd=2,lty = 2)
text(85,30000,"Start Treatment")
abline(v=180,col="purple3",lwd=2,lty = 2)
text(160,40000,"5 year")

legend("topright",c("Asymptomatic","Symptomatic"),col =c("green4","blue"),lty =c(1,1), lwd = 2,cex = 1)


#### gametocyte resistance ######
median_GIr0 <- apply(matrix_output_GIr0[4500:4740,], 1, median, na.rm = TRUE)
lower_GIr0 <- apply(matrix_output_GIr0[4500:4740,], 1, quantile, probs = 0.025, na.rm = TRUE)
upper_GIr0 <- apply(matrix_output_GIr0[4500:4740,], 1, quantile, probs = 0.975, na.rm = TRUE)

median_GIr1 <- apply(matrix_output_GIr1[4500:4740,], 1, median, na.rm = TRUE)
lower_GIr1 <- apply(matrix_output_GIr1[4500:4740,], 1, quantile, probs = 0.025, na.rm = TRUE)
upper_GIr1 <- apply(matrix_output_GIr1[4500:4740,], 1, quantile, probs = 0.975, na.rm = TRUE)

median_GIr2 <- apply(matrix_output_GIr2[4500:4740,], 1, median, na.rm = TRUE)
lower_GIr2 <- apply(matrix_output_GIr2[4500:4740,], 1, quantile, probs = 0.025, na.rm = TRUE)
upper_GIr2 <- apply(matrix_output_GIr2[4500:4740,], 1, quantile, probs = 0.975, na.rm = TRUE)

# Plot with confidence interval
plot(1:nrow(matrix_output_GIr0[4500:4740,]), median_GIr0, type = "l", col = "green", lwd = 1.5, ylim = range(matrix_output_GIr0, matrix_output_GIr1),
     ylab = "Number of People Carrying Gametocytes", xlab = "Months", main = "Gametocytes Resistance with 95% CI")


# Add shaded region for confidence interval
polygon(c(1:nrow(matrix_output_GIr0[4500:4740,]), rev(1:nrow(matrix_output_GIr0[4500:4740,]))), c(lower_GIr0, rev(upper_GIr0)), col = rgb(0.2, 0.5, 0.2, 0.25), border = NA)
lines(1:nrow(matrix_output_GIr0[4500:4740,]),lower_GIr0, col="green",lty=2)
lines(1:nrow(matrix_output_GIr0[4500:4740,]),upper_GIr0, col="green",lty=2)

polygon(c(1:nrow(matrix_output_GIr1[4500:4740,]), rev(1:nrow(matrix_output_GIr1[4500:4740,]))), c(lower_GIr1, rev(upper_GIr1)), col = rgb(0.5, 0.2, 0.2, 0.25), border = NA)
lines(1:nrow(matrix_output_GIr1[4500:4740,]),lower_GIr1, col="red",lty=2)
lines(1:nrow(matrix_output_GIr1[4500:4740,]),upper_GIr1, col="red",lty=2)

# polygon(c(1:nrow(matrix_output_GIr2[4500:4740,]), rev(1:nrow(matrix_output_inc[4500:4740,]))), c(lower_GIr2, rev(upper_GIr2)), col = rgb(0.2, 0.2, 0.5, 0.25), border = NA)
# lines(1:nrow(matrix_output_GIr2[4500:4740,]),lower_GIr2, col="blue",lty=2)
# lines(1:nrow(matrix_output_GIr2[4500:4740,]),upper_GIr2, col="blue",lty=2)

# Add median line
lines(1:nrow(matrix_output_GIr0[4500:4740,]), median_GIr0, col = "green2", lwd = 2)
lines(1:nrow(matrix_output_GIr1[4500:4740,]), median_GIr1, col = rgb(0.7, 0.1, 0.1, 1), lwd = 2)
# lines(1:nrow(matrix_output_GIr2[4500:4740,]), median_GIr2, col = rgb(0.1, 0.1, 0.7, 1), lwd = 2)

abline(v=120,col="red",lwd=2,lty = 2)
text(85,1000,"Start Treatment")
abline(v=180,col="purple3",lwd=2,lty = 2)
text(160,7000,"5 year")

legend("topright", legend = c("No treatment", "Placebo"),
       col = c(rgb(0.1, 0.7, 0.1, 1), rgb(0.7, 0.1, 0.1, 1)),
       lwd = 2, lty = 1, bg = "white",cex = 1)

#### gametocyte resistance (Asym/Sym) ######
median_Gr_sym <- apply(matrix_output_Gr_sym[4500:4740,], 1, median, na.rm = TRUE)
lower_Gr_sym <- apply(matrix_output_Gr_sym[4500:4740,], 1, quantile, probs = 0.025, na.rm = TRUE)
upper_Gr_sym <- apply(matrix_output_Gr_sym[4500:4740,], 1, quantile, probs = 0.975, na.rm = TRUE)

median_Gr_asym <- apply(matrix_output_Gr_asym[4500:4740,], 1, median, na.rm = TRUE)
lower_Gr_asym <- apply(matrix_output_Gr_asym[4500:4740,], 1, quantile, probs = 0.025, na.rm = TRUE)
upper_Gr_asym <- apply(matrix_output_Gr_asym[4500:4740,], 1, quantile, probs = 0.975, na.rm = TRUE)

# Plot with confidence interval
plot(1:nrow(matrix_output_Gr_asym[4500:4740,]), median_Gr_asym, type = "l", col = "green", lwd = 1.5,ylim = c(0,25000),
     ylab = "Number of People Carrying Gametocytes", xlab = "Months", main = "Gametocyte resistance(Asym/Sym) with 95% CI")

# Add shaded region for confidence interval
polygon(c(1:nrow(matrix_output_Gr_asym[4500:4740,]), rev(1:nrow(matrix_output_Gr_asym[4500:4740,]))), c(lower_Gr_asym, rev(upper_Gr_asym)), col = rgb(0.2, 0.5, 0.2, 0.25), border = NA)
lines(1:nrow(matrix_output_Gr_asym[4500:4740,]),lower_Gr_asym, col="green4",lty=2)
lines(1:nrow(matrix_output_Gr_asym[4500:4740,]),upper_Gr_asym, col="green4",lty=2)

polygon(c(1:nrow(matrix_output_Gr_sym[4500:4740,]), rev(1:nrow(matrix_output_Gr_sym[4500:4740,]))), c(lower_Gr_sym, rev(upper_Gr_sym)), col = rgb(0.2, 0.2, 0.5, 0.25), border = NA)
lines(1:nrow(matrix_output_Gr_sym[4500:4740,]),lower_Gr_sym, col="blue",lty=2)
lines(1:nrow(matrix_output_Gr_sym[4500:4740,]),upper_Gr_sym, col="blue",lty=2)

lines(1:nrow(matrix_output_Gr_asym[4500:4740,]), median_Gr_asym, col = "green4", lwd = 2)
lines(1:nrow(matrix_output_Gr_sym[4500:4740,]), median_Gr_sym, col = "blue", lwd = 2)

abline(v=120,col="red",lwd=2,lty = 2)
text(85,15000,"Start Treatment")
abline(v=180,col="purple3",lwd=2,lty = 2)
text(160,20000,"5 year")

legend("topright",c("Asymptomatic","Symptomatic"),col =c("green4","blue"),lty =c(1,1), lwd = 2,cex = 1)


#### total g/total inf(IS+AS+IR+AR) * 100% ######
median_G_per_inf <- apply((matrix_output_total_G[4500:4740,]/matrix_output_total_inf[4500:4740,]), 1, median, na.rm = TRUE)
lower_G_per_inf  <- apply((matrix_output_total_G[4500:4740,]/matrix_output_total_inf[4500:4740,]), 1, quantile, probs = 0.025, na.rm = TRUE)
upper_G_per_inf  <- apply((matrix_output_total_G[4500:4740,]/matrix_output_total_inf[4500:4740,]), 1, quantile, probs = 0.975, na.rm = TRUE)


# Plot with confidence interval
plot(1:nrow(matrix_output_inc[4500:4740,]), median_G_per_inf, type = "l", col = "yellow", lwd = 2, ylim = range(lower_G_per_inf, upper_G_per_inf),
     ylab = "Ratio", xlab = "Months", main = "Total Carrying Gametocytes/total Infected")

# Add shaded region for confidence interval
polygon(c(1:nrow(matrix_output_inc[4500:4740,]), rev(1:nrow(matrix_output_inc[4500:4740,]))), c(lower_G_per_inf, rev(upper_G_per_inf)), col = rgb(0.75, 0.75, 0.75, 0.5), border = NA)
lines(1:nrow(matrix_output_inc[4500:4740,]),lower_G_per_inf, col="grey75",lty=2)
lines(1:nrow(matrix_output_inc[4500:4740,]),upper_G_per_inf, col="grey75",lty=2)

# Add median line
lines(1:nrow(matrix_output_inc[4500:4740,]), median_G_per_inf, col = "yellow3", lwd = 2)
#### total gr/total gs ######
median_Gr_per_Gs <- apply((matrix_output_Gs[4500:4740,]/matrix_output_Gr[4500:4740,]), 1, median, na.rm = TRUE)
lower_Gr_per_Gs  <- apply((matrix_output_Gs[4500:4740,]/matrix_output_Gr[4500:4740,]), 1, quantile, probs = 0.025, na.rm = TRUE)
upper_Gr_per_Gs  <- apply((matrix_output_Gs[4500:4740,]/matrix_output_Gr[4500:4740,]), 1, quantile, probs = 0.975, na.rm = TRUE)


# Plot with confidence interval
plot(1:nrow(matrix_output_inc[4500:4740,]), median_Gr_per_Gs, type = "l", col = "orange", lwd = 2, ylim = range(lower_Gr_per_Gs, upper_Gr_per_Gs),
     ylab = "Ratio", xlab = "Months", main = "Total Carrying GR/Carrying GS")

# Add shaded region for confidence interval
polygon(c(1:nrow(matrix_output_inc[4500:4740,]), rev(1:nrow(matrix_output_inc[4500:4740,]))), c(lower_Gr_per_Gs, rev(upper_Gr_per_Gs)), col = rgb(0.75, 0.75, 0.75, 0.5), border = NA)
lines(1:nrow(matrix_output_inc[4500:4740,]),lower_Gr_per_Gs, col="grey75",lty=2)
lines(1:nrow(matrix_output_inc[4500:4740,]),upper_Gr_per_Gs, col="grey75",lty=2)

# Add median line
lines(1:nrow(matrix_output_inc[4500:4740,]), median_Gr_per_Gs, col = "orange", lwd = 2)

#### total gsym/total gasym ######

median_Gsym_per_Gasym <- apply((matrix_output_Gsym[4500:4740,]/matrix_output_Gasym[4500:4740,]), 1, median, na.rm = TRUE)
lower_Gsym_per_Gasym  <- apply((matrix_output_Gsym[4500:4740,]/matrix_output_Gasym[4500:4740,]), 1, quantile, probs = 0.025, na.rm = TRUE)
upper_Gsym_per_Gasym  <- apply((matrix_output_Gsym[4500:4740,]/matrix_output_Gasym[4500:4740,]), 1, quantile, probs = 0.975, na.rm = TRUE)


# Plot with confidence interval
plot(1:nrow(matrix_output_inc[4500:4740,]), median_Gsym_per_Gasym, type = "l", col = "purple2", lwd = 2, ylim = range(lower_Gsym_per_Gasym, upper_Gsym_per_Gasym),
     ylab = "Ratio", xlab = "Months", main = "Total Symptomatic Carrying Gametocytes Per\n Total Asymptomatic Carrying Gametocytes")

# Add shaded region for confidence interval
polygon(c(1:nrow(matrix_output_inc[4500:4740,]), rev(1:nrow(matrix_output_inc[4500:4740,]))), c(lower_Gsym_per_Gasym, rev(upper_Gsym_per_Gasym)), col = rgb(0.75, 0.75, 0.75, 0.5), border = NA)
lines(1:nrow(matrix_output_inc[4500:4740,]),lower_Gsym_per_Gasym, col="grey75",lty=2)
lines(1:nrow(matrix_output_inc[4500:4740,]),upper_Gsym_per_Gasym, col="grey75",lty=2)

# Add median line
lines(1:nrow(matrix_output_inc[4500:4740,]), median_Gsym_per_Gasym, col = "purple2", lwd = 2)

abline(v=120,col="red",lwd=2,lty = 2)
text(85,0.2,"Start Treatment")
abline(v=180,col="blue4",lwd=2,lty = 2)
text(160,0.35,"5 year")

