# Run Simulation
library(deSolve)

#Read model
source("model_D0_33_D1_67_baseline.r")

#
source("init_parameter.r")

#
source("function.r")

#matrix
n_run <- 2

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
  # rd1 <- runif(n=1,min = 10.9,max = 75.6) #44.1
  # rd2 <- runif(n=1,min = 1.37,max = 11.0) #4.88
  # # parameters$d0 <- 1/rd0*30
  # parameters$d1 <- 1/rd1*30
  # parameters$d2 <- 1/rd2*30
  
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



plot(out[4500:4740,"inc_s"],type="l",xlab="Months",ylab="Incidence",
     main="Incidence Sensitive/Resistant",col=1,ylim = c(0,1.2e5))
# lines(out_test[,"inc_s"],type="l",col=2)
lines(out[4500:4740,"inc_r"],type="l",col=2)
abline(v=120,col="red",lwd=2,lty = 2)
text(85,60000,"Start Treatment")
abline(v=180,col="purple3",lwd=2,lty = 2)
text(160,50000,"5 year")
legend("topright",c("Sensitive","Resistant"),col = 1:2,lty =c(1,1),lwd=3)

plot(out[4500:4740,"inc_asym"],type="l",xlab="Months",ylab="Incidence",main="Incidence Asymptomatic/Symptomatic",col=3,ylim = c(0,1.2e5))
lines(out[4500:4740,"inc_sym"],type="l",col=4)
abline(v=120,col="red",lwd=2,lty = 2)
text(85,60000,"Start Treatment")
abline(v=180,col="purple3",lwd=2,lty = 2)
text(160,50000,"5 year")
legend("topright",c("Asymptomatic","Symptomatic"),col = 3:4,lty =c(1,1),cex = 0.9,lwd=3)


plot(out[4500:4740,"GS_inc"],type="l",xlab="Months",ylab="Incidence",main="Gametocyte Infections Sensitive/Resistant",col=1,ylim = c(0,4e4))
lines(out[4500:4740,"GR_inc"],type="l",col=2)
abline(v=120,col="red",lwd=2,lty = 2)
text(85,40000,"Start Treatment")
abline(v=180,col="purple3",lwd=2,lty = 2)
text(160,30000,"5 year")
legend("topright",c("Sensitive","Resistant"),col = 1:2,lty =c(1,1),lwd=3)

plot(out[4500:4740,"Gasym_inc"],type="l",xlab="Months",ylab="Incidence",main="Gametocyte Infections Sensitive/Resistant",col=3,ylim = c(0,5e4),lwd=2)
lines(out[4500:4740,"Gsym_inc"],type="l",col=4,lwd=2)
abline(v=120,col="red",lwd=2,lty = 2)
text(85,30000,"Start Treatment")
abline(v=180,col="purple3",lwd=2,lty = 2)
text(160,20000,"5 year")
legend("topright",c("Asymptomatic","Symptomatic"),col = 3:4,lty =c(1,1),lwd=3)

plot(out[4500:4740,"Gasym_inc"],type="l",xlab="Months",ylab="Incidence",main="Gametocyte Infections Sensitive/Resistant",col=3,ylim = c(0,5e4),lwd=2)
lines(out[4500:4740,"Gsym_inc"],type="l",col=4,lwd=2)
abline(v=120,col="red",lwd=2,lty = 2)
text(85,30000,"Start Treatment")
abline(v=180,col="purple3",lwd=2,lty = 2)
text(160,20000,"5 year")
legend("topright",c("Asymptomatic","Symptomatic"),col = 3:4,lty =c(1,1),lwd=3)

plot(out[4500:4740,"GAs"],type="l",xlab="Months",ylab="people",main="Gametocyte Sensitive Compartment (Asymptomatic/Symptomatic)",col=3,ylim = c(0,9e4),lwd=2)
lines(rowSums(out[4500:4740,c("GIs0","GIs1","GIs2")]),type="l",col=4,lwd=2)
abline(v=120,col="red",lwd=2,lty = 2)
text(85,30000,"Start Treatment")
abline(v=180,col="purple3",lwd=2,lty = 2)
text(160,20000,"5 year")
legend("topright",c("Asymptomatic","Symptomatic"),col = 3:4,lty =c(1,1),lwd=3)


plot(out[4500:4740,"GAr"],type="l",xlab="Months",ylab="people",main="Gametocyte Resistant Compartment (Asymptomatic/Symptomatic)",col=3,ylim = c(0,6e4),lwd=2)
lines(rowSums(out[4500:4740,c("GIr0","GIr1","GIr2")]),type="l",col=4,lwd=2)
abline(v=120,col="red",lwd=2,lty = 2)
text(85,30000,"Start Treatment")
abline(v=180,col="purple3",lwd=2,lty = 2)
text(160,20000,"5 year")
legend("topright",c("Asymptomatic","Symptomatic"),col = 3:4,lty =c(1,1),lwd=3)
