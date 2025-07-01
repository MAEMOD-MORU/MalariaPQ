# Run Simulation
library(deSolve)

#Read model
source("model_D0_33_D1_67_baseline.r")

source("function.r")

source("init_parameter.r")

out<- ode(y = state, times = times, func = Malaria_model_with_Array_2, parms = parameters)

# total incidence
plot(out[0:4500,"inc"],type="l",xlab="Months",ylab="Incidence",main="Total Incidence (Baseline D1 66%)")
plot(0:240,out[(4500-240):4500,"inc"],type="l",xlab="Months",ylab="Incidence",main="Total Incidence (Baseline D1 66%)")

#check total population
plot(out[,"N"],type="l",xlab="Months",ylab="Total")
plot(out[,"N"],type="l",xlab="Months",ylab="Population",main="Total Population", ylim = c(70000050,70000050))

plot(out[(4500-240):4500,"inc_s"],type="l",xlab="Months",ylab="Incidence",main="Incidence Sensitive/Resistant",col=1,ylim = c(0,1.2e5))
lines(out[(4500-240):4500,"inc_r"],type="l",xlab="Months",ylab="Incidence",main="Incidence Resistant",col=2)
legend("bottomright",c("Sensitive","Resistant"),col = 1:2,lty =c(1,1))

plot(out[(4500-240):4500,"inc_asym"],type="l",xlab="Months",ylab="Incidence",main="Incidence Asymptomatic/Symptomatic",col=3,ylim = c(0,1.2e5))
lines(out[(4500-240):4500,"inc_sym"],type="l",col=4)
legend("bottomright",c("Asymptomatic","Symptomatic"),col = 3:4,lty =c(1,1),cex = 0.9)

plot(out[(4500-240):4500,"GS_inc"],type="l",xlab="Months",ylab="Incidence",main="Gametocyte Infections Sensitive/Resistant",col=1,ylim = c(0,4e4))
lines(out[(4500-240):4500,"GR_inc"],type="l",col=2)
legend("bottomright",c("Sensitive","Resistant"),col = 1:2,lty =c(1,1))

plot(out[(4500-240):4500,"Gasym_inc"],type="l",xlab="Months",ylab="Incidence",main="Gametocyte Infections Sensitive/Resistant",col=3,ylim = c(0,5e4),lwd=2)
lines(out[(4500-240):4500,"Gsym_inc"],type="l",col=4,lwd=2)
legend("bottomright",c("Asymptomatic","Symptomatic"),col = 3:4,lty =c(1,1))

plot(rowSums(out[(4500-240):4500,c("Fir1","Fir2")]),type="l",xlab="Months",ylab="Patient",main="Fail to Treatment Patient",col=2,ylim = c(0,1e3))
lines(rowSums(out[(4500-240):4500,c("Fis1","Fis2")]),type="l",col=1)
legend("bottomright",c("Sensitive","Resistant"),col = 1:2,lty =c(1,1))

