# Run Simulation
library(deSolve)

#Read model
source("model_D0_33_D1_67_baseline.r")

source("function.r")

source("init_parameter.r")

out<- ode(y = state, times = times, func = Malaria_model_with_Array_2, parms = parameters)
out_5<- ode(y = state, times = times, func = Malaria_model_with_Array_2, parms = parameters_5)
out_2<- ode(y = state, times = times, func = Malaria_model_with_Array_2, parms = parameters_2)

# 10 year
plot(out[(4500):(4500+360),"inc"],type="l",xlab="Months",ylab="Incidence",main="Total Incidence (10 year treatment)",col=1,ylim = c(0,1.6e5))
# lines(out[(3500):(3500+360),"inc_r"],type="l",col=2)
# legend("topright",c("Sensitive","Resistant"),col = 1:2,lty =c(1,1),cex = 0.9)
abline(v=120,col="green")
text(170,155000,"start treatment")
abline(v=240,col="blue")
text(300,155000,"Stop treatment")

# 5 year
plot(out_5[(4500):(4500+300),"inc"],type="l",xlab="Months",ylab="Incidence",main="Total Incidence (5 year treatment)",col=1,ylim = c(0,1.6e5))
# lines(out[(3500):(3500+360),"inc_r"],type="l",col=2)
# legend("topright",c("Sensitive","Resistant"),col = 1:2,lty =c(1,1),cex = 0.9)
abline(v=120,col="green")
text(150,155000,"start\n treatment")
abline(v=180,col="blue")
text(230,155000,"Stop treatment")

# 2 year
plot(out_2[(4500):(4500+300-36),"inc"],type="l",xlab="Months",ylab="Incidence",main="Total Incidence (2 year treatment)",col=1,ylim = c(0,1.6e5))
# lin# lin# lines(out[(3500):(3500+360),"inc_r"],type="l",col=2)
# legend("topright",c("Sensitive","Resistant"),col = 1:2,lty =c(1,1),cex = 0.9)
abline(v=120,col="green")
text(90,100000,"start\n treatment")
abline(v=180-36,col="blue")
text(180,155000,"Stop treatment")

# 10 year
plot(out[(4500):(4500+360),"inc"],type="l",xlab="Months",ylab="Incidence",main="Total Incidence (2,5,10 year treatment)",col=1,ylim = c(0,1.6e5))
# 5 year
lines(out_5[(4500):(4500+360),"inc"],type="l",col=2)
# 2 year
lines(out_2[(4500):(4500+360),"inc"],type="l",col=3)

abline(v=120,col="blue")
text(85,50000,"start\n treatment")
abline(v=180-36,col="green4")
text(160,145000,"2\n years")
abline(v=180,col="red3")
text(200,145000,"5\n years")
abline(v=240,col="black")
text(260,145000,"10\n years")

