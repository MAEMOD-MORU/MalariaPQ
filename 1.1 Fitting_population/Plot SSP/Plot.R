#detach all library
rm(list = ls())
library(deSolve)

ssp_year <- seq(2000, 2050, by = 5)
totalpop_data <- read.csv("../../data/Tanzania_pop_2000_2100.csv", header = TRUE)[1:51,2]
data_ssp <- read.csv("../../data/Tanzania_pop_SSP_2000_2100.csv", header = TRUE)[1:11,2:6]
year <- 2000:2050


# Plot the results
plot(year, totalpop_data, type = "o", col = "black",
     xlab = "Time (Year)", ylab = "Population", pch=19,
     main = "Total Population With Different SSPs",
     ylim = range(c(0,totalpop_data,data_ssp$SSP1, data_ssp$SSP2, data_ssp$SSP3)))
points(ssp_year, data_ssp$SSP1, col = "blue",type = "o", pch = 18)
points(ssp_year, data_ssp$SSP2, col = "red",type = "o", pch = 15)
points(ssp_year, data_ssp$SSP3, col =" green",type = "o", pch = 17)
grid()
legend("topleft", legend = c("SSP1 population (IIASA)", "SSP2 population (IIASA)",
                             "SSP3 population (IIASA)","WHO Population Estimate"),ncol = 1,
       col = c( "blue", "red", "green3","black"),
       pch = c( 18, 15, 17,16), lty = c(1, 1, 1, 1),lwd=3)


