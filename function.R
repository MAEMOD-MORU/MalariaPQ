run_matrix_output <- function(n_run = 1,time = 1,select = "inc"){
  times <- seq(0, time, by = 1)
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
    matrix_output[,i] <- out[,select]
  }

  return(matrix_output)
}

plot_out_line <- function(out,...){
  plot(out[,1],type = "l",...)
  for (i in 2:length(out[1,])) {
    lines(out[,i],col=i)
  }
}

plot_out_CI95 <- function(out,...){
  # Compute median, lower (2.5th), and upper (97.5th) percentiles
  median_inc <- apply(out, 1, median, na.rm = TRUE)
  lower_inc <- apply(out, 1, quantile, probs = 0.025, na.rm = TRUE)
  upper_inc <- apply(out, 1, quantile, probs = 0.975, na.rm = TRUE)

  # Plot with confidence interval
  plot(1:nrow(out), median_inc, type = "l", col = "blue", lwd = 1.5,...)

  # Add shaded region for confidence interval
  polygon(c(1:nrow(out), rev(1:nrow(out))), c(lower_inc, rev(upper_inc)), col = rgb(0.5, 0.5, 0.5, 0.75), border = 1)
  # lines(1:nrow(out),lower_inc, col="grey50",lty=2)
  # lines(1:nrow(out),upper_inc, col="grey50",lty=2)

  # Add median line
  lines(1:nrow(out), median_inc, col = "blue", lwd = 1.5)
}
