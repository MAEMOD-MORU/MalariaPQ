library(deSolve)
#Read model
source("model_D0_33_D1_67_baseline_change_birth_death_2betaSets_SS.R")

source("init_parameter_after_1st_run_2015_D1_33_SS.R")

run_malariaPQ_ode <-function(
             time =seq(1, 12*51, by = 1),
             prob_sym_s = 0.25,
             prob_sym_r = 0.25,
             mui = 0.0185585066037,
             muo = 0.0165193099545,
             mui_before = 0.0185585066037,
             muo_before = 0.0165193099545,
             flo = 0.9,
             amp = 0.1,
             period = 12,
             phase = 0,
             propIs= c(0.67, 0.33 ,0),
             propIr= c(0.67, 0.33 ,0),
             propIs_new = c(0.67, 0 ,0.33),
             propIr_new = c(0.67, 0 ,0.33),
             start_d = 6000,
             start_m = 12*1000,
             start_b = 12*1000,
             t_long = 12*100,
             Fail_rate_s = 0.025,
             Fail_rate_r = 0.13,
             beta_s = c(0.671, 0.671, 0.671),
             beta_as = 5.992,
             m_beta = 1.05,
             beta_s_2 = c(0.671, 0.671, 0.671),
             beta_as_2 = 5.992,
             d0 = (1 / 90) * 30,
             d1 = (1 / 44.1) * 30,
             d2 = (1 / 1) * 30,
             c_beta_r = 0.745,
             init_state = c(
               S      = 4.664516e+07,
               Is0    = 1.472014e+05,
               Is1    = 2.952939e+04,
               Is2    = 0.000000e+00,
               Stis1  = 4.015327e+04,
               Stis2  = 0.000000e+00,
               Fis1   = 3.833547e+02,
               Fis2   = 0.000000e+00,
               GIs0   = 5.349611e+04,
               GIs1   = 6.925867e+03,
               GIs2   = 0.000000e+00,
               GAs    = 3.248219e+05,
               As     = 5.762437e+05,
               Rs     = 4.197051e+06,
               Ir0    = 0.000000e+00,
               Ir1    = 0.000000e+00,
               Ir2    = 0.000000e+00,
               Stir1  = 0.000000e+00,
               Stir2  = 0.000000e+00,
               Fir1   = 0.000000e+00,
               Fir2   = 0.000000e+00,
               GIr0   = 0.000000e+00,
               GIr1   = 0.000000e+00,
               GIr2   = 0.000000e+00,
               GAr    = 0.000000e+00,
               Ar     = 0.000000e+00,
               Rr     = 0.000000e+00
             ),
             m_init = 0.092
             ,
             events = list(
               func = function(time, state, parameters) {
                 if (time ==  1) {
                   state["Ir0"] <- state["Ir0"] + 1
                   state["Ir1"] <- state["Ir1"] + 1
                   # state["Ir2"] <- state["Ir2"] + 1
                   state["Ar"]  <- state["Ar"]  + 1
                   
                   # reduce Is
                   state["Is0"] <- state["Is0"] - 1
                   state["Is1"] <- state["Is1"] - 1
                   # state["Is2"] <- state["Is2"] - 1
                   state["As"]  <- state["As"]  - 1
                 }
                 return(state)
               },
               time =  1  # only one time point
             )
             ) {
  parameters <- list(
    mui = mui,
    muo = muo,
    mui_before = mui_before,
    muo_before = muo_before,
    flo = flo,
    amp = amp,
    propIs=propIs,
    propIr=propIr,
    propIs_new = propIs_new,
    propIr_new = propIr_new,
    period = period,
    phase = phase,
    prob_sym_s = prob_sym_s,
    prob_sym_r = prob_sym_r,
    start_d = start_d,# 10 years after 3500 months
    start_m = start_m, # 2015
    start_b = start_b, # 2015
    t_long = t_long, # 2,5,10 years 
    Fail_rate_s = Fail_rate_s, # Fail Treatment of ACT
    Fail_rate_r = Fail_rate_r,
    
    beta_s = beta_s,
    beta_r = beta_s*m_beta,
    beta_as = beta_as,
    beta_ar = beta_as*m_beta,
    beta_s_2 =  beta_s_2,
    beta_as_2 = beta_as_2,
    
    alpha = 1.07 / 12, #  doi:10.1371/journal.pone.0001767
    
    τigs=c(0.33,0.68,6.15),
    Tas=0.33,
    Tags=0.33,
    
    τigr=c(0.33,0.68,6.15),
    Tar=0.33,
    Tagr=0.33,
    
    τntsd0= 0.33,
    τntrd0= 0.33,
    
    τnfs = rep(4.3,2) ,
    τfs = rep(0.33,2) ,
    
    τnfr=rep(4.3,2),
    τfr=rep(0.33,2),
    
    d0 = d0, # no drugs
    d1 = d1, # Gametocyte recovery rate from Palang paper, table 2 at row Placebo (DP), column microscopy
    d2 = d2, # Gametocyte recovery rate from Palang paper, table 2 at row 0.0625 mg/kg, column microscopy
    
    g_is = c(0.121,0.121,0.121), #0.1
    g_ir = c(0.242,0.242,0.242),
    g_ar = 0.372,
    g_as = 0.186,
    g_infs = c(0.019,0.019),
    g_ifs = c(1.21,1.21),  #1
    g_infr = c(0.038,0.038),
    g_ifr = c(2.42,2.42),
    
    c_beta_r = c_beta_r
  )
  out <- ode(
    y     = init_state * m_init,
    times = time,
    func  = Malaria_model_with_Array,
    parms = parameters,
    events = events
  )
  return(as.data.frame(out))
}

prep_malaria_outputs <- function(out,
                                 start_year = 2000,
                                 end_year   = 2050) {
  # Basic checks
  needed <- c("inc","inc_r","N","inc_sym","inc_sym_r","inc_asym","inc_asym_r")
  miss <- setdiff(needed, colnames(out))
  if (length(miss) > 0) stop("Missing columns in out: ", paste(miss, collapse=", "))
  
  n_months <- end_year - start_year + 1L
  if (nrow(out) < 12L * n_months) {
    stop("out has ", nrow(out), " rows; need at least ", 12L * n_months,
         " months for years ", start_year, ":", end_year)
  }
  
  # Helper: reshape a monthly vector into 12 x years
  to_12xY <- function(x) matrix(x[1:(12L*n_months)], nrow = 12L, byrow = FALSE)
  
  # Annual totals/means
  mat_inc   <- to_12xY(out[,"inc"])
  inc_tot   <- colSums(mat_inc)
  
  mat_incr  <- to_12xY(out[,"inc_r"])
  inc_r_tot <- colSums(mat_incr)
  
  mat_N     <- to_12xY(out[,"N"])
  N_mean    <- colMeans(mat_N)
  
  mat_sym_r <- to_12xY(out[,"inc_sym_r"])
  sym_r_tot <- colSums(mat_sym_r)
  
  mat_sym   <- to_12xY(out[,"inc_sym"])
  sym_tot   <- colSums(mat_sym)
  
  mat_asym_r <- to_12xY(out[,"inc_asym_r"])
  asym_r_tot <- colSums(mat_asym_r)
  
  mat_asym   <- to_12xY(out[,"inc_asym"])
  asym_tot   <- colSums(mat_asym)
  
  years <- start_year:end_year
  tiny  <- 1e-12
  
  # Ratios
  ratio_r_total   <- inc_r_tot / (inc_tot + tiny)
  ratio_sym_r_tot <- sym_r_tot / (inc_tot + tiny)
  ratio_sym_r_sym <- sym_r_tot / (sym_tot + tiny)
  ratio_asym_r_tot <- asym_r_tot / (inc_tot + tiny)
  ratio_asym_r_asym <- asym_r_tot / (asym_tot + tiny)
  
  list(
    years = years,
    annual = list(
      inc_total = inc_tot,
      inc_r     = inc_r_tot,
      N_mean    = N_mean,
      inc_sym   = sym_tot,
      inc_sym_r = sym_r_tot,
      inc_asym  = asym_tot,
      inc_asym_r= asym_r_tot
    ),
    ratio = list(
      r_over_total      = ratio_r_total,
      sym_r_over_total  = ratio_sym_r_tot,
      sym_r_over_sym    = ratio_sym_r_sym,
      asym_r_over_total = ratio_asym_r_tot,
      asym_r_over_asym  = ratio_asym_r_asym
    )
  )
}

plot_malaria_outputs <- function(prep,
                                 plot_type = c("Total Inc + Res",
                                               "Symptomatic",
                                               "Asymptomatic",
                                               "Sym and Asym"),
                                 year_plot = 2000,
                                 Tanzania_Incidence = NULL,
                                 Tanzania_k13_Allele_frequency = NULL,
                                 ref_year = 2025,
                                 ref_point = 0.4) {
  
  plot_type <- match.arg(plot_type)
  years <- prep$years
  A <- prep$annual
  R <- prep$ratio
  
  # helper for points if provided
  add_k13_points <- function() {
    if (!is.null(Tanzania_k13_Allele_frequency)) {
      points(c(2016:2022, ref_year),
             c(Tanzania_k13_Allele_frequency, ref_point),
             col = 2, lwd = 2, pch = 3)
      abline(v = ref_year, lty = 2)
    }
  }
  
  if (plot_type == "Total Inc + Res") {
    par(mfrow = c(2, 1))
    
    plot(years, A$inc_total,
         xlab="Year", ylab="Incidence",
         lwd=2, type="l",
         main="Total Incidence",
         xlim=c(year_plot, max(years)),
         ylim=c(0, max(A$inc_total)*1.1))
    
    lines(years, A$inc_r, col="red")
    lines(years, A$N_mean, col="green3")
    
    if (!is.null(Tanzania_Incidence)) {
      points(2015:2023, Tanzania_Incidence[,2], col="blue", lwd=2)
    }
    
    legend("bottomright",
           legend=c("Total Incidence", "Resistant Incidence", "Population", "Tanzania Reported Incidence"),
           col=c("black","red","green3","blue"),
           lwd=c(2,2,2,NA),
           pch=c(NA,NA,NA,1),
           cex=1)
    
    plot(years, R$r_over_total,
         type="l", col="red",
         xlab="Year", ylab="Ratio",
         xlim=c(year_plot, max(years)),
         main="Ratio Resistant Incidence / Total Incidence")
    
    add_k13_points()
    
    # annotate at ref year if in range
    if (ref_year %in% years) {
      idx <- which(years == ref_year)
      text(ref_year + 2, ref_point, round(R$r_over_total[idx], 2))
    }
    
  } else if (plot_type == "Symptomatic") {
    par(mfrow = c(2, 1))
    
    plot(years, R$sym_r_over_total,
         type="l", col="red",
         xlab="Year", ylab="Ratio",
         main="Symptomatic Resistant Incidence / Total Incidence",
         ylim=c(0,0.5),
         xlim=c(year_plot, max(years)))
    add_k13_points()
    
    plot(years, R$sym_r_over_sym,
         type="l", col="red",
         xlab="Year", ylab="Ratio",
         main="Symptomatic Resistant Incidence / Symptomatic Incidence",
         ylim=c(0,0.5),
         xlim=c(year_plot, max(years)))
    add_k13_points()
    
  } else if (plot_type == "Asymptomatic") {
    par(mfrow = c(2, 1))
    
    plot(years, R$asym_r_over_total,
         type="l", col="red",
         xlab="Year", ylab="Ratio",
         main="Asymptomatic Resistant Incidence / Total Incidence",
         ylim=c(0,0.5),
         xlim=c(year_plot, max(years)))
    add_k13_points()
    
    plot(years, R$asym_r_over_asym,
         type="l", col="red",
         xlab="Year", ylab="Ratio",
         main="Asymptomatic Resistant Incidence / Asymptomatic Incidence",
         ylim=c(0,0.5),
         xlim=c(year_plot, max(years)))
    add_k13_points()
    
  } else if (plot_type == "Sym and Asym") {
    par(mfrow = c(2, 1))
    
    plot(years, R$sym_r_over_sym,
         type="l", col="red",
         xlab="Year", ylab="Ratio",
         main="Symptomatic Resistant Incidence / Symptomatic Incidence",
         ylim=c(0,0.5),
         xlim=c(year_plot, max(years)))
    add_k13_points()
    
    plot(years, R$asym_r_over_asym,
         type="l", col="red",
         xlab="Year", ylab="Ratio",
         main="Asymptomatic Resistant Incidence / Asymptomatic Incidence",
         ylim=c(0,0.5),
         xlim=c(year_plot, max(years)))
    add_k13_points()
  }
  
  invisible(prep)
}

