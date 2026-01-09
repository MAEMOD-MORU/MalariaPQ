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
             flo = 0.9,
             amp = 0.1,
             period = 12,
             phase = 0,
             propIs= c(0.67, 0.33 ,0),
             propIr= c(0.67, 0.33 ,0),
             propIs_new = c(0.67, 0 ,0.33),
             propIr_new = c(0.67, 0 ,0.33),
             start_d = 6000,
             t_long = 12*10,
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
             # multiplicative factor for initial state (population)
             m_init = 0.092
             ,
             # for adding initial resistant infections
             # time point and function to modify state
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
    start_d = start_d,# 10 years
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
  
  # add all output form modeling into a list
  # names(out)
  # [1] "time"            "S"               "Is0"             "Is1"             "Is2"            
  # [6] "Stis1"           "Stis2"           "Fis1"            "Fis2"            "GIs0"           
  # [11] "GIs1"            "GIs2"            "GAs"             "As"              "Rs"             
  # [16] "Ir0"             "Ir1"             "Ir2"             "Stir1"           "Stir2"          
  # [21] "Fir1"            "Fir2"            "GIr0"            "GIr1"            "GIr2"           
  # [26] "GAr"             "Ar"              "Rr"              "inc"             "inc_s"          
  # [31] "inc_r"           "inc_sym"         "inc_sym_s"       "inc_sym_r"       "inc_asym"       
  # [36] "inc_asym_s"      "inc_asym_r"      "G_inc"           "GS_inc"          "GR_inc"         
  # [41] "Gsym_inc"        "Gsym_inc_s"      "Gsym_inc_r"      "Gasym_inc"       "Gasym_inc_s"    
  # [46] "Gasym_inc_r"     "N"               "F_T"             "CT_total"        "CU_total"       
  # [51] "p_resistant_inc"
  
  # Annual summaries
  mat_S     <- to_12xY(out[,"S"])
  S_total   <- colSums(mat_S)
  
  mat_Is0    <- to_12xY(out[,"Is0"])
  Is0_total  <- colSums(mat_Is0)
  
  mat_Is1    <- to_12xY(out[,"Is1"])
  Is1_total  <- colSums(mat_Is1)
  
  mat_Is2    <- to_12xY(out[,"Is2"])
  Is2_total  <- colSums(mat_Is2)
  
  mat_Is    <- to_12xY(out[,"Is0"] + out[,"Is1"] + out[,"Is2"])
  Is_total  <- colSums(mat_Is)
  
  matStis1    <- to_12xY(out[,"Stis1"])
  Stis1_total  <- colSums(matStis1)
  
  matStis2    <- to_12xY(out[,"Stis2"])
  Stis2_total  <- colSums(matStis2)
  
  matStis    <- to_12xY(out[,"Stis1"] + out[,"Stis2"])
  Stis_total  <- colSums(matStis)
  
  matFis1    <- to_12xY(out[,"Fis1"])
  Fis1_total  <- colSums(matFis1)
  
  matFis2    <- to_12xY(out[,"Fis2"])
  Fis2_total  <- colSums(matFis2)
  
  matFis    <- to_12xY(out[,"Fis1"] + out[,"Fis2"])
  Fis_total  <- colSums(matFis)
  
  mat_GIs0    <- to_12xY(out[,"GIs0"])
  GIs0_total  <- colSums(mat_GIs0)
  
  mat_GIs1    <- to_12xY(out[,"GIs1"])
  GIs1_total  <- colSums(mat_GIs1)
  
  mat_GIs2    <- to_12xY(out[,"GIs2"])
  GIs2_total  <- colSums(mat_GIs2)
  
  mat_GIs    <- to_12xY(out[,"GIs0"] + out[,"GIs1"] + out[,"GIs2"])
  GIs_total  <- colSums(mat_GIs)
  
  mat_GAs    <- to_12xY(out[,"GAs"])
  GAs_total  <- colSums(mat_GAs)
  
  mat_As    <- to_12xY(out[,"As"])
  As_total  <- colSums(mat_As)
  
  mat_Rs    <- to_12xY(out[,"Rs"])
  Rs_total  <- colSums(mat_Rs)
  
  mat_Ir0    <- to_12xY(out[,"Ir0"])
  Ir0_total  <- colSums(mat_Ir0)
  
  mat_Ir1    <- to_12xY(out[,"Ir1"])
  Ir1_total  <- colSums(mat_Ir1)
  
  mat_Ir2    <- to_12xY(out[,"Ir2"])
  Ir2_total  <- colSums(mat_Ir2)
  
  mat_Ir    <- to_12xY(out[,"Ir0"] + out[,"Ir1"] + out[,"Ir2"])
  Ir_total  <- colSums(mat_Ir)
  
  mat_Stir1    <- to_12xY(out[,"Stir1"])
  Stir1_total  <- colSums(mat_Stir1)
  
  mat_Stir2    <- to_12xY(out[,"Stir2"])
  Stir2_total  <- colSums(mat_Stir2)
  
  mat_Stir    <- to_12xY(out[,"Stir1"] + out[,"Stir2"])
  Stir_total  <- colSums(mat_Stir)
  
  mat_Fir1    <- to_12xY(out[,"Fir1"])
  Fir1_total  <- colSums(mat_Fir1)
  
  mat_Fir2    <- to_12xY(out[,"Fir2"])
  Fir2_total  <- colSums(mat_Fir2)
  
  mat_Fir    <- to_12xY(out[,"Fir1"] + out[,"Fir2"])
  Fir_total  <- colSums(mat_Fir)
  
  mat_GIr0    <- to_12xY(out[,"GIr0"])
  GIr0_total  <- colSums(mat_GIr0)
  
  mat_GIr1    <- to_12xY(out[,"GIr1"])
  GIr1_total  <- colSums(mat_GIr1)
  
  mat_GIr2    <- to_12xY(out[,"GIr2"])
  GIr2_total  <- colSums(mat_GIr2)
  
  mat_GIr    <- to_12xY(out[,"GIr0"] + out[,"GIr1"] + out[,"GIr2"])
  GIr_total  <- colSums(mat_GIr)
  
  mat_GAr    <- to_12xY(out[,"GAr"])
  GAr_total  <- colSums(mat_GAr)
  
  mat_Ar    <- to_12xY(out[,"Ar"])
  Ar_total  <- colSums(mat_Ar)
  
  mat_Rr    <- to_12xY(out[,"Rr"])
  Rr_total  <- colSums(mat_Rr)
  
  # Annual totals/means
  mat_inc   <- to_12xY(out[,"inc"])
  inc_tot   <- colSums(mat_inc)
  
  mat_inc_s   <- to_12xY(out[,"inc_s"])
  inc_s_tot   <- colSums(mat_inc_s)
  
  mat_inc_r  <- to_12xY(out[,"inc_r"])
  inc_r_tot <- colSums(mat_inc_r)
  
  mat_sym   <- to_12xY(out[,"inc_sym"])
  sym_tot   <- colSums(mat_sym)
  
  mat_sym_s   <- to_12xY(out[,"inc_sym_s"])
  sym_s_tot   <- colSums(mat_sym_s)
  
  mat_sym_r <- to_12xY(out[,"inc_sym_r"])
  sym_r_tot <- colSums(mat_sym_r)
  
  mat_asym   <- to_12xY(out[,"inc_asym"])
  asym_tot   <- colSums(mat_asym)
  
  mat_asym_s   <- to_12xY(out[,"inc_asym_s"])
  asym_s_tot   <- colSums(mat_asym_s)
  
  mat_asym_r <- to_12xY(out[,"inc_asym_r"])
  asym_r_tot <- colSums(mat_asym_r)
  
  mat_g_inc   <- to_12xY(out[,"G_inc"])
  g_inc_tot   <- colSums(mat_g_inc)
  
  mat_gs_inc   <- to_12xY(out[,"GS_inc"])
  gs_inc_tot   <- colSums(mat_gs_inc)
  
  mat_gr_inc <- to_12xY(out[,"GR_inc"])
  gr_inc_tot <- colSums(mat_gr_inc)
  
  mat_gsym_inc   <- to_12xY(out[,"Gsym_inc"])
  gsym_inc_tot   <- colSums(mat_gsym_inc)
  
  mat_gsym_s_inc   <- to_12xY(out[,"Gsym_inc_s"])
  gsym_s_inc_tot   <- colSums(mat_gsym_s_inc)
  
  mat_gsym_r_inc <- to_12xY(out[,"Gsym_inc_r"])
  gsym_r_inc_tot <- colSums(mat_gsym_r_inc)
  
  mat_gasym_inc   <- to_12xY(out[,"Gasym_inc"])
  gasym_inc_tot   <- colSums(mat_gasym_inc)
  
  mat_gasym_s_inc   <- to_12xY(out[,"Gasym_inc_s"])
  gasym_s_inc_tot   <- colSums(mat_gasym_s_inc)
  
  mat_gasym_r_inc <- to_12xY(out[,"Gasym_inc_r"])
  gasym_r_inc_tot <- colSums(mat_gasym_r_inc)
  
  mat_N     <- to_12xY(out[,"N"])
  Total_N    <- colSums(mat_N)
  
  mat_f_T     <- to_12xY(out[,"F_T"])
  f_T_mean    <- colMeans(mat_f_T)
  
  mat_cT_total     <- to_12xY(out[,"CT_total"])
  cT_total_mean    <- colMeans(mat_cT_total)
  
  mat_cU_total     <- to_12xY(out[,"CU_total"])
  cU_total_mean    <- colMeans(mat_cU_total)
  
  mat_p_resistant_inc     <- to_12xY(out[,"p_resistant_inc"])
  p_resistant_inc_mean    <- colMeans(mat_p_resistant_inc)
  
  
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
      S          = S_total,
      Is         = Is_total,
      Stis       = Stis_total,
      Fis        = Fis_total,
      GIs        = GIs_total,
      GAs        = GAs_total,
      As         = As_total,
      Rs         = Rs_total,
      Ir         = Ir_total,
      Stir       = Stir_total,
      Fir        = Fir_total,
      GIr        = GIr_total,
      GAr        = GAr_total,
      Ar         = Ar_total,
      Rr         = Rr_total,
      Total_people = Total_N,
      inc       = inc_tot,
      inc_r     = inc_r_tot,
      inc_s     = inc_s_tot,
      inc_sym   = sym_tot,
      inc_sym_s = sym_s_tot,
      inc_sym_r = sym_r_tot,
      inc_asym   = asym_tot,
      inc_asym_s = asym_s_tot,
      inc_asym_r = asym_r_tot,
      g_inc     = g_inc_tot,
      gs_inc    = gs_inc_tot,
      gr_inc    = gr_inc_tot,
      gsym_inc  = gsym_inc_tot,
      gsym_s_inc= gsym_s_inc_tot,
      gsym_r_inc= gsym_r_inc_tot,
      gasym_inc  = gasym_inc_tot,
      gasym_s_inc= gasym_s_inc_tot,
      gasym_r_inc= gasym_r_inc_tot,
      F_T       = f_T_mean,
      CT_total  = cT_total_mean,
      CU_total  = cU_total_mean,
      p_resistant_inc = p_resistant_inc_mean
      
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
  
  add_k13_points_2 <- function() {
    if (!is.null(Tanzania_k13_Allele_frequency)) {
      points(c(2016:2022, ref_year),
             c(Tanzania_k13_Allele_frequency/2, ref_point),
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
    lines(years, A$Total_people, col="green3")
    
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
         ylim=c(0,1),
         xlim=c(year_plot, max(years)))
    text(2027,0.4,round(R$sym_r_over_sym[26],2),pos=4)
    add_k13_points()
    
    # k13_points/2
    plot(years, R$asym_r_over_asym,
         type="l", col="red",
         xlab="Year", ylab="Ratio",
         main="Asymptomatic Resistant Incidence / Asymptomatic Incidence",
         ylim=c(0,1),
         xlim=c(year_plot, max(years)))
    add_k13_points_2()
    text(2027,0.4,round(R$asym_r_over_asym[26],2),pos=4)
  }
  
  invisible(prep)
}

# co
