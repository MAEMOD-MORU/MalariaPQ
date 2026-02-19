# Model Baseline (D0 33%, D1 67%) file 
##### 2 #####
Malaria_model_with_Array<- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    # Define variables
    Is <- c(Is0, Is1, Is2) # Recreate Is array from individual components
    Ir <- c(Ir0, Ir1, Ir2) # Recreate Ir array from individual components
    GIs <- c(GIs0, GIs1, GIs2)  # Gametocytes from Is
    GIr <- c(GIr0, GIr1, GIr2)  # Gametocytes from Ir
    Stis <- c(Stis1,Stis2)
    Stir <- c(Stir1,Stir2)
    Fis <- c(Fis1,Fis2)
    Fir <- c(Fir1,Fir2)
    
    tau_is <- τigs
    tau_nfs <- τnfs
    tau_fs <- τfs
    tau_as <- Tas
    tau_ags <- Tags
    tau_ntsd0 <- τntsd0
    
    tau_ir <- τigr
    tau_nfr <- τnfr
    tau_fr <- τfr
    tau_ar <- Tar
    tau_agr <- Tagr
    tau_ntrd0 <-  τntrd0
    
    gamma_is <- g_is 
    gamma_ir <- g_ir 
    gamma_ar <- g_ar 
    gamma_as <- g_as 
    gamma_infs <- g_infs 
    gamma_ifs <- g_ifs 
    gamma_infr <- g_infr 
    gamma_ifr <- g_ifr  
    
    # For uganda
    # if(t >= start_m){
    #   mui <- mui_before
    #   muo <- muo_before
    # }
    # if(t >= start_b){
    #   beta_s <- beta_s_2
    #   beta_as <- beta_as_2
    # }
    
    # Rate of gametocyte to recovery
    if(t <=start_d){
      
      gamma_is <- g_is
      # D1 66%
      gamma_infs <- c(g_infs[1],g_infs[1])
      gamma_ifs <- c(g_ifs[1],g_ifs[1])
      
      gamma_ir <- g_ir
      # D1 66%
      gamma_infr <- c(g_infr[1],g_infr[1])
      gamma_ifr <- c(g_ifr[1],g_ifr[1])
      
      
      
      tau_is <- c(d0, d1, d1)
      tau_nfs <- τnfs
      tau_fs <- τfs
      
      tau_ir <- c(d0, d1, d1)
      tau_nfr <- τnfr
      tau_fr <- τfr
      
    }else if(t <=(start_d+t_long)){
      
      # print(T)
      gamma_is <- g_is 
      gamma_ir <- g_ir 
      gamma_infs <- g_infs 
      gamma_ifs <- g_ifs 
      gamma_infr <- g_infr 
      gamma_ifr <- g_ifr  
      
      tau_is <- c(d0, d1, d2)
      tau_nfs[2] <- d2
      tau_fs <- τfs
      
      tau_ir <- c(d0, d1, d2)
      tau_nfr[2] <- d2
      tau_fr <- τfr
    }else{
      
      gamma_is <- g_is
      gamma_infs <- c(g_infs[1],g_infs[1])
      gamma_ifs <- c(g_ifs[1],g_ifs[1])
      
      gamma_ir <- g_ir
      gamma_infr <- c(g_infr[1],g_infr[1])
      gamma_ifr <- c(g_ifr[1],g_ifr[1])
      
      
      
      tau_is <- c(d0, d1, d1)
      tau_nfs <- τnfs
      tau_fs <- τfs
      
      tau_ir <- c(d0, d1, d1)
      tau_nfr <- τnfr
      tau_fr <- τfr
      
    }
    
    if (t >= start_d) { 
      propIs <- propIs_new
      propIr <- propIr_new
    }
    # print(tau_ir)
    
    N <- S + sum(Is) + sum(GIs) + sum(Stis) + sum(Fis) + GAs + As + Rs + sum(Ir) + sum(GIr) + sum(Stir) + sum(Fir) + GAr + Ar + Rr
    season <- flo + amp * sin(2 * pi / period * (phase + t))
    
    
    lam_is <- (season * (beta_s * GIs) / N)
    lam_as <- (season * (beta_as * GAs) / N)
    
    # Applying the resistant fitness cost to only one 
    # symptomatic class while leaving other symptomatic
    # classes unaffected is not biologically justified,
    # because the fitness cost reflects an intrinsic property
    # of the parasite that influences mosquito infectivity
    # regardless of the host’s clinical state or treatment
    # pathway. If resistant parasites are less transmissible,
    # this disadvantage should apply consistently across
    # all sources of transmission (symptomatic and asymptomatic),
    # rather than selectively to a single compartment.
    lam_ir <- (season * (c_beta_r * beta_r * GIr) / N)
    lam_ar <- (season * (c_beta_r * beta_ar * GAr) / N)
    
    
    # diagnostics to quantify treated-symptomatic share ----------
    tiny <- 1e-12
    
    # Treated-symptomatic (T) = D1 + D2 gametocytes
    # Sensitive treated channel contribution 
    C_T_s <- season * (sum(beta_s[2:3] * GIs[2:3])) / N
    
    # Resistant treated channel contribution
    C_T_r <- season * (sum(beta_r[2:3] * GIr[2:3])) / N
    
    # Untreated/asymptomatic (U) = D0 + asym carriers
    # Sensitive untreated/asym
    C_U_s <- season * (beta_as * GAs + beta_s[1] * GIs[1]) / N
    
    # Resistant untreated/asym
    C_U_r <- season * (beta_ar * GAr + beta_r[1] * GIr[1]) / N
    
    # Channel totals (same scale as lam_* terms)
    CT_total <- C_T_s + C_T_r
    CU_total <- C_U_s + C_U_r
    
    # Fraction of transmission coming from treated symptomatic channel
    F_T <- CT_total / (CT_total + CU_total + tiny)
    
    # Resistant fraction of incident infections (instantaneous)
    p_resistant_inc <- sum(lam_ir + lam_ar) / sum(lam_is + lam_as + lam_ir + lam_ar + tiny)
    
    # Rate of change for each state
    dS <- mui * N - muo * S - sum(lam_is * S) - (lam_as * S) - sum(lam_ir * S)  - (lam_ar * S)  + alpha * Rs + alpha * Rr
    # print(sum(lam_is * S,(lam_as * S)) * (prob_sym_s) * propIs)
    dIs <- -muo * Is + sum(lam_is * S,(lam_as * S)) * (prob_sym_s) * propIs - Is * (gamma_is)
    dIs[1] <- dIs[1] - tau_ntsd0 * Is[1]
    dIs[2:3] <- dIs[2:3] - (1)* Is[2:3] # move to STis and Fis (Succeed/Fail to Treatment)
    # print(beta_s*Is)
    dStis <- -muo * Stis + (1-Fail_rate_s)* Is[2:3] - Stis * gamma_infs - Stis * tau_nfs
    dFis  <- -muo * Fis  + Fail_rate_s  * Is[2:3] - Fis * gamma_ifs - Fis * tau_fs
    # print(sum(tau_is * GIs ))
    dGIs <- -muo * GIs + gamma_is * Is  - tau_is * GIs # Add gametocyte rate for Is
    dGIs[2:3] <- dGIs[2:3]+ Stis * gamma_infs + Fis * gamma_ifs
    
    dAs <- -muo * As + sum(lam_is * S,(lam_as * S)) * (1-prob_sym_s) - tau_as * As - gamma_as * As
    dGAs <- -muo * GAs + gamma_as * As - tau_ags * GAs # Add gametocyte rate for As
    
    dRs <- -muo * Rs + tau_as * As - alpha * Rs + sum(tau_is * GIs) + (tau_ags * GAs) + tau_ntsd0 * Is[1] + sum(Stis * tau_nfs) + sum(Fis * tau_fs)
    
    dIr <- -muo * Ir + sum(lam_ir * S,(lam_ar * S)) * (prob_sym_r) * propIr - Ir * (gamma_ir)
    dIr[1] <- dIr[1]- tau_ntrd0 * Ir[1]
    dIr[2:3] <- dIr[2:3]- (1) * Ir[2:3] # move to STir and Fir (Succeed/Fail to Treatment)
    
    dStir <- -muo * Stir + (1-Fail_rate_r) * Ir[2:3] - Stir * gamma_infr - Stir * tau_nfr
    dFir  <- -muo * Fir  + Fail_rate_r  * Ir[2:3] - Fir * gamma_ifr - Fir * tau_fr
    
    dGIr <- -muo * GIr + gamma_ir * Ir - tau_ir * GIr  # Add gametocyte rate for Ir
    dGIr[2:3] <- dGIr[2:3]+Stir * gamma_infr + Fir * gamma_ifr 
    
    dAr <- -muo * Ar + sum(lam_ir * S,(lam_ar * S)) * (1-prob_sym_r) - tau_ar * Ar - gamma_ar * Ar
    dGAr <- -muo * GAr + gamma_ar * Ar - tau_agr * GAr # Add gametocyte rate for Ar
    
    dRr <- -muo * Rr + tau_ar * Ar - alpha * Rr + sum(tau_ir * GIr) + (tau_agr * GAr) + tau_ntrd0 * Ir[1] + sum(Stir * tau_nfr) + sum(Fir * tau_fr)
    
    # Return the rate of change
    list(c(
      dS,
      
      dIs, # Unpacked back into individual 
      dStis,
      dFis,
      dGIs,
      dGAs,
      dAs,
      dRs,
      
      dIr, # Unpacked back into individual components
      dStir,
      dFir,
      dGIr,
      dGAr,
      dAr,
      dRr
    )
    ,
    inc = sum(lam_is * S) + sum(lam_as * S) + sum(lam_ir * S)  + sum(lam_ar * S) ,
    inc_s = sum(lam_is * S) + sum(lam_as * S) ,
    inc_r = sum(lam_ir * S)  + sum(lam_ar * S) ,
    inc_sym = (sum(lam_is * S) +lam_as * S)* (prob_sym_s) + (sum(lam_ir * S) +lam_ar * S)* (prob_sym_r),
    inc_sym_s = (sum(lam_is * S) +lam_as * S)* (prob_sym_s),
    inc_sym_r = (sum(lam_ir * S) +lam_ar * S)* (prob_sym_r),
    inc_asym = (sum(lam_is * S) +lam_as * S)* (1-prob_sym_s) + (sum(lam_ir * S) +lam_ar * S)* (1-prob_sym_r),
    inc_asym_s = (sum(lam_is * S) +lam_as * S)* (1-prob_sym_s),
    inc_asym_r = (sum(lam_ir * S) +lam_ar * S)* (1-prob_sym_r),
    
    # Gametocyte Infections
    G_inc = sum(gamma_is * Is) + gamma_as * As + sum(gamma_ir * Ir) + gamma_ar * Ar,
    GS_inc = sum(gamma_is * Is) + gamma_as * As,
    GR_inc = sum(gamma_ir * Ir) + gamma_ar * Ar,
    Gsym_inc = sum(gamma_is * Is) + sum(gamma_ir * Ir),
    Gsym_inc_s = sum(gamma_is * Is),
    Gsym_inc_r = sum(gamma_ir * Ir),
    Gasym_inc = gamma_as * As + gamma_ar * Ar,
    Gasym_inc_s = gamma_as * As,
    Gasym_inc_r = gamma_ar * Ar,
    N = N,
    # diagnostics 
    F_T = F_T,                         # fraction of transmission from treated symptomatic (D1+D2)
    CT_total = CT_total,               # treated channel contribution (for plotting)
    CU_total = CU_total,               # untreated/asym channel contribution (for plotting)
    p_resistant_inc = p_resistant_inc  # resistant share of incident infections
    )
  })
}