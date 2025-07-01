# Model Baseline (D0 33%, D1 67%) file 
##### 2 #####
Malaria_model_with_Array_2 <- function(t, state, parameters) {
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
    
    # Rate of gametocyte to recovery
    if(t <=start_d){

      gamma_is <- g_is 
      # gamma_as <- g_as 
      gamma_infs <- g_is[2:3]
      gamma_ifs <- g_is [2:3]
      
      gamma_ir <- g_ir 
      # gamma_ar <- g_ar 
      gamma_infr <- g_ir[2:3] 
      gamma_ifr <- g_ir[2:3]  
      
      tau_is <- c(d0, d1, d1)
      tau_nfs <- c(d1, d1)
      tau_fs <- c(d1, d1)
      # tau_as <- c(d0, d0, d0)
      # tau_ags <- c(d0, d0, d0)
      # tau_ntsd0 <- c(d0, d0, d0)
      tau_ir <- c(d0, d1, d1)
      tau_nfr <- c(d1, d1)
      tau_fr <- c(d1, d1)
      # tau_ar <- c(d0, d0, d0)
      # tau_agr <- c(d0, d0, d0)
      # tau_ntrd0 <- c(d0, d0, d0)
    }else if(t <=(start_d+t_long)){
      # print(T)
      gamma_is <- g_is 
      gamma_ir <- g_ir 
      gamma_infs <- g_infs 
      gamma_ifs <- g_ifs 
      gamma_infr <- g_infr 
      gamma_ifr <- g_ifr  
      
      tau_is <- c(d0, d1, d2)
      tau_nfs <- τnfs
      tau_fs <- τfs
      
      tau_ir <- c(d0, d1, d2)
      tau_nfr <- τnfr
      tau_fr <- τfr
    }else{
      
      gamma_is <- g_is 
      # gamma_as <- g_as 
      gamma_infs <- g_is[2:3]
      gamma_ifs <- g_is [2:3]
      
      gamma_ir <- g_ir 
      # gamma_ar <- g_ar 
      gamma_infr <- g_ir[2:3] 
      gamma_ifr <- g_ir[2:3] 
      
      tau_is <- c(d0, d1, d1)
      tau_nfs <- c(d1, d1)
      tau_fs <- c(d1, d1)

      tau_ir <- c(d0, d1, d1)
      tau_nfr <- c(d1, d1)
      tau_fr <- c(d1, d1)

    }
    
    # beta_s <- c(beta_s0, beta_s1, beta_s2) # Define beta for Is subgroups
    # beta_r <- c(beta_r0, beta_r1, beta_r2) # Define beta for Ir subgroups
    # beta_as <- beta_s0
    # beta_ar <- beta_r0
    
    N <- S + sum(Is) + sum(GIs) + sum(Stis) + sum(Fis) + GAs + As + Rs + sum(Ir) + sum(GIr) + sum(Stir) + sum(Fir) + GAr + Ar + Rr
    season <- flo + amp * sin(2 * pi / period * (phase + t))
    
    # Compute lambda for Symtomatic and Asymtomatic group
    lam <- season * (sum(beta_s * GIs) + beta_as * GAs) / N +
      season * (sum(beta_r * GIr) + beta_ar * GAr) / N
    # print(paste0("lam_s: ",lam_s))
    # print(paste0("lam_r: ",lam_r))
    
    # Rate of change for each state
    dS <- mui * N - muo * S - (prob_res*lam * S) - ((1-prob_res)*lam * S) + alpha * Rs + alpha * Rr
    
    # dIs <- -muo * Is + lam_s * S * prob_s * propIs - gamma_i * Is - tau_ntsd0 * Is[1] - (1-gamma_i-muo) * Is[2:3] # Array operation for Is
    dIs <- -muo * Is + (prob_res*lam * S) * prob_s * propIs - Is * (gamma_is)
    dIs[1] <- dIs[1] - tau_ntsd0 * Is[1]
    dIs[2:3] <- dIs[2:3] - (1) * Is[2:3] # move to STis and Fis (Succeed/Fail to Treatment)
    
    dStis <- -muo * Stis + (1-Fail_rate_s)* Is[2:3] - Stis * gamma_infs - Stis * tau_nfs
    dFis  <- -muo * Fis  + Fail_rate_s  * Is[2:3] - Fis * gamma_ifs - Fis * tau_fs
    
    dGIs <- -muo * GIs + gamma_is * Is  - tau_is * GIs # Add gametocyte rate for Is
    dGIs[2:3] <- dGIs[2:3]+ Stis * gamma_infs + Fis * gamma_ifs
    
    dAs <- -muo * As + (prob_res*lam * S) * (1 - prob_s) - tau_as * As - gamma_as * As
    dGAs <- -muo * GAs + gamma_as * As - tau_ags * GAs # Add gametocyte rate for As
    
    dRs <- -muo * Rs + tau_as * As - alpha * Rs + sum(tau_is * GIs) + (tau_ags * GAs) + tau_ntsd0 * Is[1] + sum(Stis * tau_nfs) + sum(Fis * tau_fs)
    
    # dIr <- -muo * Ir + lam_r * S * prob_r * propIr - gamma_i * Ir - tau_ntrd0 * Ir[1] - (1-gamma_i) * Ir[2:3] # Array operation for Is
    dIr <- -muo * Ir + ((1-prob_res)*lam * S) * prob_r * propIr - Ir * (gamma_ir)
    dIr[1] <- dIr[1]- tau_ntrd0 * Ir[1]
    dIr[2:3] <- dIr[2:3]- (1) * Ir[2:3] # move to STir and Fir (Succeed/Fail to Treatment)
    
    # dStir <- -muo * Stir + (1-Fail_rate)*(1-gamma_i) * Ir[2:3] - Stir * gamma_infr - Stir * tau_nfr
    # dFir  <- -muo * Fir  + Fail_rate*(1-gamma_i)  * Ir[2:3] - Fir * gamma_ifr - Fir * tau_fr
    dStir <- -muo * Stir + (1-Fail_rate_r) * Ir[2:3] - Stir * gamma_infr - Stir * tau_nfr
    dFir  <- -muo * Fir  + Fail_rate_r  * Ir[2:3] - Fir * gamma_ifr - Fir * tau_fr
    
    
    
    dGIr <- -muo * GIr + gamma_ir * Ir - tau_ir * GIr  # Add gametocyte rate for Ir
    dGIr[2:3] <- dGIr[2:3]+Stir * gamma_infr + Fir * gamma_ifr 
    
    dAr <- -muo * Ar + ((1-prob_res)*lam * S) * (1 - prob_r) - tau_ar * Ar - gamma_ar * Ar
    dGAr <- -muo * GAr + gamma_ar * Ar - tau_agr * GAr # Add gametocyte rate for Ar
    
    dRr <- -muo * Rr + tau_ar * Ar - alpha * Rr + sum(tau_ir * GIr) + (tau_agr * GAr) + tau_ntrd0 * Ir[1] + sum(Stir * tau_nfs) + sum(Fir * tau_fs)
    
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
    ),
    inc = (prob_res*lam * S) + ((1-prob_res)*lam * S),
    inc_s = (prob_res*lam * S) ,
    inc_r = ((1-prob_res)*lam * S) ,
    inc_sym = ((prob_res*lam * S) * (prob_s)) + (((1-prob_res)*lam * S) * (prob_s)),
    inc_sym_s = ((prob_res*lam * S) * (prob_s)),
    inc_sym_r =(((1-prob_res)*lam * S) * (prob_s)),
    inc_asym = ((prob_res*lam * S) * (1 - prob_s)) + (((1-prob_res)*lam * S) * (1 - prob_r)),
    inc_asym_s =((prob_res*lam * S) * (1 - prob_s)),
    inc_asym_r =(((1-prob_res)*lam * S) * (1 - prob_r)),
    total_inf = sum(Is) + sum(Stis) + sum(Fis) + As + sum(Ir) + sum(Stir) + sum(Fir) + Ar,
    total_G = sum(GIs) + GAs +sum(GIr) + GAr,
    Gs = sum(GIs) + GAs,
    Gs_sym = sum(GIs),
    Gs_asym = GAs,
    Gr = sum(GIr) + GAr,
    Gr_sym = sum(GIr),
    Gr_asym = GAr,
    Gsym = sum(GIs)+sum(GIr),
    Gasym = GAs+GAr,
    N = N,
    malaria_case = sum(Is) + sum(GIs) + sum(Stis) + sum(Fis) + GAs + As + sum(Ir) + sum(GIr) + sum(Stir) + sum(Fir) + GAr + Ar,
    GS_inc = sum(gamma_is * Is) + sum(Stis * gamma_infs) + sum(Fis * gamma_ifs)+gamma_as * As,
    GR_inc = sum(gamma_ir * Ir) + sum(Stir * gamma_infr) + sum(Fir * gamma_ifr)+ gamma_ar * Ar,
    Gsym_inc = sum(gamma_is * Is) + sum(Stis * gamma_infs) + sum(Fis * gamma_ifs)+ sum(gamma_ir * Ir) + sum(Stir * gamma_infr) + sum(Fir * gamma_ifr),
    Gsym_inc_s = sum(gamma_is * Is) + sum(Stis * gamma_infs) + sum(Fis * gamma_ifs),
    Gsym_inc_r = sum(gamma_ir * Ir) + sum(Stir * gamma_infr) + sum(Fir * gamma_ifr),
    Gasym_inc = gamma_as * As+ gamma_ar * Ar,
    Gasym_inc_s = gamma_as * As,
    Gasym_inc_r = gamma_ar * Ar
    
    )
  })
}



