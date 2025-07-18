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
    
    if(t >= start_m){
      mui <- mui_2025
      muo <- muo_2025
    }
    
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
      tau_nfs <- c(d1, d1)
      tau_fs <- c(d1, d1)

      tau_ir <- c(d0, d1, d1)
      tau_nfr <- c(d1, d1)
      tau_fr <- c(d1, d1)

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
      gamma_infs <- c(g_infs[1],g_infs[1])
      gamma_ifs <- c(g_ifs[1],g_ifs[1])
      
      gamma_ir <- g_ir
      gamma_infr <- c(g_infr[1],g_infr[1])
      gamma_ifr <- c(g_ifr[1],g_ifr[1])
      
      
      
      tau_is <- c(d0, d1, d1)
      tau_nfs <- c(d1, d1)
      tau_fs <- c(d1, d1)
      tau_ir <- c(d0, d1, d1)
      tau_nfr <- c(d1, d1)
      tau_fr <- c(d1, d1)
      
    }
    
    N <- S + sum(Is) + sum(GIs) + sum(Stis) + sum(Fis) + GAs + As + Rs + sum(Ir) + sum(GIr) + sum(Stir) + sum(Fir) + GAr + Ar + Rr
    season <- flo + amp * sin(2 * pi / period * (phase + t))
    
    # Compute lambda
    lam_is <- sum(season * (beta_s * GIs) / N)
    lam_as <- sum(season * (beta_as * GAs) / N)
    lam_ir <- sum(season * (beta_r * GIr) / N)
    lam_ar <- sum(season * (beta_ar * GAr) / N)
    
    # Rate of change for each state
    dS <- mui * N - muo * S - sum(lam_is * S) - sum(lam_as * S) - sum(lam_ir * S)  - sum(lam_ar * S)  + alpha * Rs + alpha * Rr
    
    dIs <- -muo * Is + (lam_is * S) * propIs - Is * (gamma_is)
    dIs[1] <- dIs[1] - tau_ntsd0 * Is[1]
    dIs[2:3] <- dIs[2:3] - (1) * Is[2:3] # move to STis and Fis (Succeed/Fail to Treatment)
    
    dStis <- -muo * Stis + (1-Fail_rate_s)* Is[2:3] - Stis * gamma_infs - Stis * tau_nfs
    dFis  <- -muo * Fis  + Fail_rate_s  * Is[2:3] - Fis * gamma_ifs - Fis * tau_fs
    
    dGIs <- -muo * GIs + gamma_is * Is  - tau_is * GIs # Add gametocyte rate for Is
    dGIs[2:3] <- dGIs[2:3]+ Stis * gamma_infs + Fis * gamma_ifs

    dAs <- -muo * As + (lam_as * S) - tau_as * As - gamma_as * As
    dGAs <- -muo * GAs + gamma_as * As - tau_ags * GAs # Add gametocyte rate for As
    
    dRs <- -muo * Rs + tau_as * As - alpha * Rs + sum(tau_is * GIs) + (tau_ags * GAs) + tau_ntsd0 * Is[1] + sum(Stis * tau_nfs) + sum(Fis * tau_fs)
    
    dIr <- -muo * Ir + (lam_ir * S) * propIr - Ir * (gamma_ir)
    dIr[1] <- dIr[1]- tau_ntrd0 * Ir[1]
    dIr[2:3] <- dIr[2:3]- (1) * Ir[2:3] # move to STir and Fir (Succeed/Fail to Treatment)
    
    dStir <- -muo * Stir + (1-Fail_rate_r) * Ir[2:3] - Stir * gamma_infr - Stir * tau_nfr
    dFir  <- -muo * Fir  + Fail_rate_r  * Ir[2:3] - Fir * gamma_ifr - Fir * tau_fr
    
    dGIr <- -muo * GIr + gamma_ir * Ir - tau_ir * GIr  # Add gametocyte rate for Ir
    dGIr[2:3] <- dGIr[2:3]+Stir * gamma_infr + Fir * gamma_ifr 
    
    dAr <- -muo * Ar + (lam_ar * S) - tau_ar * Ar - gamma_ar * Ar
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
    inc_sym = sum(lam_is * S) + sum(lam_ir * S),
    inc_sym_s = sum(lam_is * S),
    inc_sym_r = sum(lam_ir * S),
    inc_asym = sum(lam_as * S) + sum(lam_ar * S),
    inc_asym_s = sum(lam_as * S),
    inc_asym_r = sum(lam_ar * S),
    # Gametocyte Infections
    GS_inc = sum(lam_is * GIs) + sum(lam_as * GAs),
    GR_inc = sum(lam_ir * GIr) + sum(lam_ar * GAr),
    Gsym_inc = sum(lam_is * GIs) + sum(lam_ir * GIr),
    Gasym_inc = sum(lam_as * GAs) + sum(lam_ar * GAr),
    N = N
    
    
    )
  })
}



