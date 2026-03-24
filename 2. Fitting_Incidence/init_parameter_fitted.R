# Initial Population
initP <- 34260139

total_inc <- 1000

prob_sym_res <-0

mui <- 0.0226    # form 1.1 Fitting_population
muo <- 0.0203    # form 1.1 Fitting_population
mui_before <- 0.0226    # form 1.1 Fitting_population
muo_before <- 0.0203    # form 1.1 Fitting_population

# Initial States from UI Inputs
init_total_IAs <- total_inc*(1-prob_sym_res)
init_total_IAr <- total_inc*(prob_sym_res)

prob_sym_s = 0.25 # Ratio of infection on sym/asym sensitive group
prob_sym_r = 0.25 # Ratio of infection on sym/asym resistance group

Is_total <- init_total_IAs * prob_sym_s
As_total <- init_total_IAs * (1-prob_sym_s)

Ir_total <- init_total_IAr * prob_sym_r
Ar_total <- init_total_IAr * (1-prob_sym_r)

init_D1s <- 0.33
init_D2s <- 0

init_D1r <- 0.33
init_D2r <- 0

propIs <- c(1 - (init_D1s + init_D2s), init_D1s, init_D2s) # ratio of
propIr <- c(1 - (init_D1r + init_D2r), init_D1r, init_D2r)

initIs <- Is_total * propIs  # Sensitive subgroups
initAs <- As_total
initGIs <- c(0, 0, 0)
initGAs <- 0
initRs <- 0

initIr <- Ir_total * propIr  # Resistant subgroups
initAr <- Ar_total
initGIr <- c(0, 0, 0)
initGAr <- 0
initRr <- 0

initStis <- c(0,0) # Succeed to Treatment (Sensitive) (D1,D2)
initFis <- c(0,0) # Fail to Treatment (Sensitive) (D1,D2)

initStir <- c(0,0) # Succeed to Treatment (Resistant) (D1,D2)
initFir <- c(0,0) # Fail to Treatment (Resistant) (D1,D2)

# Calculate initial susceptible population
initS <- initP - sum(initIs) - sum(initGIs) - initGAs - initAs - initRs -
  sum(initIr) - sum(initGIr) - initGAr - initAr - initRr

state <- c(
  S = initS,
  Is0 = initIs[1], Is1 = initIs[2], Is2 = initIs[3], # Flatten Is array into individual components
  Stis1 = initStis[1],Stis2=initStis[2],
  Fis1 = initFis[1],Fis2=initFis[2],
  GIs0 = initGIs[1], GIs1 = initGIs[2], GIs2 = initGIs[3],  # Gametocyte from Is
  GAs = initGAs,  # Gametocyte from As
  As = initAs,
  Rs = initRs,
  Ir0 = initIr[1], Ir1 = initIr[2], Ir2 = initIr[3], # Flatten Ir array into individual components
  Stir1 = initStir[1],Stir2 = initStir[2],
  Fir1 = initFir[1],Fir2 = initFir[2],
  GIr0 = initGIr[1], GIr1 = initGIr[2], GIr2 = initGIr[3],  # Gametocyte from Ir
  GAr = initGAr,  # Gametocyte from Ar
  Ar = initAr,
  Rr = initRr
)


# Dynamic tau from UI
τnfs <- rep(4.3,2) # Succeed to Treatment group to Recovered (Sensitive)
τigs <- c(0.33,0.68,6.15) # Symptomatic Gametocytes to Recovered (Sensitive)
τfs <-rep(0.33,2) # Fail to Treatment group to Recovered (Sensitive)
τntsd0 <- 0.33 # D0 to Recovered (Sensitive)
Tas <- 0.33 # Asymptomatic to Recovered (Sensitive)
Tags <- 0.33 # Asymptomatic Gametocytes to Recovered (Sensitive)

τnfr <- rep(4.3,2)  # Succeed to Treatment group to Recovered (Resistant)
τigr <- c(0.33,0.68,6.15) # Symptomatic Gametocytes to Recovered (Resistant)
τfr <- rep(0.33,2) # Fail to Treatment group to Recovered (Resistant)
τntrd0 <- 0.33 # D0 to Recovered (Resistant)
Tar <- 0.33 # Asymptomatic to Recovered (Resistant)
Tagr <- 0.33 # Asymptomatic Gametocytes to Recovered (Resistant)

g_is = c(0.121,0.121,0.121) #0.1
g_ir = c(0.242,0.242,0.242)

g_as = 0.186
g_ar = 0.372

g_infs = c(0.019,0.019)
g_ifs = c(1.21,1.21)  #1

g_infr = c(0.038,0.038)
g_ifr = c(2.42,2.42) 

beta_is = 1.436          
beta_as = 1.436       
beta_is_2 = 1.436        
beta_as_2 = 1.436       
beta_r = 1.436 *1.05
beta_ar = 1.436 *1.05 
beta_ir_2 = 1.436 *1.05   
beta_ar_2 = 1.436 *1.05

# mui <- 0.01
# muo <- 0.01
c_beta_r <- 0.7622   # resistant parasites transmit 20% less outside treatment


parameters <- list(
  mui = mui,
  muo = muo,
  mui_before = mui_before,
  muo_before = muo_before,
  flo = 0.9,
  amp = 0.1,
  propIs=propIs,
  propIr=propIr,
  propIs_new = c(0.67, 0 ,0.33),
  propIr_new = c(0.67, 0 ,0.33),
  period = 12,
  phase = 0,
  prob_sym_s = prob_sym_s,
  prob_sym_r = prob_sym_r,
  start_d = 6000,# 10 years after 3500 months
  start_m = 12*1000, # 2015
  start_b = 12*1000, # 2015
  t_long = 12*100, # 2,5,10 years 
  Fail_rate_s = 0.025, # Fail Treatment of ACT
  Fail_rate_r = 0.025,
  
  beta_s = c(beta_is, beta_is, beta_is),
  beta_r = c(beta_r, beta_r,beta_r),
  beta_as = beta_as,
  beta_ar = beta_ar,
  beta_s_2 =  c(beta_is_2, beta_is_2, beta_is_2),
  beta_as_2 = beta_as_2,
  beta_ir_2 = c(beta_ir_2,beta_ir_2,beta_ir_2),
  beta_ar_2 = beta_ar_2,
  
  alpha = 1.07 / 12, #  doi:10.1371/journal.pone.0001767
  τigs=τigs,
  Tas=Tas,
  Tags=Tags,
  
  τigr=τigr,
  Tar=Tar,
  Tagr=Tagr,
  
  τntsd0= τntsd0,
  τntrd0= τntrd0,
  
  τnfs = τnfs ,
  τfs = τfs ,
  
  τnfr=τnfr,
  τfr=τfr,
  
  d0 = (1 / 90) * 30, # no drugs
  d1 = (1 / 44.1) * 30, # Gametocyte recovery rate from Palang paper, table 2 at row Placebo (DP), column microscopy
  d2 = (1 / 1) * 30, # Gametocyte recovery rate from Palang paper, table 2 at row 0.0625 mg/kg, column microscopy
  
  g_is = g_is, #0.1
  g_ir = g_ir,
  g_ar = g_ar,
  g_as = g_as,
  g_infs = g_infs,
  g_ifs = g_ifs,  #1
  g_infr = g_infr,
  g_ifr = g_ifr,

  c_beta_r = c_beta_r
)

times <- seq(0, 6000, by = 1)
