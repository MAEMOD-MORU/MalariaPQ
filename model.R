# Model file
Malaria_model_with_Array <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Define variables
    Is <- c(Is0, Is1, Is2) # Recreate Is array from individual components
    Ir <- c(Ir0, Ir1, Ir2) # Recreate Ir array from individual components
    GIs <- c(GIs0, GIs1, GIs2)  # Gametocytes from Is
    GIr <- c(GIr0, GIr1, GIr2)  # Gametocytes from Ir

    # Rate of gametocyte to recovery
    drug <- c(d0, d1, d2)

    # beta_s <- c(beta_s0, beta_s1, beta_s2) # Define beta for Is subgroups
    # beta_r <- c(beta_r0, beta_r1, beta_r2) # Define beta for Ir subgroups
    # beta_as <- beta_s0
    # beta_ar <- beta_r0

    N <- S + sum(Is) + sum(GIs) + GAs + As + Rs + sum(Ir) + sum(GIr) + GAr + Ar + Rr
    season <- flo + amp * sin(2 * pi / period * (phase + t))

    # Compute lambda for Symtomatic and Asymtomatic group
    lam_s <- season * (sum(beta_s * GIs) + beta_as * GAs) / N
    lam_r <- season * (sum(beta_r * GIr) + beta_ar * GAr) / N

    # Compute tau separately for each subgroup
    tau_is <- tau[c(2,3,4)]
    tau_as <- tau[1]
    tau_ir <- tau[c(6,7,8)]
    tau_ar <- tau[5]

    # Rate of change for each state
    dS <- mui * N - muo * S - lam_s * S - lam_r * S + alpha * Rs + alpha * Rr

    dIs <- -muo * Is + lam_s * S * prob_s * propIs - tau_is * Is - gamma_i * Is # Array operation for Is
    dGIs <- -muo * GIs + gamma_i * Is - drug * GIs # Add gametocyte rate for Is
    dAs <- -muo * As + lam_s * S * (1 - prob_s) - tau_as * As - gamma_a * As
    dGAs <- -muo * GAs + gamma_a * As - d0 * GAs # Add gametocyte rate for As
    dRs <- -muo * Rs + sum(tau_is * Is) + tau_as * As - alpha * Rs + sum(drug * GIs) + sum(d0 * GAs)

    dIr <- -muo * Ir + lam_r * S * prob_r * propIr - tau_ir * Ir - gamma_i * Ir # Array operation for Ir
    dGIr <- -muo * GIr + gamma_i * Ir - drug * GIr # Add gametocyte rate for Ir
    dAr <- -muo * Ar + lam_r * S * (1 - prob_r) - tau_ar * Ar - gamma_a * Ar
    dGAr <- -muo * GAr + gamma_a * Ar - d0 * GAr # Add gametocyte rate for Ar
    dRr <- -muo * Rr + sum(tau_ir * Ir) + tau_ar * Ar - alpha * Rr + sum(drug * GIr) + sum(d0 * GAr)

    # Return the rate of change
    list(c(
      dS,
      dIs, # Unpacked back into individual components
      dGIs,
      dGAs,
      dAs,
      dRs,
      dIr, # Unpacked back into individual components
      dGIr,
      dGAr,
      dAr,
      dRr
    ),
    inc = lam_s * S + lam_r * S + lam_s * Rs + lam_r * Rr,
    inc_s = lam_s * S + lam_s * Rs,
    inc_r = lam_r * S + lam_r * Rr,
    Gs = sum(GIs) + GAs,
    Gr = sum(GIr) + GAr,
    N = N
    )
  })
}
