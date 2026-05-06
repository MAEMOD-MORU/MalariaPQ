# 2. Fitting Incidence

Fits the malaria transmission model to observed **incidence** and **K13 allele frequency** data simultaneously using **Bayesian MCMC** (Differential Evolution with snooker update, DEzs). This step estimates five parameters governing malaria transmission and drug resistance dynamics in Tanzania (2000–2025).

> **Requires outputs from Steps 1.1 and 1.2:**
> - `../1.1 Fitting_population/calibration_population_results.rds` — fitted demographic parameters
> - `init_state_2000.rds` — equilibrium initial state from Step 1.2

---

## Files

| File | Description |
|---|---|
| `init_parameter_optim_fitted.R` | Loads fitted parameters from Steps 1.1 & 1.2, sets up model parameters and events |
| `1.fitting_model_MCMC_pois_change_distribution.R` | Runs MCMC fitting — defines likelihood, priors, and sampler; saves posterior samples |
| `init_state_2000.rds` | Equilibrium initial state from Step 1.2 (read-only input) |
| `MCMC_out.rds` | Saved MCMC output — posterior samples from all chains |

---

## Workflow

```
                                             init_state_2000.rds and 
                              init_parameter_optim_fitted.R(birth rate from 1.1 and beta from 1.2)
                                                  │
                                                  ▼
                                         1.fitting_model_MCMC_pois_change_distribution.R
                                                  │  Define likelihood (Poisson + Beta)
                                                  │  Define priors
                                                  │  Run DEzs MCMC (3 chains)
                                                  │  Check convergence (Gelman-Rubin)
                                                  ▼
                                             MCMC_out.rds
```

---

## Dependencies

```r
library(deSolve)        # ODE solver
library(BayesianTools)  # MCMC framework (DEzs sampler)
library(scales)         # Plot utilities
```

Model sourced from:
```r
source("../model/model_D0_33_D1_67_baseline_change_birth_death_2betaSets_SS.R")
```

Data sourced from:
```r
read.csv("../data/Reported malaria cases by method of confirmation.csv")  # rows 6–14, col 4
read.csv("../data/Tanzania_K13_Allele_frequency.csv")                     # col 2
read.csv("../data/Tanzania_pop_2000_2100.csv")
```


## Parameters Estimated

| Parameter | Description | Prior |
|---|---|---|
| `beta_s` | Transmission rate — sensitive strain, season 1 | Uniform(0.01, 5) |
| `beta_s_2` | Transmission rate — sensitive strain, season 2 | Uniform(0.01, 5) |
| `beta_r_multiply` | Resistant strain transmission multiplier (≥ 1) | Uniform(1, 1.5) |
| `c_beta_r` | Relative competitive fitness of resistant strain | Normal(0.75, 0.1) |
| `φ` (phi) | Precision parameter for Beta likelihood (K13) | Uniform(1, 25) |

> Resistant strain transmission rates are derived as:
> `beta_r = beta_s × beta_r_multiply`

---

## Likelihood Function

The log-likelihood combines two components:

**1. Poisson — Annual Malaria Incidence (2015–2023)**
```r
ll_inc <- sum(dpois(x = obs_inc, lambda = inc_model[2015:2023], log = TRUE))
```

**2. Beta — K13 Allele Frequency (2016–2022, 2025)**
```r
shape1 <- model_ratio * phi
shape2 <- (1 - model_ratio) * phi
ll_k13 <- sum(dbeta(x = obs_k13, shape1 = shape1, shape2 = shape2, log = TRUE))
```

Total: `ll = ll_inc + ll_k13`

---

## MCMC Settings

| Setting | Value |
|---|---|
| Sampler | DEzs (Differential Evolution with snooker update) |
| Chains | 3 |
| Burn-in | 1,500 iterations |
| Sampling | 7,500 iterations |
| Total per chain | 9,000 iterations |

Just a exmaple. you can change number of iterations and Chains.

### Convergence Check
Convergence is assessed using the **Gelman-Rubin diagnostic** (multivariate PSRF):

- ✅ **PSRF < 1.1** — chains have converged; proceed to plotting
- ⚠️ **PSRF ≥ 1.1** — chains have **not** converged; run additional iterations:

```r
# Continue from existing chains with new settings
MCMC_out_2nd <- runMCMC(MCMC_out, settings = list(iterations = 15000, burnin = 5000))
saveRDS(MCMC_out_2nd, "MCMC_out_2nd.rds")
```
---

## Output

| File | Used by | Description |
|---|---|---|
| `MCMC_out.rds` | `3. Results/` | Full posterior samples from 3 MCMC chains |

---

