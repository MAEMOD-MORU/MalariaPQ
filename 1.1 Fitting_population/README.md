# 1.1 Fitting Population

Calibrates the demographic parameters of the malaria transmission model to Tanzania population data (2000–2050) using numerical optimisation. This is the **first step** of the pipeline — its outputs feed into all downstream steps.

---

## Files

| File | Description |
|---|---|
| `init_pop_calibration.R` | Defines the initial compartment state vector and all fixed model parameters |
| `1.calibration_pop.R` | Fits birth rate (`mui`) and death rate (`muo`) to population data via L-BFGS-B optimisation |
| `2.plot_Fitted_result_CI95.R` | Loads saved results and plots the fitted population with 95% confidence intervals |
| `calibration_population_results.rds` | Saved output from optimisation (parameter estimates + Hessian) |
| `Plot SSP/` | (Archived) SSP scenario population plots — not used in current analysis |

---

## Workflow

```
init_pop_calibration.R          # 1. Define initial state & parameters
        │
        ▼
1.calibration_pop.R             # 2. Optimise mui, muo → save .rds
        │
        ▼
2.plot_Fitted_result_CI95.R     # 3. Load .rds → plot fit + CI95
```

---

## Scripts

### `init_pop_calibration.R`
Sets up everything needed before running the model:

- **Initial population**: Tanzania year-2000 population (`initP = 34,260,139`)
- **Compartment state vector** (`state`): Susceptible (`S`), Symptomatic infected sensitive/resistant (`Is`, `Ir`) across treatment day strata (D0/D1/D2), Asymptomatic (`As`, `Ar`), Gametocyte carriers (`GIs`, `GIr`, `GAs`, `GAr`), Treatment success/failure groups (`Stis`, `Fis`, `Stir`, `Fir`), Recovered (`Rs`, `Rr`)

### `1.calibration_pop.R`
Fits monthly **birth rate** (`mui`) and **death rate** (`muo`) by minimising RMSD between modelled and observed total population:

```r
# Objective function
rmsd <- function(totalpop_model) { ... }

# Optimisation (L-BFGS-B)
op <- optim(c(0.035, 0.00825), pop_run, hessian = TRUE)
```

Saves the full optimisation result (estimates + Hessian) to:
```
calibration_population_results.rds
```

### `2.plot_Fitted_result_CI95.R`
Loads `calibration_population_results.rds` and constructs 95% CI from the Hessian:

```r
vcov  <- -solve(op$hessian)
se    <- sqrt(diag(vcov))
lower <- op$par - 1.96 * se
upper <- op$par + 1.96 * se
```

CI bounds are applied **asymmetrically** to construct population envelopes:

| Bound | Birth rate `mui` | Death rate `muo` |
|---|---|---|
| Upper population | Upper CI | Lower CI |
| Lower population | Lower CI | Upper CI |

Produces two plots:
1. **Fitted population** — model line vs. observed data points
2. **Population with 95% CI** — shaded uncertainty band (2000–2050)

---

## Output

| Object | Description |
|---|---|
| `op$par[1]` — `mui` | Fitted monthly birth rate |
| `op$par[2]` — `muo` | Fitted monthly death rate |
| `calibration_population_results.rds` | Full `optim()` object (used by downstream steps) |

---

## Dependencies

```r
library(deSolve)   # ODE solver
library(dplyr)     # Data manipulation
```

Model ODE sourced from:
```r
source("../Model/model_D0_33_D1_67_baseline_change_birth_death_2betaSets_SS.R")
```

Data sourced from:
```r
read.csv("../data/Tanzania_pop_2000_2100.csv")  # rows 1–51 (year 2000–2050)
```
