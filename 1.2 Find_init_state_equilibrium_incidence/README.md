# 1.2 Find Initial State — Equilibrium Incidence

Finds the **equilibrium initial state** of the malaria transmission model by fitting transmission parameters (`beta_s`, `beta_s_2`) to observed incidence data under a **fixed (constant) population size**. The resulting equilibrium state vector is saved as the initial condition for Step 2 (incidence fitting).

> **Why fixed population?**
> Holding population size constant removes demographic drift from the system, making it straightforward to solve for the true disease equilibrium. The demographic parameters (`mui`, `muo`) fitted in Step 1.1 are used but population is held static during this calibration.

---

## Files

| File | Description |
|---|---|
| `init_parameter_calibration_infection.R` | Defines parameters and initial state inherited from Step 1.1, with fixed population |
| `1.calibration_infection.R` | Optimises transmission parameters to match observed incidence at equilibrium |
| `2.plot_equilibrium.R` | Loads saved results and plots the equilibrium incidence fit |
| `init_state_2000.rds` | Saved equilibrium compartment state vector (initial condition for Step 2) |
| `optimized_params_equilibrium.rds` | Saved optimisation output (parameter estimates + Hessian) |

---

## Workflow

```
init_parameter_calibration_infection.R  # 1. Load params, fix population size
        │
        ▼
1.calibration_infection.R               # 2. Fit beta_s, beta_s_2 → equilibrium state
        │                                      saves: init_state_2000.rds
        │                                             optimized_params_equilibrium.rds
        ▼
2.plot_equilibrium.R                    # 3. Plot equilibrium incidence fit
```

---

## Scripts

### `init_parameter_calibration_infection.R`
- Loads fitted demographic parameters (`mui`, `muo`) from `../1.1 Fitting_population/calibration_population_results.rds`
- Fixes total population to the year-2000 value (`initP = 34,260,139`) — population does **not** change during this optimisation
- Sets initial compartment state and all other fixed model parameters (transmission rates, recovery rates, gametocyte rates, treatment failure rates)

### `1.calibration_infection.R`
Fits transmission rate parameters by minimising the difference between modelled and observed annual incidence:

**Parameters optimised:**

| Parameter | Description |
|---|---|
| `beta_s` | Transmission rate — Before 2020 |
| `beta_s_2` | Transmission rate — After 2020 |

The resistant strain transmission rates are derived from these via the relative fitness parameter `c_beta_r`.

Runs the ODE model to equilibrium and saves:
```r
saveRDS(init_state, file = "init_state_2000.rds")
saveRDS(op,         file = "optimized_params_equilibrium.rds")
```

### `2.plot_equilibrium.R`
- Loads `optimized_params_equilibrium.rds`
- Re-runs the model at fitted parameters
- Plots modelled equilibrium incidence against observed incidence data

---

## Output

| File | Used by | Description |
|---|---|---|
| `init_state_2000.rds` | Step 2 | Equilibrium compartment state — initial condition for incidence fitting |
| `optimized_params_equilibrium.rds` | Step 2 | Fitted transmission parameters at equilibrium |

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
read.csv("../data/Tanzania_inc1000_pop.csv")   # Observed incidence per 1,000
read.csv("../data/Tanzania_pop_2000_2100.csv") # Population (for fixed initial size)
```
