# 3. Results

Contains plotting scripts and output figures for the **MalariaPQ** project. Each subfolder corresponds to one publication figure and includes the R script used to generate it along with any output image files.

> **Requires output from Step 2:**
> - `../2. Fitting_Incidence/MCMC_out.rds` — posterior samples from MCMC fitting

---

## Files

| File | Description |
|---|---|
| `init_parameter_MCMC.R` | extracts fitted parameter estimates from MCMC |
| `init_state_2000.rds` | Equilibrium initial state from Step 1.2 |
| `Figure 1 CI95/` |  plots with 95% credible intervals |
| `Figure 2/` | Scenario comparison D2 treatment scenarios |
| `Figure 3/` | Scenario comparison start Ratio Resistance and D2 |

---

## Workflow

```
../2. Fitting_Incidence/MCMC_out.rds
        │
        ▼
init_parameter_MCMC.R       # extracts fitted parameters
and init_state_2000.rds
        │
        ├──► Figure 1 CI95/  # Incidence + K13 resistance ratio with 95% CI
        ├──► Figure 2/        # D2 scenario comparisons
        └──► Figure 3/        # start Ratio Resistance scenario comparisons
```

---

## Figure Descriptions

### `Figure 1 CI95/` — Posterior Predictive with 95% Credible Intervals
Draws posterior predictive samples from `MCMC_out.rds` and plots:

### `Figure 2/` — Scenario Comparison (D2 Scenarios)
Compares model projections under different treatment coverage or D2 scenario assumptions.

### `Figure 3/` — Additional Figures (Start Ratio Resistance and D2)
Compares model projections across different combinations of initial resistance ratio and D2 treatment scenarios.

---
