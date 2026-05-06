# Figure 1 CI95 — Posterior Predictive Plot with 95% Credible Intervals

Generates two-panel posterior predictive plots from MCMC samples:
1. **Malaria Incidence** (log scale, 2000–2050)
2. **K13 Resistance Ratio** (2000–2050)

---

## Files

| File | Description |
|---|---|
| `Plot_CI95.R` | Main plotting script — runs posterior predictive simulation and generates CI95 figures |
| `MCMC_out_example.rds` | Example MCMC output for testing and reproducing the plot |

---

## Workflow

```
../init_parameter_MCMC.R          # Fitted MCMC parameters
../init_state_2000.rds            # Equilibrium initial state (from Step 1.2)
MCMC_out_example.rds              # Posterior samples
        │
        ▼
Plot_CI95.R
   │
   ├── 1. Sample 90 parameter sets from posterior (start = 1000)
   │
   ├── 2. Run ODE for each parameter set
   │       ├── Annual incidence  → ysample   (Poisson posterior predictive)
   │       └── K13 ratio         → ysample2  (Beta posterior predictive, precision = φ)
   │
   ├── 3. Compute 95% credible intervals via getCredibleIntervals() (2.5%, 50%, 97.5%)
   │       ├── CI95             for incidence
   │       └── CI95_2           for resistance ratio
   │
   ├── 4. Plot Panel 1 — Malaria Incidence (log scale)
   └── 5. Plot Panel 2 — K13 Resistance Ratio
```

---

## Dependencies

```r
library(deSolve)        # ODE solver
library(BayesianTools)  # getSample(), getCredibleIntervals()
library(scales)         # alpha() for transparent colours
```

---

## Output Plots

### Panel 1 — Malaria Incidence
| Element | Description |
|---|---|
| Blue line | Model predicted median (50th percentile) |
| Grey band | 95% credible interval — Poisson posterior predictive |
| Open circles | Confirmed reported cases (2015–2023) |

### Panel 2 — K13 Resistance Ratio
| Element | Description |
|---|---|
| Red line | Model predicted median (50th percentile) |
| Grey band | 95% credible interval — Beta posterior predictive (parameterised by precision φ) |
| Open circles | Observed K13 allele frequency (2016–2022, 2025) |

---

## Key Parameters

| Parameter | Value | Description |
|---|---|---|
| `numSamples` | 90 | Posterior samples drawn per run (more samples = smoother CI) |
| `start` | 1000 | Burn-in — samples drawn after iteration 1000 |
| `times_fit` | `seq(1, 12×51)` | Monthly time steps, 51 years (2000–2050) |
| `start_b` | `12 × 21 = 252` | Month at which resistant strain begins competing |
| `phi` | `paramsample_test[i, 5]` | Beta precision — 5th column of MCMC chain |

Note: numSamples and start can be changed depending on the size of the MCMC output and convergence quality.
---