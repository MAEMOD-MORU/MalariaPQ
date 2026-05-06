# Figure 3 — Scenario Comparison: Initial Resistance Ratio and D2 Coverage

Simulates malaria transmission under different **initial resistance ratios** and **D2 treatment coverages**, starting from population states where the resistance ratio has already reached a specified level (1%, 5%, 10%, 15%). Results are saved as Excel workbooks and plotted as resistance ratio trajectories over 10 years.

---

## Files

| File | Description |
|---|---|
| `find_inits_ratio.R` | Main script — finds initial states at target resistance ratios, runs D2 scenarios, plots and saves results |
| `results_example_0_xx.xlsx` | Example output — simulations starting from resistance ratio ≈ 1%,5%,10%,15% |

---

## Workflow

```
../init_parameter_MCMC.R        # Fitted model parameters from MCMC
../init_state_2000.rds          # Equilibrium initial state (from Step 1.2)
        │
        ▼
find_inits_ratio.R
   │
   ├── 1. Run ODE from year 2000 → compute rolling 12-month resistance ratio
   │
   ├── 2. Find positions closest to target ratios (1%, 5%, 10%, 15%)
   │       └── refine within ±5 positions (range_search = 5)
   │
   ├── 3. redistribute_init() — move asymptomatic resistant (Ar, GAr)
   │       into symptomatic compartments (Ir, GIr) proportionally
   │
   ├── 4. run_simulation_set() — for each initial state:
   │       ├── Baseline ( no D2)
   │       └── D2 ∈ {20%, 30%, 40%, 50%, 60%, 70%, 80%}
   │               └── 10-year projection
   │
   ├── 5. Plot resistance ratio trajectories (D2 30% vs 80%, all 4 init ratios)
   │
   └── 6. Save 4 Excel workbooks (one per initial resistance ratio)
```

---

## D2 Scenario Grid

| Parameter | Values |
|---|---|
| **Initial resistance ratio** | 1%, 5%, 10%, 15% |
| **D2 coverage** | 20%, 30%, 40%, 50%, 60%, 70%, 80% |
| **D0** | `1 − D2` |
| **D1** | Fixed at 0 |
| **Baseline** | D0 = 67%, D1 = 33%, D2 = 0% |
| **Simulation length** | 10 years (120 months) |

---

## Output Excel Files

Each workbook (`results_0_xx.xlsx`, etc.) contains **4 sheets**:

| Sheet | Description |
|---|---|
| `inc_r_model` | Annual resistant incidence |
| `inc_asym_r_model` | Annual asymptomatic resistant incidence |
| `inc_GR_model` | Annual resistant gametocyte incidence |
| `ratio_model` | Annual mean resistance ratio |

Each sheet has columns: `year`, `c_beta`, `Baseline (D1 33%)`, `D2 20%`, `D2 30%`, …, `D2 80%`

---

## Output Plot

Resistance ratio trajectories over 10 years for **D2 30%** and **D2 80%** across all 4 initial states:

| Colour | Initial resistance ratio |
|---|---|
| Black | 1% |
| Red | 5% |
| Green | 10% |
| Blue | 15% |

Dashed horizontal line at resistance ratio = 0.4 as a reference threshold.

---

## Dependencies

```r
library(deSolve)    # ODE solver
library(scales)     # alpha() for transparent colours
library(openxlsx)   # Excel workbook creation
```