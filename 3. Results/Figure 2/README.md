# Figure 2 — D2 Scenario Comparison Table

Generates a multi-sheet Excel workbook comparing model projections across **D2 treatment coverage scenarios** and **asymptomatic proportions**, starting from year 2025. Each sheet reports the impact on resistant incidence and resistant gametocyte carriers at 1-year and 3-year windows after D2 introduction.

---

## Files

| File | Description |
|---|---|
| `create_table_D2.R` | Runs ODE simulations across D2 and asymptomatic scenarios; writes results to Excel |
| `Scenario_tables_by_D2.xlsx` | Output workbook — one sheet per D2 value (D2_0.0 to D2_1.0) |

---

## Workflow

```
                ../init_parameter_MCMC.R       # Fitted model parameters from MCMC
                ../init_state_2000.rds         # Equilibrium initial state (from Step 1.2)
                        │
                        ▼
                create_table_D2.R
                        │
                        ▼
                Scenario_tables_by_D2.xlsx
```

---

## Scenario Grid

| Axis | Values |
|---|---|
| **D2** (proportion treated with D2 regimen) | 0.0, 0.1, 0.2, …, 1.0 |
| **Asym** (% asymptomatic infections) | 0, 10, 20, …, 100 |
| **D1** | Fixed at 0 |
| **D0** | `1 − D2` |
| **D2 start year** | 2025 (`start_d = 12 × 26`) |

---

## Output Columns (per sheet)

| Column | Description |
|---|---|
| `Asym` | % asymptomatic (`"Baseline"` for first row) |
| `c_beta` | Competitive fitness of resistant strain |
| `D0`, `D1`, `D2` | Treatment proportion for each day-regimen |
| `Incidence_Res_year_1` | Annual resistant incidence at year 2026 |
| `Incidence_Res_year_3` | Annual resistant incidence at year 2028 |
| `red_inc_r_1y` | % reduction in resistant incidence vs baseline (1-year) |
| `red_inc_r_3y` | % reduction in resistant incidence vs baseline (3-year) |
| `Gametocyte_Res_year_1` | Resistant gametocyte incidence at year 2026 |
| `Gametocyte_Res_year_3` | Resistant gametocyte incidence at year 2028 |
| `red_GR_inc_1y` | % reduction in resistant gametocytes vs baseline (1-year) |
| `red_GR_inc_3y` | % reduction in resistant gametocytes vs baseline (3-year) |

> The **first row of each sheet is the baseline** (D2 = 0, no D2 treatment) for reference comparison.

---

## Dependencies

```r
library(deSolve)    # ODE solver
library(openxlsx)   # Excel workbook creation
```

---



