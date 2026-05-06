# MalariaPQ — Single Low-Dose Primaquine Model

![R](https://img.shields.io/badge/R-%3E%3D4.3.2-276DC3?style=flat&logo=r&logoColor=white)

An R-based pipeline for fitting a compartmental epidemiological model of malaria transmission with **single low-dose primaquine (SLD-PQ)** intervention.

---

## Repository Structure

```
.
├── 1.1 Fitting_population/                      # Fit model to population data
├── 1.2 Find_init_state_equilibrium_incidence/   # Find equilibrium initial state
├── 2. Fitting_Incidence/                        # Fit model to observed incidence
├── 3. Plot Result/                              # Generate output figures
├── Model/                                       # ODE system and parameter definitions
└── data/                                        # Input datasets
```

---

## Workflow Overview

Outputs from Steps 1.1 and 1.2 feed directly into Step 2.

```
Tanzania Population Data ──► 1.1 Fit Population 一 一 一 一 一 一 一 一 一 一
                                                                           │
                                                                           │
                                                                           │
Tanzania Incidence Data  ──► 1.2 Find Equilibrium Initial State            │
                             (fixed population size for tractable          │ fitted parameters
                              equilibrium incidence calculation)           │
                                      │                                    │
                                      │  equilibrium initial state         │
                                      ▼                                    │
Tanzania Incidence Data  ──►  2. Fit Incidence ◄一 一 一 一 一 一 一 一 一 一                            
and K13 Allele Frequency              │
                                      │
                                      ▼
                               3. Plot Results
```

### Step 1.1 — Fitting Population (`1.1 Fitting_population/`)
Fits the model to observed Tanzania population data (`Tanzania_pop_2000_2100.csv`) to estimate demographic parameters.

### Step 1.2 — Equilibrium Initial State (`1.2 Find_init_state_equilibrium_incidence/`)
Uses a **fixed (constant) population size** to analytically find the model's equilibrium state. Fixing population size simplifies solving for equilibrium incidence and avoids demographic drift confounding the calibration. The resulting equilibrium state vector is passed to Step 2 as the initial condition.

### Step 2 — Fitting Incidence (`2. Fitting_Incidence/`)
Fits the calibrated model to observed malaria incidence data (`Reported malaria cases by method of confirmation.csv`) and reported case data, incorporating K13 allele frequency trends (`Tanzania_K13_Allele_frequency.csv`). Generates parameter estimates and 95% confidence intervals (`2.Plot_CI95.R`).

### Step 3 — Plot Results (`3. Plot Result/`)
Produces publication-ready figures from model fits and projections.

---

## Data

| File | Description |
|---|---|
| `Reported malaria cases by method of confirmation.csv` | Observed malaria case counts stratified by diagnostic method |
| `Tanzania_K13_Allele_frequency.csv` | K13 propeller mutation allele frequencies (artemisinin partial resistance marker) |
| `Tanzania_pop_2000_2100.csv` | Tanzania population estimates and projections 2000–2100 |

> **Note:** SSP (Shared Socioeconomic Pathway) scenario files are included for demonstration purposes only.

---

## Requirements

- R >= 4.3.2 — "Eye Holes" (2023-10-31 ucrt) or later

Install required packages:

```r
install.packages(c("deSolve","dplyr", "ggplot2", "scales", "BayesianTools","openxlsx"))
```

---

## Getting Started
Clone the repository:
   ```bash
   git clone https://github.com/MAEMOD-MORU/MalariaPQ.git
   cd MalariaPQ
   ```

## Output

- Fitted demographic and transmission parameters
- Equilibrium model state (initial conditions for incidence fitting)
- Incidence fits with 95% confidence intervals
- Model projections with and without SLD-PQ intervention

---

## Related

- **MAEMOD-MORU** — [github.com/MAEMOD-MORU](https://github.com/MAEMOD-MORU)

---

## License

[MIT](LICENSE)
