# MalariaPQ

An interactive R Shiny dashboard for simulating *Plasmodium falciparum* transmission dynamics and evaluating the impact of artemisinin-based combination therapy (ACT) with and without single low-dose primaquine (PQ) on gametocyte prevalence and artemisinin resistance propagation.

---

## Overview

MalariaPQ implements a deterministic compartmental ODE model that stratifies the human population by:

- **Clinical presentation** — asymptomatic vs. symptomatic infection
- **Treatment pathway** — no treatment (D0), ACT alone (D1), or ACT + low-dose primaquine (D2)
- **Parasite phenotype** — artemisinin-sensitive vs. artemisinin-resistant

The model explicitly tracks gametocyte-producing compartments to capture the infectious reservoir driving onward transmission, and incorporates seasonal forcing of the force of infection.

---

## App Structure

```
.
├── app.R                          # Main Shiny application (UI + server)
├── model/
│   └── model_D0_33_D1_67_baseline_change_birth_death_2betaSets_SS.R
├── data/
│   ├── init_state_d1_33.rds
│   ├── init_parameter_equilibrium_birth_2_rate_2beta_res_fitted_Tanzania.R
│   ├── Tanzania_incidence.csv
│   ├── Tanzania_pop_2000_2100.csv
│   ├── Tanzania_pop_SSP_2000_2100.csv
│   ├── Tanzania_inc1000_pop.csv
│   ├── Tanzania_K13_Allele_frequency.csv
│   └── Reported malaria cases by method of confirmation.csv
└── www/
    ├── css/style.css
    ├── js/toTheTop.js
    └── img/
        ├── figure1.png
        └── figure2.png
```

---

## Features

### Introduction tab
- Full model description with compartment diagrams
- ODE equations rendered with MathJax
- Parameter table

### Simulation tab
Interactive plots (baseline vs. D2 treatment scenario) covering:

| Plot | Description |
|------|-------------|
| Total Incidence & Population | Annual incidence fitted to WHO data; total population vs. UN projections |
| Fraction of Transmission & Resistant Fraction | F_T and proportion of resistant incident infections over time |
| Symptomatic / Asymptomatic Incidence | Absolute and proportional breakdown |
| Sensitive / Resistant Incidence | Absolute and proportional breakdown |
| Gametocyte Infections (Sym vs. Asym) | Gametocyte reservoir by clinical status |
| Gametocyte Infections (Sensitive vs. Resistant) | Gametocyte reservoir by parasite type |
| Treated / Untreated Channel Contribution | CT_total and CU_total over time |

### Sidebar controls
All parameters are adjustable via sliders:

- **Transmissibility** — β (symptomatic & asymptomatic) for two time periods; resistant multiplier; fitness cost
- **Symptomatic proportion**
- **Resistance** — start year of resistance emergence
- **Gametocyte rates** — γ for each compartment (sensitive & resistant)
- **Treatment** — D2 start year; clearance times for D0/D1/D2; % allocation across D0/D1/D2 (auto-constrained to sum to 100%); treatment failure rates

---

## Installation

```r
# Install required packages
install.packages(c(
  "shiny",
  "shinydashboard",
  "htmltools",
  "shinyBS",
  "deSolve",
  "shinyjs",
  "bsplus",
  "shinycssloaders"
))
```

Then run:

```r
shiny::runApp("path/to/MalariaPQ")
```

Or access the deployed version on shinyapps.io:  
**https://moru.shinyapps.io/MalariaPQ/**

---

## License

For academic and research use. Please contact the development team before redistribution.
