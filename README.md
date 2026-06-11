# MalariaPQ
![R](https://img.shields.io/badge/R-%3E%3D4.3.2-276DC3?style=flat&logo=r&logoColor=white)

An interactive R Shiny dashboard for simulating *Plasmodium falciparum* transmission dynamics and evaluating the impact of artemisinin-based combination therapy (ACT) with and without single low-dose primaquine (PQ) on gametocyte prevalence and artemisinin resistance propagation.

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
        └── figure1.png
```
---

## Simulation tab
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

## Sidebar controls
All parameters are adjustable via sliders:

- **Transmissibility** — β (symptomatic & asymptomatic) for two time periods; resistant multiplier; fitness cost
- **Symptomatic proportion**
- **Resistance** — start year of resistance emergence
- **Gametocyte rates** — γ for each compartment (sensitive & resistant)
- **Treatment** — D2 start year; clearance times for D0/D1/D2; % allocation across D0/D1/D2 (auto-constrained to sum to 100%); treatment failure rates

---

## Requirements

- R >= 4.3.2 — "Eye Holes" (2023-10-31 ucrt) or later

Install required packages:
```r
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
