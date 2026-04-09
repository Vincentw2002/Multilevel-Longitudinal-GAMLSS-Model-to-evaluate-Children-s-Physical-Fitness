# Children's Fitness Assessment — Multilevel GAMLSS Pipeline & Dashboard

A complete R toolkit for analysing longitudinal children's fitness test data using **Generalised Additive Models for Location, Scale and Shape (GAMLSS)** in a three-level hierarchical structure (school → class/cohort → child). The repository contains a simulation engine for generating synthetic data and an interactive R Shiny dashboard aimed at school teachers and PE curriculum planners.

---

## Repository structure

| File | Role |
|---|---|
| `GAMLSS Model Architecture.R` | GAMLSS model architecture — fitting functions for all four fitness tests |
| `Helper Functions.R` | Diagnostic helpers sourced by the Shiny dashboard |
| `Synthetic Data Generator.R` | Synthetic data generator — baseline design and seven stress-test scenarios |
| `Shiny Dashboard.R` | Interactive dashboard for uploading real data, fitting models, and exploring results |

---

## Statistical approach

### Outcome models

Each fitness test is modelled with a distribution chosen to match the data-generating mechanism:

| Fitness test | Distribution | Rationale |
|---|---|---|
| **Sit-and-Reach** (cm) | Two-part hurdle: logistic + Log-Normal | Structural zeros (children who cannot reach past their toes) require a separate probability model |
| **Rope Skipping** (count) | Negative-Binomial (NBI) | Overdispersed count data with potential zero-inflation |
| **50m Dash** (seconds) | Log-Normal (LOGNO) | Right-skewed continuous positive outcome |
| **Vital Capacity** (mL) | Log-Normal (LOGNO) | Right-skewed continuous positive outcome |

### Multilevel random effects

All models share the same three-level random-effects structure on the **location (µ)** parameter:

```
School (intercept)
  └── Cohort / Class (intercept)
        └── Child (intercept + time slope)
```

The **scale (σ)** parameter additionally includes a school-level random intercept, allowing heteroskedasticity to vary across schools. Fixed effects on µ include a penalised spline of age (`pb(age)`), sex (0/1), and optionally BMI.

### Intraclass Correlation Coefficients (ICC)

The dashboard reports empirical (moment-based) ICCs at the school and school+cohort levels, decomposing total variance into the three hierarchical layers.

---

## Synthetic data generator (`Synthetic Data Generator.R`)

The script defines functions only — sourcing it loads them into the R environment without producing any data. To generate a dataset, call the functions explicitly after sourcing:

```r
source("Synthetic Data Generator.R")

# Generate baseline simulated data
df <- simulate_fitness_baseline(J = 10, H = 60, n = 50, T = 3, seed = 1)

# Optionally apply a stress-test scenario
df <- apply_scenarios_fitness(df, scenario = "S3", seed = 1)

# Save to CSV for use in the dashboard
write.csv(df, "simulated_data.csv", row.names = FALSE)
```

### Default design

| Parameter | Default | Description |
|---|---|---|
| J | 10 | Number of schools |
| H | 60 | Cohorts (classes) per school |
| n | 50 | Children per cohort |
| T | 3 | Annual measurement waves (time = 0, 1, 2) |

Children's age is drawn uniformly from [6, 10] at baseline and incremented by 1 year per wave, giving an observed age range of approximately 6–12 years. Measurements are **annual**.

### Shared latent fitness predictor

A shared latent variable η drives all four outcomes, partitioned via target ICCs:

- **ICC_school** = 0.01 (1% of variance between schools)
- **ICC_cohort** = 0.05 (5% between cohorts within school)
- **ICC_child** = 0.50 (50% between children within cohort)

Children also have a random slope on time correlated with the random intercept (default correlation −0.20), capturing regression-to-the-mean in repeated measurements.

### Stress-test scenarios

Seven scenarios probe model robustness under common real-world violations:

| Scenario | Description |
|---|---|
| **S1** | Baseline — no perturbation |
| **S2** | Strong floor effects in Sit-and-Reach (more structural zeros) |
| **S3** | Rope-skipping zero-inflation + higher overdispersion (σ = 1.2) |
| **S4** | Measurement heaping — 50m Dash rounded to nearest 0.5 s; Vital Capacity to nearest 50 mL |
| **S5** | Informative dropout — weaker children more likely to drop out at follow-up |
| **S6** | Cohort-level intervention shock — 10% of cohorts receive a performance boost from time 1 |
| **S7** | Device/administrator calibration bias across schools + 100% rounding |

---

## Shiny Dashboard (`Shiny Dashboard.R`)

### Setup

Place all four files in the **same folder**, rename `Shiny Dashboard.R` to `app.R`, then run:

```r
shiny::runApp()
```

#### Install required packages (run once)

```r
install.packages(c(
  "shiny", "shinydashboard", "shinyWidgets",
  "gamlss", "gamlss.dist",
  "dplyr", "ggplot2", "DT", "tidyr", "scales"
))
```

### Input CSV format

| Column | Type | Required | Notes |
|---|---|---|---|
| `school_id` | text | ✅ | School identifier |
| `cohort_id` | text | ✅ | Class / cohort identifier |
| `child_id` | text | ✅ | Unique child identifier |
| `time` | integer or label | ✅ | Wave index. Accepts `0, 1, 2` or strings like `"Year_1"`, `"Wave 2"`, `"T3"` — auto-converted to 0-based integers |
| `age` | numeric | ✅ | Age in years (e.g. `8.5`) |
| `sex` | 0/1 or text | ✅ | `1` / `"Male"` / `"M"` / `"Boy"` = male; `0` / `"Female"` / `"F"` = female |
| `bmi` | numeric | optional | Used as covariate if present and enabled in sidebar |
| `sitreach` | numeric | at least one outcome required | Sit-and-reach in cm; `0` = floor value |
| `ropeskip` | integer | | Rope-skipping jump count |
| `dash50` | numeric | | 50 m dash time in seconds |
| `vitalcap` | numeric | | Vital capacity in mL |

Multiple CSV files can be uploaded simultaneously and are row-bound automatically. A downloadable template is available inside the dashboard sidebar.

### Dashboard tabs

#### Tab 1 — Data Overview
- Data preview table and quality summary (N valid, % zeros, mean ± SD, range)
- Children observed per time point (missingness bar chart)
- Score distributions (histogram + kernel density per test)
- Age & sex profile
- **Score vs Age (by Sex)** — filterable by school and class; scatter of individual measurements with thin grey per-child connecting lines (longitudinal trajectories) overlaid with loess smooth ± 95% CI separately for males and females

#### Tab 2 — Model Fit Metrics
*(requires clicking "Fit Models" in the sidebar)*
- Convergence status and AIC per model
- Quantile residual summary table
- Prediction interval coverage vs nominal levels
- Centile calibration plot (empirical vs target — points on the diagonal indicate a well-calibrated model)
- Quantile residual histograms overlaid with the standard Normal curve

#### Tab 3 — Heterogeneity
- **Faceted boxplot** — top-N and bottom-N scoring classes per school, with a dashed school-average reference line; reveals within-school inequality at a glance
- School mean comparison dot plot (± 95% CI vs grand average)
- **Variance decomposition bar chart** — school-level (red) and class-level (orange) between-group variance plotted from zero, with individual-child variance annotated as grey text
- ICC summary table (empirical, unadjusted)

#### Tab 4 — Summary for Teachers
- **School selector** — all outputs in this tab filter reactively to the chosen school
- Header banner with key statistics (mean, SD, 10th/90th percentiles) per fitness test
- **Score vs Age (by Sex)** — same scatter + loess representation, with an optional class filter and a toggle to connect each child's repeated measurements
- **Class-level boxplots** — all classes in the selected school sorted by median, with a dashed school-wide average reference line
- **At-risk student table** — flags children whose score is >1.5 SD below their class mean (or above for 50m Dash, where lower is better); sortable and filterable, colour-coded by severity
- **Individual student lookup** — searchable dropdown of child IDs; highlights the selected child's measurements as large annotated points on a trajectory, overlaid on the school-wide scatter and loess context
- **Data-driven PE recommendations** — generated from the actual data; each recommendation (heterogeneity alert, flexibility deficit, coordination deficit, running speed, lung capacity, longitudinal monitoring, curriculum design) only appears when the corresponding threshold is exceeded in the selected school

---

## Intended users

| User | Primary use |
|---|---|
| **Statistical researchers** | Simulation engine for generating synthetic datasets under controlled design parameters and data-quality scenarios |
| **School teachers / PE staff** | Dashboard Tab 4 — plain-language summaries, at-risk identification, actionable PE recommendations |
| **School administrators** | Dashboard Tabs 1 & 3 — cross-school and cross-class comparisons, variance decomposition |

---

## Notes and limitations

- **Fitting time**: GAMLSS model fitting can take 5–20 minutes for datasets with many random-effect levels. Use the "Max GAMLSS cycles" slider to reduce iterations during exploratory analysis.
- **At-risk flagging**: The >1.5 SD threshold is unadjusted for age and sex. For formal referral decisions, age- and sex-standardised z-scores from the fitted GAMLSS model are preferable.
- **Empirical ICC**: The variance decomposition shown in Tab 3 is a moment-based approximation without covariate adjustment. Model-based ICCs derived from the fitted random effects are more precise but require the fitting step to be run first.
- **Quantile residuals**: For Log-Normal models (50m Dash, Vital Capacity, Sit-and-Reach), quantile residuals are computed analytically as (log y − µ̂) / σ̂. Because GAMLSS jointly optimises both µ and σ, these will always appear approximately N(0,1) by construction. The centile calibration plot is the more informative diagnostic for these outcomes.

---
