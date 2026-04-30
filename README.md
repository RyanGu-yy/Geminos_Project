# Nonlinear DID: Methods & Theory Extension

**Geminos Team · April 2026**

This repository accompanies the report *Nonlinear DID: Methods & Theory Extension* and provides fully reproducible code for generating synthetic panel data and evaluating five difference-in-differences estimators under nonlinear data-generating processes.

---

## Repository Structure

```
Geminos_repo/
├── simulate_data.R              # Generate synthetic datasets
├── run_models.R                 # Estimate ETT using all DID methods
├── synthetic_data.csv           # Scenario A output (produced by simulate_data.R)
├── stress_data.csv              # Scenario B output (produced by simulate_data.R)
└── Geminos_NonlinearDID_Extended.qmd   # Full report (Quarto source)
```

---

## How to Run

**Step 1 — Generate synthetic data**

```r
source("simulate_data.R")
```

Produces `synthetic_data.csv` (Scenario A) and `stress_data.csv` (Scenario B).

**Step 2 — Run DID estimators**

```r
source("run_models.R")
```

Prints results tables and displays comparison plots for both scenarios.

---

## Data Description

Both datasets contain **n = 2000** simulated individuals with the following columns:

| Column | Type | Description |
|--------|------|-------------|
| `age` | numeric | Age, drawn from Uniform(20, 30) |
| `genre_pref` | character | Categorical covariate: Action / Drama / Comedy / SciFi / Romance |
| `major_release` | character | Categorical covariate: sports / politics / celebrity / tech / news |
| `genre_num` | integer | Integer encoding of `genre_pref` |
| `release_num` | integer | Integer encoding of `major_release` |
| `D` | integer | Binary treatment indicator (1 = treated, 0 = control) |
| `ps_true` | numeric | True propensity score P(D=1 \| X) |
| `Y_pre1` | numeric | Outcome at t = −1 (pre-pre period) |
| `Y_pre2` | numeric | Outcome at t = 0 (pre period) |
| `Y_post` | numeric | Outcome at t = 1 (post period) |
| `diff` | numeric | First difference: Y_post − Y_pre2 |
| `tau` | numeric | Individual treatment effect τᵢ |

The estimand of interest is the **Effect of Treatment on the Treated (ETT)**:

$$\text{ETT} = \mathbb{E}[\tau_i \mid D_i = 1]$$

---

## DGP Design

### Scenario A — Moderate Nonlinearity (`synthetic_data.csv`)

- Propensity score is nonlinear in covariates: `plogis(0.5·sin(age/3) + ...)`
- Outcome trends are nonlinear: `Y_pre ~ sin(age) + genre_num^1.5 + ...`
- Treatment effect is homogeneous: τᵢ = 1.5 for all units
- Unconditional parallel trends holds approximately but not exactly

### Scenario B — Strong Selection + Interactions (`stress_data.csv`)

- Stronger propensity score with age × genre interaction terms
- Highly nonlinear untreated counterfactual trend (large sin/cos components)
- Designed so that naive unconditional DID fails visibly, making the advantage of conditional/IPW/DR methods clear

---

## Methods

### 1. Unconditional DID (Naive)

Simple two-way difference with no covariate adjustment:

$$\widehat{\text{ETT}}_\text{naive} = (\bar{Y}_{\text{post},1} - \bar{Y}_{\text{pre},1}) - (\bar{Y}_{\text{post},0} - \bar{Y}_{\text{pre},0})$$

Biased when treated and control groups differ in covariates that affect outcome trends.

### 2. Linear Conditional DID

Fits an OLS model for E[ΔY | X] on controls, then subtracts predicted counterfactual trends from the treated mean. Misspecified when the true trend μ₀(X) is nonlinear.

### 3. GAM Conditional DID

Replaces the linear model with a Generalized Additive Model (GAM):

$$\mathbb{E}[\Delta Y_i \mid X_i, D_i=0] = \alpha + f_1(\text{age}_i) + f_2(g_i) + f_3(r_i)$$

where f₁ is estimated as a cubic regression spline. Consistent under the Nonlinear Conditional Parallel Trends (NCPT) assumption.

### 4. IPW-DID — Abadie (2005)

Estimates the propensity score π(X) = P(D=1|X) via GAM, then applies the Abadie (2005, eq. 10) inverse-probability-weighted estimator:

$$\widehat{\text{ETT}}_\text{IPW} = \frac{1}{\hat{P}(D=1)} \mathbb{E}\left[\Delta Y_i \cdot \frac{D_i - \hat{\pi}(X_i)}{1 - \hat{\pi}(X_i)}\right]$$

Propensity scores are trimmed to [0.02, 0.98]. Consistent if the propensity score model is correctly specified; no outcome model needed.

### 5. Doubly Robust DID

Combines outcome regression and IPW. Consistent if **either** the propensity score model **or** the outcome model is correctly specified (but not necessarily both):

$$\widehat{\text{ETT}}_\text{DR} = \frac{1}{n}\sum_i (w^{\text{treat}}_i - w^{\text{ctrl}}_i)(\Delta Y_i - \hat{m}_0(X_i))$$

where $w^{\text{treat}}_i = D_i / \bar{D}$ and $w^{\text{ctrl}}_i$ are normalized IPW weights.

### 6. Matching DID (Mahalanobis)

Matches each treated unit to its nearest control unit using Mahalanobis distance on (age, genre_num, release_num), 1:1 without replacement. ETT is estimated as the average difference in first-differences across matched pairs.

---

## Dependencies

```r
install.packages(c("ggplot2", "tidyverse", "mgcv", "MatchIt", "patchwork", "knitr"))
```

---

## Key Result

Under both scenarios, **GAM Conditional DID**, **IPW-DID**, and **Doubly Robust DID** substantially outperform naive and linear DID. The advantage is most pronounced in Scenario B, where the untreated counterfactual trend is highly nonlinear and the covariate distributions are strongly unbalanced between treated and controls.

---

## References

- Abadie, A. (2005). Semiparametric difference-in-differences estimators. *Review of Economic Studies*, 72(1), 1–19.
- Callaway, B., & Sant'Anna, P. H. C. (2021). Difference-in-differences with multiple time periods. *Journal of Econometrics*, 225(2), 200–230.
- Chernozhukov, V., et al. (2018). Double/debiased machine learning for treatment and structural parameters. *The Econometrics Journal*, 21(1), C1–C68.
- Sant'Anna, P. H. C., & Zhao, J. (2020). Doubly robust difference-in-differences estimators. *Journal of Econometrics*, 219(1), 101–122.
- Roth, J., Sant'Anna, P. H. C., Bilinski, A., & Poe, J. (2023). What's trending in difference-in-differences? *Journal of Econometrics*, 235(2), 2218–2244.
