# =============================================================================
# run_models.R
# Nonlinear DID Project — Geminos Team
#
# PURPOSE: Load synthetic data and estimate the Effect of Treatment on the
#          Treated (ETT) using five DID methods. Compares performance across
#          two DGP scenarios (Scenario A: moderate nonlinearity;
#          Scenario B: strong selection + interactions).
#
# PREREQUISITES: Run simulate_data.R first to produce:
#   synthetic_data.csv   — Scenario A dataset
#   stress_data.csv      — Scenario B dataset
#
# METHODS IMPLEMENTED:
#   1. Unconditional DID       — naive two-way difference, no covariates
#   2. Linear Conditional DID  — OLS outcome regression on controls (misspecified)
#   3. GAM Conditional DID     — nonparametric outcome regression via GAM
#   4. IPW-DID (Abadie 2005)   — inverse-probability weighting on propensity score
#   5. Doubly Robust DID       — combines outcome regression + IPW
#   6. Matching DID            — Mahalanobis nearest-neighbour matching
#
# NOTATION:
#   D          : binary treatment indicator (1 = treated, 0 = control)
#   X          : baseline covariates (age, genre_pref, major_release)
#   delta_Y    : first difference Y_post - Y_pre2  (stored as `diff` in data)
#   pi(X)      : propensity score P(D=1 | X), estimated via GAM
#   pi_trim    : trimmed propensity score, clipped to [0.02, 0.98]
#   m0(X)      : E[delta_Y | X, D=0], estimated outcome model for controls
#   ETT_true   : average treatment effect on treated, computed from tau column
#   ETT_*      : estimated ETT under each method
# =============================================================================

library(ggplot2)
library(tidyverse)
library(mgcv)      # GAM estimation
library(MatchIt)   # nearest-neighbour matching
library(patchwork)
library(knitr)

# -----------------------------------------------------------------------------
# Load data
# -----------------------------------------------------------------------------
data        <- read.csv("synthetic_data.csv")
stress_data <- read.csv("stress_data.csv")


# -----------------------------------------------------------------------------
# Helper: estimate_did_suite(data)
#
# Applies all six estimators to a single dataset and returns a tibble of
# (Estimator, Estimate, Bias) rows. Used for both Scenario A and Scenario B.
# -----------------------------------------------------------------------------
estimate_did_suite <- function(data) {

  ETT_true <- mean(data$tau[data$D == 1])

  # ── 1. Unconditional DID ─────────────────────────────────────────────────
  # Two-way difference with no covariate adjustment.
  # Biased when treated/control covariate distributions differ.
  ETT_naive <-
    (mean(data$Y_post[data$D == 1]) - mean(data$Y_pre2[data$D == 1])) -
    (mean(data$Y_post[data$D == 0]) - mean(data$Y_pre2[data$D == 0]))

  # ── 2. Linear Conditional DID ────────────────────────────────────────────
  # Fit a linear model for E[delta_Y | X] on controls, then subtract the
  # predicted counterfactual trend from the treated mean difference.
  # Misspecified here because the true trend mu_0(X) is nonlinear.
  ctrl_linear <- lm(
    diff ~ age + genre_pref + major_release,
    data = subset(data, D == 0)
  )
  cf_linear  <- predict(ctrl_linear, newdata = subset(data, D == 1))
  ETT_linear <- mean(data$diff[data$D == 1]) - mean(cf_linear)

  # ── 3. GAM Conditional DID ───────────────────────────────────────────────
  # Replace the linear model with a GAM:
  #   E[delta_Y | X, D=0] = alpha + f1(age) + f2(genre) + f3(release)
  # where f1 is a cubic regression spline. Consistent under NCPT.
  ctrl_gam <- gam(
    diff ~ s(age) + genre_pref + major_release,
    data = subset(data, D == 0)
  )
  cf_gam  <- predict(ctrl_gam, newdata = subset(data, D == 1))
  ETT_gam <- mean(data$diff[data$D == 1]) - mean(cf_gam)

  # ── 4. IPW-DID (Abadie 2005) ─────────────────────────────────────────────
  # Estimate propensity score pi(X) = P(D=1|X) via GAM, then apply
  # the Abadie (2005, eq. 10) IPW estimator:
  #
  #   ETT_IPW = (1 / E[D]) * E[ delta_Y * (D - pi(X)) / (1 - pi(X)) ]
  #
  # Propensity scores are trimmed to [0.02, 0.98] for numerical stability.
  ps_gam <- gam(
    D ~ s(age) + genre_pref + major_release,
    family = binomial,
    data   = data
  )
  ps_hat  <- predict(ps_gam, type = "response")
  ps_trim <- pmax(0.02, pmin(0.98, ps_hat))
  pD_bar  <- mean(data$D)

  ETT_ipw <- mean(
    (data$diff / pD_bar) *
    (data$D - ps_trim) /
    (1 - ps_trim)
  )

  # ── 5. Doubly Robust DID ─────────────────────────────────────────────────
  # Combines outcome regression (GAM) and IPW. Consistent if EITHER the
  # propensity score model OR the outcome model is correctly specified.
  #
  # Normalized DR estimator for panel ETT:
  #   w_treat_i  = D_i / D_bar
  #   w_ctrl_i   = [(1-D_i)*pi(X_i)/(1-pi(X_i))] / mean(same)
  #   ETT_DR     = mean( (w_treat - w_ctrl) * (delta_Y - m0_hat) )
  m0_hat     <- predict(ctrl_gam, newdata = data)   # reuse GAM from step 3
  w_treat    <- data$D / pD_bar
  w_ctrl_raw <- (1 - data$D) * ps_trim / (1 - ps_trim)
  w_ctrl     <- w_ctrl_raw / mean(w_ctrl_raw)

  ETT_dr <- mean((w_treat - w_ctrl) * (data$diff - m0_hat))

  # ── 6. Matching DID (Mahalanobis) ────────────────────────────────────────
  # Match each treated unit to its nearest control unit using Mahalanobis
  # distance on (age, genre_num, release_num), 1:1 without replacement.
  m_out <- matchit(
    D ~ age + genre_num + release_num,
    data     = data,
    method   = "nearest",
    distance = "mahalanobis",
    ratio    = 1
  )
  matched_data <- match.data(m_out)
  ETT_match    <- with(matched_data,
    mean(diff[D == 1]) - mean(diff[D == 0])
  )

  # ── Results tibble ────────────────────────────────────────────────────────
  tibble(
    Estimator = c(
      "True ETT",
      "Unconditional DID",
      "Linear Conditional DID",
      "GAM Conditional DID",
      "IPW-DID (Abadie 2005)",
      "Doubly Robust DID",
      "Matching DID (Mahalanobis)"
    ),
    Estimate = c(ETT_true, ETT_naive, ETT_linear, ETT_gam,
                 ETT_ipw,  ETT_dr,    ETT_match),
    Bias = c(0,
             ETT_naive  - ETT_true,
             ETT_linear - ETT_true,
             ETT_gam    - ETT_true,
             ETT_ipw    - ETT_true,
             ETT_dr     - ETT_true,
             ETT_match  - ETT_true)
  )
}


# -----------------------------------------------------------------------------
# Scenario A results
# -----------------------------------------------------------------------------
cat("Running Scenario A estimators...\n")
results_a <- estimate_did_suite(data) %>%
  mutate(Scenario = "Scenario A: moderate nonlinearity")

cat("\nScenario A results:\n")
print(knitr::kable(select(results_a, Estimator, Estimate, Bias),
                   caption = "Scenario A — ETT Estimates"))

# ── Single-scenario bar chart ────────────────────────────────────────────────
results_a %>%
  filter(Estimator != "True ETT") %>%
  mutate(Estimator = fct_reorder(Estimator, abs(Bias))) %>%
  ggplot(aes(x = Estimator, y = Estimate, fill = abs(Bias) < 0.1)) +
  geom_col(alpha = 0.85) +
  geom_hline(
    yintercept = results_a$Estimate[results_a$Estimator == "True ETT"],
    linetype = "dashed", color = "red", linewidth = 1
  ) +
  scale_fill_manual(values = c("#CC3311", "#009988"),
                    labels = c("|Bias| >= 0.1", "|Bias| < 0.1")) +
  coord_flip() +
  labs(
    title    = "Scenario A: DID Estimator Comparison",
    subtitle = paste("True ETT =",
                     round(results_a$Estimate[results_a$Estimator == "True ETT"], 4)),
    y = "ETT Estimate", x = NULL, fill = ""
  ) +
  theme_minimal()


# -----------------------------------------------------------------------------
# Scenario B results
# -----------------------------------------------------------------------------
cat("\nRunning Scenario B estimators (stress test)...\n")
results_b <- estimate_did_suite(stress_data) %>%
  mutate(Scenario = "Scenario B: strong selection + interactions")

cat("\nScenario B results:\n")
print(knitr::kable(select(results_b, Estimator, Estimate, Bias),
                   caption = "Scenario B — ETT Estimates"))


# -----------------------------------------------------------------------------
# Combined scenario comparison
# -----------------------------------------------------------------------------
scenario_compare <- bind_rows(results_a, results_b) %>%
  mutate(
    Estimate = round(Estimate, 4),
    Bias     = round(Bias, 4)
  )

cat("\nCombined scenario comparison:\n")
print(knitr::kable(scenario_compare,
                   caption = "Two DGPs: Baseline vs Stress-Test Scenario"))

# ── Side-by-side bias plot ───────────────────────────────────────────────────
scenario_compare %>%
  filter(Estimator != "True ETT") %>%
  ggplot(aes(x = Estimator, y = Bias, fill = Scenario)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray30") +
  coord_flip() +
  labs(
    title    = "Bias Comparison Across Two Nonlinear DID Scenarios",
    subtitle = "Scenario B is designed so that unconditional DID fails more visibly",
    x = NULL,
    y = "Bias (Estimate - True ETT)"
  ) +
  theme_minimal()
