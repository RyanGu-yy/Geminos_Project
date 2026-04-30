# =============================================================================
# simulate_data.R
# Nonlinear DID Project — Geminos Team
#
# PURPOSE: Generate synthetic panel datasets used in the DID method comparison.
#          Produces two scenarios:
#            - Scenario A: moderate nonlinearity (saved as synthetic_data.csv)
#            - Scenario B: strong selection + covariate interactions (stress test)
#                          (saved as stress_data.csv)
#
# OUTPUTS:
#   synthetic_data.csv   — baseline DGP (n = 2000)
#   stress_data.csv      — stress-test DGP (n = 2000)
#
# NOTATION:
#   n            : sample size
#   age          : continuous covariate, Uniform(20, 30)
#   genre_pref   : categorical covariate (Action / Drama / Comedy / SciFi / Romance)
#   major_release: categorical covariate (sports / politics / celebrity / tech / news)
#   genre_num    : integer encoding of genre_pref
#   release_num  : integer encoding of major_release
#   ps_true      : true propensity score P(D=1 | X), nonlinear in X
#   D            : binary treatment indicator, D ~ Bernoulli(ps_true)
#   Y_pre1       : outcome at t = -1 (pre-pre period)
#   Y_pre2       : outcome at t =  0 (pre period)
#   Y_post       : outcome at t =  1 (post period)
#   diff         : first difference Y_post - Y_pre2
#   tau          : individual treatment effect (tau_i = 1.5 in Scenario A)
#   ETT_true     : average tau_i among treated units — the estimand of interest
# =============================================================================

library(ggplot2)
library(tidyverse)
library(patchwork)

set.seed(42)


# -----------------------------------------------------------------------------
# SCENARIO A: Moderate nonlinearity
# -----------------------------------------------------------------------------

n <- 2000

# ── Covariates ──────────────────────────────────────────────────────────────
age           <- runif(n, 20, 30)
genre_pref    <- sample(c("Action", "Drama", "Comedy", "SciFi", "Romance"),
                        n, replace = TRUE)
major_release <- sample(c("sports", "politics", "celebrity", "tech", "news"),
                        n, replace = TRUE)

genre_num   <- as.numeric(factor(genre_pref))   # 1–5 integer encoding
release_num <- as.numeric(factor(major_release)) # 1–5 integer encoding

# ── Nonlinear propensity score P(D=1 | X) ───────────────────────────────────
# Logistic model with nonlinear and interaction terms.
# plogis() is the inverse-logit (sigmoid) function.
linpred <-
  0.5  * sin(age / 3) +
  0.3  * (genre_num^2) / 10 +
 -0.4  * log(release_num + 1) +
  0.2  * (age * genre_num) / 30

ps_true <- plogis(linpred)
D       <- rbinom(n, 1, ps_true)

# ── Pre-period outcomes (nonlinear in X) ────────────────────────────────────
Y_pre1 <-
   2   * sin(age / 4) +
   0.5 * genre_num^1.5 +
  -0.3 * release_num^2 / 5 +
   rnorm(n, 0, 2) + 30

Y_pre2 <-
   3   * sin(age / 5) +
   0.7 * genre_num^1.5 +
  -0.4 * release_num^2 / 5 +
   0.5 * age +
   rnorm(n, 0, 2) + 30

# ── Homogeneous treatment effect tau_i = 1.5 ────────────────────────────────
# Coefficients on sin(age/2) and genre indicators are set to 0 here;
# they can be changed to induce heterogeneous effects.
tau <- 1.5 +
  0 * sin(age / 2) +
  0 * (genre_num == 1) +   # Action
  0 * (genre_num == 2)     # Drama

ETT_true <- mean(tau[D == 1])
cat("Scenario A — True ETT (treated units only):", round(ETT_true, 4), "\n")

# ── Post-period outcome ──────────────────────────────────────────────────────
# The nonlinear term 5*sin(age/3) represents time-varying confounding that
# affects both groups equally — it does NOT violate conditional parallel trends.
Y_post <-
  Y_pre2 +
  tau * D +
  5 * sin(age / 3) +
  rnorm(n, 0, 2)

# ── Assemble and save ────────────────────────────────────────────────────────
data <- data.frame(
  age, genre_pref, major_release,
  genre_num, release_num,
  D, ps_true,
  Y_pre1, Y_pre2, Y_post,
  diff = Y_post - Y_pre2,
  tau
)

write.csv(data, "synthetic_data.csv", row.names = FALSE)
cat("Saved: synthetic_data.csv\n\n")

# ── Exploratory plots ────────────────────────────────────────────────────────
p1 <- ggplot(data, aes(x = age, y = Y_post, color = factor(D))) +
  geom_point(alpha = 0.15, size = 0.8) +
  geom_smooth(se = FALSE, linewidth = 1.2) +
  scale_color_manual(values = c("#4477AA", "#EE6677"),
                     labels = c("Control", "Treated")) +
  labs(title = "Scenario A: Post-outcome vs Age", color = "Group") +
  theme_minimal()

p2 <- ggplot(data, aes(x = age, y = diff, color = factor(D))) +
  geom_point(alpha = 0.15, size = 0.8) +
  geom_smooth(se = FALSE, linewidth = 1.2) +
  scale_color_manual(values = c("#4477AA", "#EE6677"),
                     labels = c("Control", "Treated")) +
  labs(title = "First Difference (Y_post - Y_pre2) vs Age", color = "Group") +
  theme_minimal()

p3 <- ggplot(data, aes(x = age, y = tau)) +
  geom_point(alpha = 0.2, size = 0.8, color = "#228833") +
  geom_smooth(se = FALSE, color = "#228833", linewidth = 1.2) +
  labs(title = "True Individual Treatment Effect tau_i vs Age") +
  theme_minimal()

print(p1 + p2 + p3)


# -----------------------------------------------------------------------------
# SCENARIO B: Strong selection + nonlinear covariate interactions (stress test)
# -----------------------------------------------------------------------------
# This DGP amplifies violations of unconditional parallel trends so that the
# performance gap between naive DID and conditional/IPW/DR methods is starker.

n_stress <- 2000

age_s           <- runif(n_stress, 20, 30)
genre_pref_s    <- sample(c("Action", "Drama", "Comedy", "SciFi", "Romance"),
                          n_stress, replace = TRUE)
major_release_s <- sample(c("sports", "politics", "celebrity", "tech", "news"),
                          n_stress, replace = TRUE)

genre_num_s   <- as.numeric(factor(genre_pref_s))
release_num_s <- as.numeric(factor(major_release_s))

# Stronger, more complex propensity score with interaction terms
linpred_s <-
  -1.2 +
   1.1 * cos(age_s / 2) +
   0.9 * (genre_pref_s == "Action") +
   0.7 * (genre_pref_s == "SciFi") -
   0.8 * (major_release_s == "politics") +
   1.0 * ((age_s < 24) & (genre_pref_s == "Action")) +
   0.8 * ((age_s > 27) & (major_release_s == "tech"))

ps_true_s <- plogis(linpred_s)
D_s       <- rbinom(n_stress, 1, ps_true_s)

Y_pre1_s <-
  1.5 * sin(age_s / 2.5) +
  0.7 * genre_num_s -
  0.4 * release_num_s +
  0.9 * sin(age_s * genre_num_s / 8) +
  rnorm(n_stress, 0, 2) + 28

Y_pre2_s <-
  Y_pre1_s +
  1.2 * cos(age_s / 2.2) +
  0.5 * (genre_pref_s == "Action") -
  0.6 * (major_release_s == "politics") +
  0.8 * sin(age_s * release_num_s / 7)

tau_s <-
  1.5 +
  0 * sin(age_s / 1.8) +
  0 * (genre_pref_s == "Action") -
  0 * (genre_pref_s == "Drama")

# Highly nonlinear untreated counterfactual trend (violates linear parallel trends)
m0_stress <-
  2.8 * sin(age_s / 1.7) +
  1.1 * (genre_pref_s == "Action") -
  1.0 * (genre_pref_s == "Drama") +
  1.3 * (major_release_s == "tech") -
  1.1 * (major_release_s == "politics") +
  1.4 * sin(age_s * genre_num_s / 3.5) +
  1.0 * ((age_s > 26) & (major_release_s == "tech")) +
  rnorm(n_stress, 0, 1.5)

Y_post_s <- Y_pre2_s + tau_s * D_s + m0_stress

stress_data <- data.frame(
  age          = age_s,
  genre_pref   = genre_pref_s,
  major_release = major_release_s,
  genre_num    = genre_num_s,
  release_num  = release_num_s,
  D            = D_s,
  ps_true      = ps_true_s,
  Y_pre1       = Y_pre1_s,
  Y_pre2       = Y_pre2_s,
  Y_post       = Y_post_s,
  diff         = Y_post_s - Y_pre2_s,
  tau          = tau_s
)

write.csv(stress_data, "stress_data.csv", row.names = FALSE)
cat("Saved: stress_data.csv\n")
