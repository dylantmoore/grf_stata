#!/usr/bin/env Rscript
# generate_reference.R -- Generate reference data for grf Stata plugin tests
#
# Requires: grf (pinned version), MASS
# Install: install.packages("grf", repos="https://cloud.r-project.org")
#
# Usage: Rscript generate_reference.R

suppressPackageStartupMessages({
  library(grf)
  library(MASS)
})

cat("grf version:", as.character(packageVersion("grf")), "\n")

set.seed(42)
outdir <- "ref"
dir.create(outdir, showWarnings = FALSE)

# ============================================================
# Helper: save data frame as CSV
# ============================================================
save_csv <- function(df, name) {
  fpath <- file.path(outdir, paste0(name, ".csv"))
  write.csv(df, fpath, row.names = FALSE)
  cat("  Saved:", fpath, "(", nrow(df), "x", ncol(df), ")\n")
}

# ============================================================
# 1. REGRESSION FOREST
# ============================================================
cat("\n=== Regression Forest ===\n")

n <- 500
p <- 5
X <- matrix(rnorm(n * p), n, p)
colnames(X) <- paste0("x", 1:p)
Y <- 2 * X[, 1] + X[, 2]^2 + 0.5 * rnorm(n)

# Save input data
input_df <- data.frame(X, y = Y)
save_csv(input_df, "regression_input")

# Fit regression forest
t0 <- proc.time()
rf <- regression_forest(X, Y,
  num.trees = 2000,
  seed = 42,
  honesty = TRUE,
  honesty.fraction = 0.5,
  min.node.size = 5,
  alpha = 0.05,
  ci.group.size = 2
)
t1 <- proc.time()
cat("  R training time:", (t1 - t0)[3], "seconds\n")

# OOB predictions with variance
pred <- predict(rf, estimate.variance = TRUE)
output_df <- data.frame(
  prediction = pred$predictions,
  variance = pred$variance.estimates
)
save_csv(output_df, "regression_output")

# Variable importance
vi <- variable_importance(rf)
vi_df <- data.frame(variable = paste0("x", 1:p), importance = as.vector(vi))
save_csv(vi_df, "regression_variable_importance")

cat("  Prediction mean:", mean(pred$predictions), "\n")
cat("  Prediction SD:", sd(pred$predictions), "\n")
cat("  Correlation(Y, pred):", cor(Y, pred$predictions), "\n")

# ============================================================
# 2. CAUSAL FOREST
# ============================================================
cat("\n=== Causal Forest ===\n")

n <- 500
p <- 5
X <- matrix(rnorm(n * p), n, p)
colnames(X) <- paste0("x", 1:p)
W <- rbinom(n, 1, 0.5)
tau <- X[, 1]  # heterogeneous treatment effect
Y <- 2 * X[, 1] + tau * W + 0.5 * rnorm(n)

input_df <- data.frame(X, y = Y, w = W)
save_csv(input_df, "causal_input")

t0 <- proc.time()
cf <- causal_forest(X, Y, W,
  num.trees = 2000,
  seed = 42,
  honesty = TRUE,
  min.node.size = 5
)
t1 <- proc.time()
cat("  R training time:", (t1 - t0)[3], "seconds\n")

# CATE predictions
pred_cf <- predict(cf, estimate.variance = TRUE)
output_df <- data.frame(
  cate = pred_cf$predictions,
  variance = pred_cf$variance.estimates
)
save_csv(output_df, "causal_output")

# ATE
ate <- average_treatment_effect(cf)
ate_df <- data.frame(estimate = ate[1], std_err = ate[2])
save_csv(ate_df, "causal_ate")

# Variable importance
vi_cf <- variable_importance(cf)
vi_df <- data.frame(variable = paste0("x", 1:p), importance = as.vector(vi_cf))
save_csv(vi_df, "causal_variable_importance")

# Best linear projection
# Returns a coeftest object (matrix with cols: Estimate, Std. Error, t value, Pr(>|t|))
blp <- best_linear_projection(cf, X)
blp_df <- data.frame(
  variable = c("intercept", paste0("x", 1:p)),
  estimate = blp[, 1],
  std_error = blp[, 2],
  t_value = blp[, 3],
  p_value = blp[, 4]
)
save_csv(blp_df, "causal_blp")

# Save X from section 2 before it gets overwritten by later sections
X_causal <- X

# Test calibration
tc <- test_calibration(cf)
tc_df <- data.frame(
  term = rownames(tc),
  estimate = tc[, 1],
  std_error = tc[, 2],
  t_value = tc[, 3],
  p_value = tc[, 4]
)
save_csv(tc_df, "causal_test_calibration")

cat("  CATE mean:", mean(pred_cf$predictions), "\n")
cat("  ATE:", ate[1], "(SE:", ate[2], ")\n")
cat("  Correlation(tau, CATE):", cor(tau, pred_cf$predictions), "\n")

# ============================================================
# 3. QUANTILE FOREST
# ============================================================
cat("\n=== Quantile Forest ===\n")

n <- 500
p <- 5
X <- matrix(rnorm(n * p), n, p)
colnames(X) <- paste0("x", 1:p)
Y <- 2 * X[, 1] + X[, 2] + rnorm(n) * (1 + abs(X[, 3]))  # heteroscedastic

input_df <- data.frame(X, y = Y)
save_csv(input_df, "quantile_input")

quantiles <- c(0.1, 0.25, 0.5, 0.75, 0.9)

t0 <- proc.time()
qf <- quantile_forest(X, Y,
  quantiles = quantiles,
  num.trees = 2000,
  seed = 42,
  honesty = TRUE,
  min.node.size = 5
)
t1 <- proc.time()
cat("  R training time:", (t1 - t0)[3], "seconds\n")

pred_qf <- predict(qf, quantiles = quantiles)
output_df <- as.data.frame(pred_qf$predictions)
colnames(output_df) <- gsub("\\.", "_", paste0("q", quantiles))
save_csv(output_df, "quantile_output")

cat("  Median prediction mean:", mean(output_df$q0.5), "\n")

# ============================================================
# 4. INSTRUMENTAL FOREST
# ============================================================
cat("\n=== Instrumental Forest ===\n")

n <- 500
p <- 5
X <- matrix(rnorm(n * p), n, p)
colnames(X) <- paste0("x", 1:p)
Z <- rbinom(n, 1, 0.5)                     # instrument
W <- Z + 0.2 * rnorm(n) > 0.5             # endogenous treatment
W <- as.numeric(W)
tau <- X[, 1]
Y <- 2 * X[, 1] + tau * W + 0.5 * rnorm(n)

input_df <- data.frame(X, y = Y, w = W, z = Z)
save_csv(input_df, "instrumental_input")

t0 <- proc.time()
ivf <- instrumental_forest(X, Y, W, Z,
  num.trees = 2000,
  seed = 42,
  honesty = TRUE,
  min.node.size = 5
)
t1 <- proc.time()
cat("  R training time:", (t1 - t0)[3], "seconds\n")

pred_ivf <- predict(ivf, estimate.variance = TRUE)
output_df <- data.frame(
  late = pred_ivf$predictions,
  variance = pred_ivf$variance.estimates
)
save_csv(output_df, "instrumental_output")

cat("  LATE mean:", mean(pred_ivf$predictions), "\n")
cat("  Correlation(tau, LATE):", cor(tau, pred_ivf$predictions), "\n")

# ============================================================
# 5. PROBABILITY FOREST
# ============================================================
cat("\n=== Probability Forest ===\n")

n <- 500
p <- 5
X <- matrix(rnorm(n * p), n, p)
colnames(X) <- paste0("x", 1:p)
probs <- plogis(2 * X[, 1] + X[, 2])
Y_class <- factor(rbinom(n, 2, probs), levels = 0:2)

input_df <- data.frame(X, y = as.numeric(as.character(Y_class)))
save_csv(input_df, "probability_input")

t0 <- proc.time()
pf <- probability_forest(X, Y_class,
  num.trees = 2000,
  seed = 42,
  honesty = TRUE,
  min.node.size = 5
)
t1 <- proc.time()
cat("  R training time:", (t1 - t0)[3], "seconds\n")

pred_pf <- predict(pf)
output_df <- as.data.frame(pred_pf$predictions)
colnames(output_df) <- paste0("class", 0:2)
save_csv(output_df, "probability_output")

cat("  Class 0 prob mean:", mean(output_df$class0), "\n")

# ============================================================
# 6. SURVIVAL FOREST
# ============================================================
cat("\n=== Survival Forest ===\n")

n <- 500
p <- 5
X <- matrix(rnorm(n * p), n, p)
colnames(X) <- paste0("x", 1:p)
T_true <- rexp(n, rate = exp(0.5 * X[, 1]))
C <- rexp(n, rate = 0.5)
Y_time <- pmin(T_true, C)
D <- as.numeric(T_true <= C)  # 1 = event, 0 = censored

input_df <- data.frame(X, time = Y_time, status = D)
save_csv(input_df, "survival_input")

t0 <- proc.time()
sf <- survival_forest(X, Y_time, D,
  num.trees = 2000,
  seed = 42,
  honesty = TRUE,
  min.node.size = 5
)
t1 <- proc.time()
cat("  R training time:", (t1 - t0)[3], "seconds\n")

pred_sf <- predict(sf)
# Survival predictions are a matrix: rows=obs, cols=failure times
failure_times <- sf$failure.times
output_mat <- pred_sf$predictions

# Save first few failure time columns and the failure times
ft_df <- data.frame(failure_time = failure_times)
save_csv(ft_df, "survival_failure_times")

# Save predictions for first 20 failure times (or all if fewer)
ncols_save <- min(ncol(output_mat), 20)
output_df <- as.data.frame(output_mat[, 1:ncols_save])
colnames(output_df) <- paste0("t", 1:ncols_save)
save_csv(output_df, "survival_output")

cat("  Number of unique failure times:", length(failure_times), "\n")
cat("  Prediction matrix:", nrow(output_mat), "x", ncol(output_mat), "\n")

# ============================================================
# 7. CAUSAL FOREST: ATE TARGET SAMPLES
# ============================================================
cat("\n=== Causal ATE Target Samples ===\n")

# Use the causal forest (cf) fitted in section 2
# Reuse the same X, Y, W, cf objects from section 2

# ATE (treated)
ate_treated <- average_treatment_effect(cf, target.sample = "treated")
ate_treated_df <- data.frame(estimate = ate_treated[1], std_err = ate_treated[2])
save_csv(ate_treated_df, "causal_ate_treated")
cat("  ATE (treated):", ate_treated[1], "(SE:", ate_treated[2], ")\n")

# ATE (control)
ate_control <- average_treatment_effect(cf, target.sample = "control")
ate_control_df <- data.frame(estimate = ate_control[1], std_err = ate_control[2])
save_csv(ate_control_df, "causal_ate_control")
cat("  ATE (control):", ate_control[1], "(SE:", ate_control[2], ")\n")

# ATE (overlap)
ate_overlap <- average_treatment_effect(cf, target.sample = "overlap")
ate_overlap_df <- data.frame(estimate = ate_overlap[1], std_err = ate_overlap[2])
save_csv(ate_overlap_df, "causal_ate_overlap")
cat("  ATE (overlap):", ate_overlap[1], "(SE:", ate_overlap[2], ")\n")

# ============================================================
# 8. CAUSAL FOREST: AVERAGE PARTIAL EFFECT (continuous W)
# ============================================================
cat("\n=== Average Partial Effect ===\n")

# APE requires continuous treatment — create a new forest
n_ape <- 500
set.seed(42)
X_ape <- matrix(rnorm(n_ape * p), n_ape, p)
colnames(X_ape) <- paste0("x", 1:p)
W_cont <- rnorm(n_ape)
Y_cont <- 2 * X_ape[, 1] + X_ape[, 1] * W_cont + 0.5 * rnorm(n_ape)

# Save APE input data
ape_input_df <- data.frame(X_ape, y = Y_cont, w = W_cont)
save_csv(ape_input_df, "causal_ape_input")

cf_cont <- causal_forest(X_ape, Y_cont, W_cont,
  num.trees = 2000,
  seed = 42,
  honesty = TRUE,
  min.node.size = 5
)

ape <- average_treatment_effect(cf_cont)
ape_df <- data.frame(estimate = ape[1], std_err = ape[2])
save_csv(ape_df, "causal_ape")
cat("  APE:", ape[1], "(SE:", ape[2], ")\n")

# ============================================================
# 9. CAUSAL FOREST: BLP WITH EXPLICIT VCOV TYPES
# ============================================================
cat("\n=== BLP with HC0 and HC3 ===\n")

# BLP with HC0
# Use X_causal (saved from section 2) to avoid using overwritten X from later sections
blp_hc0 <- best_linear_projection(cf, X_causal, vcov.type = "HC0")
blp_hc0_df <- data.frame(
  variable = c("intercept", paste0("x", 1:p)),
  estimate = blp_hc0[, 1],
  std_error = blp_hc0[, 2],
  t_value = blp_hc0[, 3],
  p_value = blp_hc0[, 4]
)
save_csv(blp_hc0_df, "causal_blp_hc0")
cat("  BLP HC0 intercept:", blp_hc0[1, 1], "\n")

# BLP with HC3
# Use X_causal (saved from section 2) to avoid using overwritten X from later sections
blp_hc3 <- best_linear_projection(cf, X_causal, vcov.type = "HC3")
blp_hc3_df <- data.frame(
  variable = c("intercept", paste0("x", 1:p)),
  estimate = blp_hc3[, 1],
  std_error = blp_hc3[, 2],
  t_value = blp_hc3[, 3],
  p_value = blp_hc3[, 4]
)
save_csv(blp_hc3_df, "causal_blp_hc3")
cat("  BLP HC3 intercept:", blp_hc3[1, 1], "\n")

# ============================================================
# 10. CAUSAL FOREST: DR SCORES
# ============================================================
cat("\n=== DR Scores ===\n")

scores <- get_scores(cf)
scores_df <- data.frame(score = as.vector(scores))
save_csv(scores_df, "causal_scores")
cat("  Mean DR score:", mean(scores), "\n")
cat("  SD DR score:", sd(scores), "\n")

# ============================================================
# 11. SURVIVAL FOREST: EXPECTED SURVIVAL TIME
# ============================================================
cat("\n=== Expected Survival Time ===\n")

# Compute E[T|X] via trapezoidal integration of survival curves
pred_sf_full <- predict(sf)
surv_times <- sf$failure.times
surv_preds <- pred_sf_full$predictions

# Trapezoidal integration: E[T|X] = integral_0^t_max S(t|X) dt
# Assume S(0) = 1
E_T <- apply(surv_preds, 1, function(s) {
  n_t <- length(surv_times)
  # First interval [0, t_1]: 0.5 * (1 + S(t_1)) * t_1
  area <- 0.5 * (1 + s[1]) * surv_times[1]
  # Subsequent intervals
  if (n_t > 1) {
    for (j in 2:n_t) {
      area <- area + 0.5 * (s[j - 1] + s[j]) * (surv_times[j] - surv_times[j - 1])
    }
  }
  area
})

expected_surv_df <- data.frame(expected_time = E_T)
save_csv(expected_surv_df, "survival_expected")
cat("  Mean E[T|X]:", mean(E_T), "\n")
cat("  SD E[T|X]:", sd(E_T), "\n")

# ============================================================
# 12. RATE (Rank-Average Treatment Effect) — AUTOC
# ============================================================
cat("\n=== RATE AUTOC ===\n")

# Use the causal forest (cf) and its DR scores from section 2
scores_cf <- get_scores(cf)
pred_cf_cate <- predict(cf)$predictions

rate_out <- rank_average_treatment_effect(cf, pred_cf_cate)
rate_df <- data.frame(
  estimate = rate_out$estimate,
  std_err  = rate_out$std.err
)
save_csv(rate_df, "causal_rate_autoc")
cat("  AUTOC:", rate_out$estimate, "(SE:", rate_out$std.err, ")\n")

# ============================================================
# 13. LOCAL LINEAR REGRESSION FOREST
# ============================================================
cat("\n=== Local Linear Regression Forest ===\n")

set.seed(42)
n_ll <- 500
p_ll <- 5
X_ll <- matrix(rnorm(n_ll * p_ll), n_ll, p_ll)
colnames(X_ll) <- paste0("x", 1:p_ll)
Y_ll <- 2 * X_ll[, 1] + X_ll[, 2]^2 + 0.5 * rnorm(n_ll)

input_ll <- data.frame(X_ll, y = Y_ll)
save_csv(input_ll, "ll_regression_input")

t0 <- proc.time()
llf <- ll_regression_forest(X_ll, Y_ll,
  num.trees = 2000,
  seed = 42,
  honesty = TRUE,
  min.node.size = 5
)
t1 <- proc.time()
cat("  R training time:", (t1 - t0)[3], "seconds\n")

pred_ll <- predict(llf)
output_ll <- data.frame(prediction = pred_ll$predictions)
save_csv(output_ll, "ll_regression_output")
cat("  Prediction mean:", mean(pred_ll$predictions), "\n")

# ============================================================
# 14. LM FOREST (Linear Model Forest)
# ============================================================
cat("\n=== LM Forest ===\n")

# Use W separate from X to avoid degenerate nuisance (E[W|X]=W when W=X).
# DGP: Y = (1 + x1) * w + x1 + noise, so the heterogeneous coefficient of w
# is h(x) = 1 + x1.  Both R and Stata use regression_forest for W.hat,
# ensuring comparable nuisance estimates and high OOB prediction correlation.
set.seed(42)
n_lm <- 500
p_lm <- 5
X_lm <- matrix(rnorm(n_lm * p_lm), n_lm, p_lm)
colnames(X_lm) <- paste0("x", 1:p_lm)
W_lm <- rnorm(n_lm)                             # single treatment, not in X
tau_lm <- 1 + X_lm[, 1]                         # true heterogeneous coef
Y_lm <- tau_lm * W_lm + X_lm[, 1] + 0.5 * rnorm(n_lm)

input_lm <- data.frame(X_lm, w = W_lm, y = Y_lm)
save_csv(input_lm, "lm_forest_input")

t0 <- proc.time()
lmf <- lm_forest(X_lm, Y_lm, matrix(W_lm, ncol = 1),
  num.trees = 2000,
  seed = 42,
  honesty = TRUE,
  min.node.size = 5
)
t1 <- proc.time()
cat("  R training time:", (t1 - t0)[3], "seconds\n")

pred_lm <- predict(lmf)
# lm_forest returns one coefficient column per W column; here ncol(W)=1
output_lm <- as.data.frame(pred_lm$predictions)
colnames(output_lm) <- paste0("coef", seq_len(ncol(pred_lm$predictions)))
save_csv(output_lm, "lm_forest_output")
cat("  Number of coefficient columns:", ncol(pred_lm$predictions), "\n")

# ============================================================
# 15. MULTI-ARM CAUSAL FOREST
# ============================================================
cat("\n=== Multi-Arm Causal Forest ===\n")

set.seed(42)
n_ma <- 500
p_ma <- 5
X_ma <- matrix(rnorm(n_ma * p_ma), n_ma, p_ma)
colnames(X_ma) <- paste0("x", 1:p_ma)
W_ma <- factor(sample(c("control", "treat1", "treat2"), n_ma, replace = TRUE))
tau1 <- X_ma[, 1]       # effect of treat1 vs control
tau2 <- 0.5 * X_ma[, 2] # effect of treat2 vs control
Y_ma <- X_ma[, 1] +
  ifelse(W_ma == "treat1", tau1, 0) +
  ifelse(W_ma == "treat2", tau2, 0) +
  0.5 * rnorm(n_ma)

input_ma <- data.frame(X_ma, y = Y_ma, w = as.numeric(W_ma) - 1)
save_csv(input_ma, "multi_arm_input")

t0 <- proc.time()
maf <- multi_arm_causal_forest(X_ma, Y_ma, W_ma,
  num.trees = 2000,
  seed = 42,
  honesty = TRUE,
  min.node.size = 5
)
t1 <- proc.time()
cat("  R training time:", (t1 - t0)[3], "seconds\n")

pred_ma <- predict(maf)
# predictions is a 3D-ish structure: n x (num_arms - 1) contrasts
output_ma <- as.data.frame(pred_ma$predictions)
colnames(output_ma) <- paste0("contrast", seq_len(ncol(output_ma)))
save_csv(output_ma, "multi_arm_output")
cat("  Number of contrast columns:", ncol(output_ma), "\n")

# ============================================================
# 16. MULTI-REGRESSION FOREST
# ============================================================
cat("\n=== Multi-Regression Forest ===\n")

set.seed(42)
n_mr <- 500
p_mr <- 5
X_mr <- matrix(rnorm(n_mr * p_mr), n_mr, p_mr)
colnames(X_mr) <- paste0("x", 1:p_mr)
Y_mr <- cbind(
  2 * X_mr[, 1] + 0.5 * rnorm(n_mr),
  X_mr[, 2]^2 + 0.5 * rnorm(n_mr)
)
colnames(Y_mr) <- c("y1", "y2")

input_mr <- data.frame(X_mr, Y_mr)
save_csv(input_mr, "multi_regression_input")

t0 <- proc.time()
mrf <- multi_regression_forest(X_mr, Y_mr,
  num.trees = 2000,
  seed = 42,
  honesty = TRUE,
  min.node.size = 5
)
t1 <- proc.time()
cat("  R training time:", (t1 - t0)[3], "seconds\n")

pred_mr <- predict(mrf)
output_mr <- as.data.frame(pred_mr$predictions)
colnames(output_mr) <- paste0("pred_y", seq_len(ncol(output_mr)))
save_csv(output_mr, "multi_regression_output")
cat("  Number of outcome columns:", ncol(output_mr), "\n")

# ============================================================
# 17. CAUSAL SURVIVAL FOREST
# ============================================================
cat("\n=== Causal Survival Forest ===\n")

set.seed(42)
n_cs <- 500
p_cs <- 5
X_cs <- matrix(rnorm(n_cs * p_cs), n_cs, p_cs)
colnames(X_cs) <- paste0("x", 1:p_cs)
W_cs <- rbinom(n_cs, 1, 0.5)
T_true_cs <- rexp(n_cs, rate = exp(0.3 * X_cs[, 1] + 0.2 * W_cs))
C_cs <- rexp(n_cs, rate = 0.3)
Y_cs <- pmin(T_true_cs, C_cs)
D_cs <- as.numeric(T_true_cs <= C_cs)

input_cs <- data.frame(X_cs, time = Y_cs, status = D_cs, w = W_cs)
save_csv(input_cs, "causal_survival_input")

# causal_survival_forest requires a horizon
horizon_val <- median(Y_cs[D_cs == 1])

t0 <- proc.time()
csf <- causal_survival_forest(X_cs, Y_cs, W_cs, D_cs,
  horizon = horizon_val,
  num.trees = 2000,
  seed = 42,
  honesty = TRUE,
  min.node.size = 15
)
t1 <- proc.time()
cat("  R training time:", (t1 - t0)[3], "seconds\n")

pred_cs <- predict(csf)
output_cs <- data.frame(cate = pred_cs$predictions)
save_csv(output_cs, "causal_survival_output")
cat("  CATE mean:", mean(pred_cs$predictions), "\n")
cat("  Horizon used:", horizon_val, "\n")

# Save horizon for Stata test to use
horizon_df <- data.frame(horizon = horizon_val)
save_csv(horizon_df, "causal_survival_horizon")

# ============================================================
# 18. BOOSTED REGRESSION FOREST
# ============================================================
cat("\n=== Boosted Regression Forest ===\n")

set.seed(42)
n_br <- 500
p_br <- 5
X_br <- matrix(rnorm(n_br * p_br), n_br, p_br)
colnames(X_br) <- paste0("x", 1:p_br)
Y_br <- 2 * X_br[, 1] + X_br[, 2]^2 + sin(X_br[, 3]) + 0.5 * rnorm(n_br)

input_br <- data.frame(X_br, y = Y_br)
save_csv(input_br, "boosted_regression_input")

t0 <- proc.time()
brf <- boosted_regression_forest(X_br, Y_br,
  num.trees = 2000,
  seed = 42,
  honesty = TRUE,
  min.node.size = 5
)
t1 <- proc.time()
cat("  R training time:", (t1 - t0)[3], "seconds\n")

pred_br <- predict(brf)
output_br <- data.frame(prediction = pred_br$predictions)
save_csv(output_br, "boosted_regression_output")
cat("  Prediction mean:", mean(pred_br$predictions), "\n")
cat("  Number of boosts:", brf$num.boosts, "\n")

# ============================================================
# 19. TMLE ATE/ATT/ATC
# ============================================================
cat("\n=== TMLE ===\n")

# Use the causal forest from section 2
ate_tmle <- average_treatment_effect(cf, method = "TMLE")
ate_tmle_df <- data.frame(estimate = ate_tmle[1], std_err = ate_tmle[2])
save_csv(ate_tmle_df, "causal_ate_tmle")
cat("  ATE (TMLE):", ate_tmle[1], "(SE:", ate_tmle[2], ")\n")

att_tmle <- average_treatment_effect(cf, method = "TMLE", target.sample = "treated")
att_tmle_df <- data.frame(estimate = att_tmle[1], std_err = att_tmle[2])
save_csv(att_tmle_df, "causal_att_tmle")
cat("  ATT (TMLE):", att_tmle[1], "(SE:", att_tmle[2], ")\n")

atc_tmle <- average_treatment_effect(cf, method = "TMLE", target.sample = "control")
atc_tmle_df <- data.frame(estimate = atc_tmle[1], std_err = atc_tmle[2])
save_csv(atc_tmle_df, "causal_atc_tmle")
cat("  ATC (TMLE):", atc_tmle[1], "(SE:", atc_tmle[2], ")\n")

# ============================================================
# 20. CLUSTERED INFERENCE
# ============================================================
cat("\n=== Clustered Inference ===\n")

set.seed(42)
n_cl <- 500
p_cl <- 5
X_cl <- matrix(rnorm(n_cl * p_cl), n_cl, p_cl)
colnames(X_cl) <- paste0("x", 1:p_cl)
cl <- rep(1:50, each = 10)  # 50 clusters of size 10
W_cl <- rbinom(n_cl, 1, 0.5)
tau_cl <- X_cl[, 1]
Y_cl <- 2 * X_cl[, 1] + tau_cl * W_cl + 0.5 * rnorm(n_cl)

input_cl <- data.frame(X_cl, y = Y_cl, w = W_cl, cluster = cl)
save_csv(input_cl, "causal_clustered_input")

cf_cl <- causal_forest(X_cl, Y_cl, W_cl, clusters = cl,
  num.trees = 2000, seed = 42, honesty = TRUE, min.node.size = 5)

ate_cl <- average_treatment_effect(cf_cl)
ate_cl_df <- data.frame(estimate = ate_cl[1], std_err = ate_cl[2])
save_csv(ate_cl_df, "causal_clustered_ate")
cat("  Clustered ATE:", ate_cl[1], "(SE:", ate_cl[2], ")\n")

blp_cl <- best_linear_projection(cf_cl, X_cl[, 1:2])
blp_cl_df <- data.frame(
  variable = c("intercept", "x1", "x2"),
  estimate = blp_cl[, 1],
  std_error = blp_cl[, 2],
  t_value = blp_cl[, 3],
  p_value = blp_cl[, 4]
)
save_csv(blp_cl_df, "causal_clustered_blp")
cat("  Clustered BLP intercept:", blp_cl[1, 1], "\n")

pred_cl <- predict(cf_cl)
output_cl <- data.frame(cate = pred_cl$predictions)
save_csv(output_cl, "causal_clustered_output")

# ============================================================
# 21. DEBIASING WEIGHTS
# ============================================================
cat("\n=== Debiasing Weights ===\n")

# Use the causal forest from section 2 with custom debiasing weights
dbw <- abs(X[, 1]) + 0.5
ate_dbw <- average_treatment_effect(cf, debiasing.weights = dbw)
ate_dbw_df <- data.frame(estimate = ate_dbw[1], std_err = ate_dbw[2])
save_csv(ate_dbw_df, "causal_ate_debiasing")
cat("  ATE (debiasing weights):", ate_dbw[1], "(SE:", ate_dbw[2], ")\n")

# Also save the weights for Stata to use
dbw_df <- data.frame(debiasing_weight = dbw)
save_csv(dbw_df, "causal_debiasing_weights")

# ============================================================
# 22. EQUALIZE CLUSTER WEIGHTS
# ============================================================
cat("\n=== Equalize Cluster Weights ===\n")

cf_eq <- causal_forest(X_cl, Y_cl, W_cl, clusters = cl,
  equalize.cluster.weights = TRUE,
  num.trees = 2000, seed = 42, honesty = TRUE, min.node.size = 5)

ate_eq <- average_treatment_effect(cf_eq)
ate_eq_df <- data.frame(estimate = ate_eq[1], std_err = ate_eq[2])
save_csv(ate_eq_df, "causal_equalize_cluster_ate")
cat("  Equalized ATE:", ate_eq[1], "(SE:", ate_eq[2], ")\n")

pred_eq <- predict(cf_eq)
output_eq <- data.frame(cate = pred_eq$predictions)
save_csv(output_eq, "causal_equalize_cluster_output")

# ============================================================
# 23. USER-SUPPLIED NUISANCE
# ============================================================
cat("\n=== User-Supplied Nuisance ===\n")

# Fit nuisance models externally, then pass to causal_forest
# Using the same data from section 2
rf_y <- regression_forest(X, Y, num.trees = 500, seed = 42)
Y_hat <- predict(rf_y)$predictions
rf_w <- regression_forest(X, W, num.trees = 500, seed = 42)
W_hat <- predict(rf_w)$predictions

cf_nuis <- causal_forest(X, Y, W,
  Y.hat = Y_hat, W.hat = W_hat,
  num.trees = 2000, seed = 42, honesty = TRUE, min.node.size = 5)

ate_nuis <- average_treatment_effect(cf_nuis)
ate_nuis_df <- data.frame(estimate = ate_nuis[1], std_err = ate_nuis[2])
save_csv(ate_nuis_df, "causal_nuisance_ate")
cat("  Nuisance ATE:", ate_nuis[1], "(SE:", ate_nuis[2], ")\n")

nuisance_df <- data.frame(y_hat = Y_hat, w_hat = W_hat)
save_csv(nuisance_df, "causal_nuisance_inputs")

pred_nuis <- predict(cf_nuis)
output_nuis <- data.frame(cate = pred_nuis$predictions)
save_csv(output_nuis, "causal_nuisance_output")

# ============================================================
# 24. QUANTILE FOREST WITH REGRESSION SPLITTING
# ============================================================
cat("\n=== Quantile Forest + Regression Splitting ===\n")

set.seed(42)
X_qrs <- matrix(rnorm(500 * 5), 500, 5)
colnames(X_qrs) <- paste0("x", 1:5)
Y_qrs <- 2 * X_qrs[, 1] + X_qrs[, 2] + rnorm(500) * (1 + abs(X_qrs[, 3]))

qf_rs <- quantile_forest(X_qrs, Y_qrs,
  regression.splitting = TRUE,
  quantiles = c(0.1, 0.5, 0.9),
  num.trees = 2000, seed = 42, honesty = TRUE, min.node.size = 5)

pred_qrs <- predict(qf_rs, quantiles = c(0.1, 0.5, 0.9))
output_qrs <- as.data.frame(pred_qrs$predictions)
colnames(output_qrs) <- c("q0_1", "q0_5", "q0_9")

input_qrs <- data.frame(X_qrs, y = Y_qrs)
save_csv(input_qrs, "quantile_regsplit_input")
save_csv(output_qrs, "quantile_regsplit_output")
cat("  Median prediction mean:", mean(output_qrs$q0_5), "\n")

# ============================================================
# 25. VARIABLE IMPORTANCE WITH DIFFERENT DECAY EXPONENT
# ============================================================
cat("\n=== Variable Importance (decay exponent) ===\n")

# Use the regression forest from section 1
vi_d1 <- variable_importance(rf, decay.exponent = 1.0)
vi_d1_df <- data.frame(variable = paste0("x", 1:p), importance = as.vector(vi_d1))
save_csv(vi_d1_df, "regression_vi_decay1")
cat("  Decay=1.0:", as.vector(vi_d1), "\n")

vi_d3 <- variable_importance(rf, decay.exponent = 3.0)
vi_d3_df <- data.frame(variable = paste0("x", 1:p), importance = as.vector(vi_d3))
save_csv(vi_d3_df, "regression_vi_decay3")
cat("  Decay=3.0:", as.vector(vi_d3), "\n")

# ============================================================
# Summary
# ============================================================
cat("\n=== Reference Data Summary ===\n")
files <- list.files(outdir, pattern = "\\.csv$")
cat("Generated", length(files), "reference files in", outdir, "/\n")
for (f in files) {
  cat("  ", f, "\n")
}
cat("\nDone.\n")
