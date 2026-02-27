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
ate_df <- data.frame(estimate = ate[1], std.err = ate[2])
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
  std.error = blp[, 2],
  t.value = blp[, 3],
  p.value = blp[, 4]
)
save_csv(blp_df, "causal_blp")

# Test calibration
tc <- test_calibration(cf)
tc_df <- data.frame(
  term = rownames(tc),
  estimate = tc[, 1],
  std.error = tc[, 2],
  t.value = tc[, 3],
  p.value = tc[, 4]
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
colnames(output_df) <- paste0("q", quantiles)
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
ate_treated_df <- data.frame(estimate = ate_treated[1], std.err = ate_treated[2])
save_csv(ate_treated_df, "causal_ate_treated")
cat("  ATE (treated):", ate_treated[1], "(SE:", ate_treated[2], ")\n")

# ATE (control)
ate_control <- average_treatment_effect(cf, target.sample = "control")
ate_control_df <- data.frame(estimate = ate_control[1], std.err = ate_control[2])
save_csv(ate_control_df, "causal_ate_control")
cat("  ATE (control):", ate_control[1], "(SE:", ate_control[2], ")\n")

# ATE (overlap)
ate_overlap <- average_treatment_effect(cf, target.sample = "overlap")
ate_overlap_df <- data.frame(estimate = ate_overlap[1], std.err = ate_overlap[2])
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

ape <- average_partial_effect(cf_cont)
ape_df <- data.frame(estimate = ape[1], std.err = ape[2])
save_csv(ape_df, "causal_ape")
cat("  APE:", ape[1], "(SE:", ape[2], ")\n")

# ============================================================
# 9. CAUSAL FOREST: BLP WITH EXPLICIT VCOV TYPES
# ============================================================
cat("\n=== BLP with HC0 and HC3 ===\n")

# BLP with HC0
blp_hc0 <- best_linear_projection(cf, X, vcov.type = "HC0")
blp_hc0_df <- data.frame(
  variable = c("intercept", paste0("x", 1:p)),
  estimate = blp_hc0[, 1],
  std.error = blp_hc0[, 2],
  t.value = blp_hc0[, 3],
  p.value = blp_hc0[, 4]
)
save_csv(blp_hc0_df, "causal_blp_hc0")
cat("  BLP HC0 intercept:", blp_hc0[1, 1], "\n")

# BLP with HC3
blp_hc3 <- best_linear_projection(cf, X, vcov.type = "HC3")
blp_hc3_df <- data.frame(
  variable = c("intercept", paste0("x", 1:p)),
  estimate = blp_hc3[, 1],
  std.error = blp_hc3[, 2],
  t.value = blp_hc3[, 3],
  p.value = blp_hc3[, 4]
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
# Summary
# ============================================================
cat("\n=== Reference Data Summary ===\n")
files <- list.files(outdir, pattern = "\\.csv$")
cat("Generated", length(files), "reference files in", outdir, "/\n")
for (f in files) {
  cat("  ", f, "\n")
}
cat("\nDone.\n")
