#!/usr/bin/env Rscript
# R Reference Script for BLP, Test Calibration, Get Scores
# Generates all reference results for R vs Stata fidelity comparison

library(grf)
library(jsonlite)

cat("grf version:", as.character(packageVersion("grf")), "\n")

# Check what target.sample values are valid in this grf version
cat("BLP target.sample options (checking grf source):\n")
tryCatch({
  formals_blp <- formals(best_linear_projection)
  cat("target.sample default:", deparse(formals_blp$target.sample), "\n")
}, error = function(e) cat("Could not inspect formals\n"))

# Helper to extract from coeftest object
extract_blp <- function(blp) {
  list(
    coefs      = as.vector(blp[, 1]),
    ses        = as.vector(blp[, 2]),
    tvals      = as.vector(blp[, 3]),
    pvals      = as.vector(blp[, 4]),
    coef_names = rownames(blp)
  )
}

# ============================================================
# DGP (fixed seed)
# ============================================================
set.seed(42)
n <- 1000
p <- 5
X <- matrix(rnorm(n * p), n, p)
W <- rbinom(n, 1, 0.5)
tau <- X[, 1] + X[, 2]
Y <- X[, 1] + tau * W + rnorm(n)

# Cluster variable (100 clusters of 10)
clusters <- rep(1:100, each = 10)

# Custom debiasing weights (random positive weights)
set.seed(123)
dbw <- abs(rnorm(n)) + 0.5

# ============================================================
# Fit causal forest
# ============================================================
cat("\nFitting causal forest...\n")
cf <- causal_forest(X, Y, W,
  num.trees = 2000,
  seed = 42
)
cat("Forest fitted.\n")

# ============================================================
# Save data and forest predictions for Stata
# ============================================================
data_df <- data.frame(
  Y = Y, W = W,
  X1 = X[, 1], X2 = X[, 2], X3 = X[, 3],
  X4 = X[, 4], X5 = X[, 5],
  tau_true = tau,
  clusters = clusters,
  dbw = dbw
)
write.csv(data_df,
  "/tmp/grf_stata/tests/fidelity_reports/14_blp_calibration/data.csv",
  row.names = FALSE)
cat("Data saved.\n")

predictions_df <- data.frame(
  tau_hat = cf$predictions,
  yhat    = cf$Y.hat,
  what    = cf$W.hat
)
write.csv(predictions_df,
  "/tmp/grf_stata/tests/fidelity_reports/14_blp_calibration/predictions.csv",
  row.names = FALSE)
cat("Predictions saved.\n")

# ============================================================
# BLP tests
# ============================================================
results <- list()

cat("\n=== BLP Tests ===\n")

# Test 1: BLP default HC3, all covariates
cat("Test 1: BLP HC3 (default), all covariates\n")
blp1 <- best_linear_projection(cf, X)
print(blp1)
results$blp_hc3_all <- extract_blp(blp1)

# Test 2: BLP HC0
cat("Test 2: BLP HC0\n")
blp2 <- best_linear_projection(cf, X, vcov.type = "HC0")
print(blp2)
results$blp_hc0 <- extract_blp(blp2)

# Test 3: BLP HC1
cat("Test 3: BLP HC1\n")
blp3 <- best_linear_projection(cf, X, vcov.type = "HC1")
print(blp3)
results$blp_hc1 <- extract_blp(blp3)

# Test 4: BLP HC2
cat("Test 4: BLP HC2\n")
blp4 <- best_linear_projection(cf, X, vcov.type = "HC2")
print(blp4)
results$blp_hc2 <- extract_blp(blp4)

# Test 5: BLP subset (X1, X2 only)
cat("Test 5: BLP subset X1, X2\n")
blp5 <- best_linear_projection(cf, X[, 1:2])
print(blp5)
results$blp_subset_x1x2 <- extract_blp(blp5)

# Test 6: BLP target.sample=treated (check if supported)
cat("Test 6: BLP target.sample=treated\n")
blp6_result <- tryCatch({
  blp6 <- best_linear_projection(cf, X, target.sample = "treated")
  print(blp6)
  c(extract_blp(blp6), list(supported = TRUE))
}, error = function(e) {
  cat("  Not supported in grf", as.character(packageVersion("grf")), ":", conditionMessage(e), "\n")
  list(supported = FALSE, error = conditionMessage(e))
})
results$blp_treated <- blp6_result

# Test 7: BLP target.sample=control (check if supported)
cat("Test 7: BLP target.sample=control\n")
blp7_result <- tryCatch({
  blp7 <- best_linear_projection(cf, X, target.sample = "control")
  print(blp7)
  c(extract_blp(blp7), list(supported = TRUE))
}, error = function(e) {
  cat("  Not supported in grf", as.character(packageVersion("grf")), ":", conditionMessage(e), "\n")
  list(supported = FALSE, error = conditionMessage(e))
})
results$blp_control <- blp7_result

# Test 8: BLP target.sample=overlap
cat("Test 8: BLP target.sample=overlap\n")
blp8 <- best_linear_projection(cf, X, target.sample = "overlap")
print(blp8)
results$blp_overlap <- extract_blp(blp8)

# Test 9: BLP with debiasing weights
cat("Test 9: BLP with debiasing weights\n")
blp9 <- best_linear_projection(cf, X, debiasing.weights = dbw)
print(blp9)
results$blp_dbw <- extract_blp(blp9)

# Test 10: BLP with clusters
cat("Test 10: BLP with clusters\n")
cf_cl <- causal_forest(X, Y, W,
  num.trees = 2000,
  clusters  = clusters,
  seed      = 42
)
blp10 <- best_linear_projection(cf_cl, X)
print(blp10)
results$blp_clusters <- extract_blp(blp10)

# Save clustered forest predictions
pred_cl <- data.frame(
  tau_hat = cf_cl$predictions,
  yhat    = cf_cl$Y.hat,
  what    = cf_cl$W.hat
)
write.csv(pred_cl,
  "/tmp/grf_stata/tests/fidelity_reports/14_blp_calibration/predictions_cl.csv",
  row.names = FALSE)

# ============================================================
# Test Calibration tests
# ============================================================
cat("\n=== Test Calibration ===\n")

# Test 13: Basic calibration
cat("Test 13: Basic calibration\n")
cal1 <- test_calibration(cf)
print(cal1)
results$calibration_basic <- list(
  b_mean  = unname(cal1["mean.forest.prediction", "Estimate"]),
  se_mean = unname(cal1["mean.forest.prediction", "Std. Error"]),
  t_mean  = unname(cal1["mean.forest.prediction", "t value"]),
  p_mean  = unname(cal1["mean.forest.prediction", "Pr(>t)"]),
  b_diff  = unname(cal1["differential.forest.prediction", "Estimate"]),
  se_diff = unname(cal1["differential.forest.prediction", "Std. Error"]),
  t_diff  = unname(cal1["differential.forest.prediction", "t value"]),
  p_diff  = unname(cal1["differential.forest.prediction", "Pr(>t)"])
)

# Test 14: Strong heterogeneity
cat("Test 14: Calibration with strong heterogeneity\n")
set.seed(42)
tau_strong <- 3 * X[, 1] + 3 * X[, 2]
Y_strong   <- X[, 1] + tau_strong * W + rnorm(n)
cf_strong  <- causal_forest(X, Y_strong, W, num.trees = 2000, seed = 42)
# Save data_strong for Stata
data_strong_df <- cbind(data_df[, c("W","X1","X2","X3","X4","X5","clusters","dbw")],
  Y_strong = Y_strong)
write.csv(data_strong_df,
  "/tmp/grf_stata/tests/fidelity_reports/14_blp_calibration/data_strong.csv",
  row.names = FALSE)
cal2 <- test_calibration(cf_strong)
print(cal2)
results$calibration_strong_het <- list(
  b_mean  = unname(cal2["mean.forest.prediction", "Estimate"]),
  se_mean = unname(cal2["mean.forest.prediction", "Std. Error"]),
  t_mean  = unname(cal2["mean.forest.prediction", "t value"]),
  p_mean  = unname(cal2["mean.forest.prediction", "Pr(>t)"]),
  b_diff  = unname(cal2["differential.forest.prediction", "Estimate"]),
  se_diff = unname(cal2["differential.forest.prediction", "Std. Error"]),
  t_diff  = unname(cal2["differential.forest.prediction", "t value"]),
  p_diff  = unname(cal2["differential.forest.prediction", "Pr(>t)"])
)
pred_strong <- data.frame(
  tau_hat = cf_strong$predictions,
  yhat    = cf_strong$Y.hat,
  what    = cf_strong$W.hat
)
write.csv(pred_strong,
  "/tmp/grf_stata/tests/fidelity_reports/14_blp_calibration/predictions_strong.csv",
  row.names = FALSE)

# Test 15: No heterogeneity (constant tau)
cat("Test 15: Calibration with no heterogeneity\n")
set.seed(42)
tau_const <- rep(1.5, n)
Y_const   <- X[, 1] + tau_const * W + rnorm(n)
cf_const  <- causal_forest(X, Y_const, W, num.trees = 2000, seed = 42)
# Save data_const for Stata
data_const_df <- cbind(data_df[, c("W","X1","X2","X3","X4","X5","clusters","dbw")],
  Y_const = Y_const)
write.csv(data_const_df,
  "/tmp/grf_stata/tests/fidelity_reports/14_blp_calibration/data_const.csv",
  row.names = FALSE)
cal3 <- test_calibration(cf_const)
print(cal3)
results$calibration_no_het <- list(
  b_mean  = unname(cal3["mean.forest.prediction", "Estimate"]),
  se_mean = unname(cal3["mean.forest.prediction", "Std. Error"]),
  t_mean  = unname(cal3["mean.forest.prediction", "t value"]),
  p_mean  = unname(cal3["mean.forest.prediction", "Pr(>t)"]),
  b_diff  = unname(cal3["differential.forest.prediction", "Estimate"]),
  se_diff = unname(cal3["differential.forest.prediction", "Std. Error"]),
  t_diff  = unname(cal3["differential.forest.prediction", "t value"]),
  p_diff  = unname(cal3["differential.forest.prediction", "Pr(>t)"])
)
pred_const <- data.frame(
  tau_hat = cf_const$predictions,
  yhat    = cf_const$Y.hat,
  what    = cf_const$W.hat
)
write.csv(pred_const,
  "/tmp/grf_stata/tests/fidelity_reports/14_blp_calibration/predictions_const.csv",
  row.names = FALSE)

# ============================================================
# Get Scores tests
# ============================================================
cat("\n=== Get Scores ===\n")

# Test 17-20: DR scores
cat("Test 17: DR scores from causal forest\n")
dr_scores <- get_scores(cf)
cat("Mean DR scores:", mean(dr_scores), "\n")
cat("SD DR scores:", sd(dr_scores), "\n")
cat("Min DR scores:", min(dr_scores), "\n")
cat("Max DR scores:", max(dr_scores), "\n")
results$dr_scores_summary <- list(
  mean = mean(dr_scores),
  sd   = sd(dr_scores),
  min  = min(dr_scores),
  max  = max(dr_scores),
  n    = length(dr_scores)
)

scores_df <- data.frame(dr_score_r = as.vector(dr_scores))
write.csv(scores_df,
  "/tmp/grf_stata/tests/fidelity_reports/14_blp_calibration/dr_scores.csv",
  row.names = FALSE)
cat("DR scores saved.\n")

# ============================================================
# Save all results to JSON
# ============================================================
write_json(results,
  "/tmp/grf_stata/tests/fidelity_reports/14_blp_calibration/r_results.json",
  auto_unbox = TRUE, digits = 12)

cat("\nR reference script completed successfully.\n")
