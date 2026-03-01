#!/usr/bin/env Rscript
# Test 15: Variable Importance, RATE, and Tune Parameters
# R reference script

library(grf)
set.seed(42)

cat("=== GRF Test 15: VI, RATE, Tune ===\n")
cat("R version:", R.version$major, ".", R.version$minor, "\n", sep="")
cat("grf version:", as.character(packageVersion("grf")), "\n\n")

# ---- DGP ----
set.seed(42)
n <- 500
p <- 10

X <- matrix(rnorm(n * p), n, p)
colnames(X) <- paste0("x", 1:p)

# Regression DGP: only X1 and X2 matter
Y <- 3 * X[, 1] + 2 * X[, 2] + rnorm(n)

# Causal DGP
W <- rbinom(n, 1, 0.5)
tau <- X[, 1] + X[, 2]
Y_causal <- X[, 1] + tau * W + rnorm(n)

# Cluster variable (5 clusters)
cluster <- rep(1:5, each = 100)

# Weights
weights_vec <- runif(n, 0.5, 1.5)

# Held-out test set for tune MSE test
set.seed(123)
n_test <- 200
X_test <- matrix(rnorm(n_test * p), n_test, p)
Y_test <- 3 * X_test[, 1] + 2 * X_test[, 2] + rnorm(n_test)

# ---- Helper function to write results ----
write_result <- function(test_id, test_name, values, extra = "") {
  fname <- sprintf("/tmp/grf_stata/tests/fidelity_reports/15_vi_rate_tune/r_test%02d.txt", test_id)
  lines <- c(
    sprintf("TEST %d: %s", test_id, test_name),
    sprintf("STATUS: OK"),
    sprintf("N: %d", n),
    extra
  )
  if (is.numeric(values) && !is.null(names(values))) {
    for (nm in names(values)) {
      lines <- c(lines, sprintf("%s: %.8f", nm, values[nm]))
    }
  } else if (is.numeric(values)) {
    for (i in seq_along(values)) {
      lines <- c(lines, sprintf("val%d: %.8f", i, values[i]))
    }
  }
  writeLines(lines, fname)
  cat(sprintf("  [DONE] Test %d written to %s\n", test_id, fname))
}

# ---- Base regression forest ----
cat("Fitting base regression forest for VI tests...\n")
rf_base <- regression_forest(X, Y, num.trees = 500, seed = 42)

# ============================================================
# VARIABLE IMPORTANCE TESTS (1-11)
# ============================================================
cat("\n--- Variable Importance Tests ---\n")

# Test 1: Default (decay=2, depth=4)
cat("Test 1: Default VI (decay=2, depth=4)\n")
vi1 <- variable_importance(rf_base, decay.exponent = 2, max.depth = 4)
names(vi1) <- colnames(X)
write_result(1, "VI default decay=2 depth=4", vi1)

# Test 2: decay.exponent=1.0
cat("Test 2: VI decay=1.0\n")
vi2 <- variable_importance(rf_base, decay.exponent = 1.0, max.depth = 4)
names(vi2) <- colnames(X)
write_result(2, "VI decay=1.0", vi2)

# Test 3: decay.exponent=3.0
cat("Test 3: VI decay=3.0\n")
vi3 <- variable_importance(rf_base, decay.exponent = 3.0, max.depth = 4)
names(vi3) <- colnames(X)
write_result(3, "VI decay=3.0", vi3)

# Test 4: decay.exponent=0.5
cat("Test 4: VI decay=0.5\n")
vi4 <- variable_importance(rf_base, decay.exponent = 0.5, max.depth = 4)
names(vi4) <- colnames(X)
write_result(4, "VI decay=0.5", vi4)

# Test 5: max.depth=2
cat("Test 5: VI max.depth=2\n")
vi5 <- variable_importance(rf_base, decay.exponent = 2, max.depth = 2)
names(vi5) <- colnames(X)
write_result(5, "VI max.depth=2", vi5)

# Test 6: max.depth=8
cat("Test 6: VI max.depth=8\n")
vi6 <- variable_importance(rf_base, decay.exponent = 2, max.depth = 8)
names(vi6) <- colnames(X)
write_result(6, "VI max.depth=8", vi6)

# Test 7: With clusters
cat("Test 7: VI with clusters\n")
rf_cluster <- regression_forest(X, Y, num.trees = 500, seed = 42,
                                 clusters = cluster)
vi7 <- variable_importance(rf_cluster, decay.exponent = 2, max.depth = 4)
names(vi7) <- colnames(X)
write_result(7, "VI with cluster", vi7)

# Test 8: With weights
cat("Test 8: VI with weights\n")
rf_weighted <- regression_forest(X, Y, num.trees = 500, seed = 42,
                                  sample.weights = weights_vec)
vi8 <- variable_importance(rf_weighted, decay.exponent = 2, max.depth = 4)
names(vi8) <- colnames(X)
write_result(8, "VI with weights", vi8)

# Test 9: ntrees=1000
cat("Test 9: VI ntrees=1000\n")
rf_1000 <- regression_forest(X, Y, num.trees = 1000, seed = 42)
vi9 <- variable_importance(rf_1000, decay.exponent = 2, max.depth = 4)
names(vi9) <- colnames(X)
write_result(9, "VI ntrees=1000", vi9)

# Test 10: VI ranking — top-2 should be x1, x2
cat("Test 10: VI ranking (top-2 = x1, x2)\n")
vi10 <- variable_importance(rf_base, decay.exponent = 2, max.depth = 4)
names(vi10) <- colnames(X)
ranking <- order(vi10, decreasing = TRUE)
top2 <- colnames(X)[ranking[1:2]]
top2_correct <- all(c("x1", "x2") %in% top2)
cat(sprintf("  Top-2 variables: %s, %s\n", top2[1], top2[2]))
cat(sprintf("  x1 rank: %d, x2 rank: %d\n", which(ranking == 1), which(ranking == 2)))
lines10 <- c(
  sprintf("TEST 10: VI ranking top-2 should be x1 x2"),
  sprintf("STATUS: %s", ifelse(top2_correct, "PASS", "FAIL")),
  sprintf("top1: %s", top2[1]),
  sprintf("top2: %s", top2[2]),
  sprintf("x1_rank: %d", which(ranking == 1)),
  sprintf("x2_rank: %d", which(ranking == 2)),
  sprintf("top2_correct: %d", as.integer(top2_correct))
)
for (nm in names(vi10)) {
  lines10 <- c(lines10, sprintf("%s: %.8f", nm, vi10[nm]))
}
writeLines(lines10, "/tmp/grf_stata/tests/fidelity_reports/15_vi_rate_tune/r_test10.txt")
cat("  [DONE] Test 10 written\n")

# Test 11: Many irrelevant vars (p=20)
cat("Test 11: VI with p=20 (many irrelevant vars)\n")
set.seed(42)
p20 <- 20
X20 <- matrix(rnorm(n * p20), n, p20)
colnames(X20) <- paste0("x", 1:p20)
Y20 <- 3 * X20[, 1] + 2 * X20[, 2] + rnorm(n)
rf20 <- regression_forest(X20, Y20, num.trees = 500, seed = 42)
vi11 <- variable_importance(rf20, decay.exponent = 2, max.depth = 4)
names(vi11) <- colnames(X20)
ranking11 <- order(vi11, decreasing = TRUE)
top2_11 <- colnames(X20)[ranking11[1:2]]
top2_correct11 <- all(c("x1", "x2") %in% top2_11)
lines11 <- c(
  sprintf("TEST 11: VI p=20 top-2 should be x1 x2"),
  sprintf("STATUS: %s", ifelse(top2_correct11, "PASS", "FAIL")),
  sprintf("top1: %s", top2_11[1]),
  sprintf("top2: %s", top2_11[2]),
  sprintf("x1_rank: %d", which(ranking11 == 1)),
  sprintf("x2_rank: %d", which(ranking11 == 2)),
  sprintf("top2_correct: %d", as.integer(top2_correct11))
)
for (nm in names(vi11)) {
  lines11 <- c(lines11, sprintf("%s: %.8f", nm, vi11[nm]))
}
writeLines(lines11, "/tmp/grf_stata/tests/fidelity_reports/15_vi_rate_tune/r_test11.txt")
cat("  [DONE] Test 11 written\n")

# ============================================================
# RATE TESTS (12-17)
# ============================================================
cat("\n--- RATE Tests ---\n")

# Fit causal forest for RATE tests
cat("Fitting causal forest for RATE tests...\n")
set.seed(42)
cf <- causal_forest(X, Y_causal, W, num.trees = 500, seed = 42)
tau_hat <- predict(cf)$predictions
DR_scores <- get_scores(cf)

# Test 12: AUTOC with CATE predictions as priorities
cat("Test 12: RATE AUTOC with CATE priorities\n")
set.seed(42)
rate12 <- rank_average_treatment_effect.fit(DR_scores, tau_hat,
                                             target = "AUTOC",
                                             R = 200)
lines12 <- c(
  "TEST 12: RATE AUTOC CATE priorities",
  "STATUS: OK",
  sprintf("estimate: %.8f", rate12$estimate),
  sprintf("std_err: %.8f", rate12$std.err),
  sprintf("target: %s", rate12$target)
)
writeLines(lines12, "/tmp/grf_stata/tests/fidelity_reports/15_vi_rate_tune/r_test12.txt")
cat(sprintf("  AUTOC estimate=%.4f, SE=%.4f\n", rate12$estimate, rate12$std.err))
cat("  [DONE] Test 12 written\n")

# Test 13: QINI
cat("Test 13: RATE QINI\n")
set.seed(42)
rate13 <- rank_average_treatment_effect.fit(DR_scores, tau_hat,
                                             target = "QINI",
                                             R = 200)
lines13 <- c(
  "TEST 13: RATE QINI CATE priorities",
  "STATUS: OK",
  sprintf("estimate: %.8f", rate13$estimate),
  sprintf("std_err: %.8f", rate13$std.err),
  sprintf("target: %s", rate13$target)
)
writeLines(lines13, "/tmp/grf_stata/tests/fidelity_reports/15_vi_rate_tune/r_test13.txt")
cat(sprintf("  QINI estimate=%.4f, SE=%.4f\n", rate13$estimate, rate13$std.err))
cat("  [DONE] Test 13 written\n")

# Test 14: Custom priorities = X1 (drives tau)
cat("Test 14: RATE AUTOC custom priorities = X1\n")
set.seed(42)
rate14 <- rank_average_treatment_effect.fit(DR_scores, X[, 1],
                                             target = "AUTOC",
                                             R = 200)
lines14 <- c(
  "TEST 14: RATE AUTOC custom priorities X1",
  "STATUS: OK",
  sprintf("estimate: %.8f", rate14$estimate),
  sprintf("std_err: %.8f", rate14$std.err),
  sprintf("target: %s", rate14$target)
)
writeLines(lines14, "/tmp/grf_stata/tests/fidelity_reports/15_vi_rate_tune/r_test14.txt")
cat(sprintf("  AUTOC (X1 priority) estimate=%.4f, SE=%.4f\n", rate14$estimate, rate14$std.err))
cat("  [DONE] Test 14 written\n")

# Test 15: Random priorities (RATE near 0)
cat("Test 15: RATE AUTOC random priorities (should be near 0)\n")
set.seed(99)
rand_priority <- rnorm(n)
set.seed(42)
rate15 <- rank_average_treatment_effect.fit(DR_scores, rand_priority,
                                             target = "AUTOC",
                                             R = 200)
lines15 <- c(
  "TEST 15: RATE AUTOC random priorities near 0",
  "STATUS: OK",
  sprintf("estimate: %.8f", rate15$estimate),
  sprintf("std_err: %.8f", rate15$std.err),
  sprintf("z_stat: %.8f", rate15$estimate / rate15$std.err),
  sprintf("target: %s", rate15$target)
)
writeLines(lines15, "/tmp/grf_stata/tests/fidelity_reports/15_vi_rate_tune/r_test15.txt")
cat(sprintf("  AUTOC (random priority) estimate=%.4f, SE=%.4f\n", rate15$estimate, rate15$std.err))
cat("  [DONE] Test 15 written\n")

# Test 16: Bootstrap=500
cat("Test 16: RATE AUTOC bootstrap=500\n")
set.seed(42)
rate16 <- rank_average_treatment_effect.fit(DR_scores, tau_hat,
                                             target = "AUTOC",
                                             R = 500)
lines16 <- c(
  "TEST 16: RATE AUTOC bootstrap=500",
  "STATUS: OK",
  sprintf("estimate: %.8f", rate16$estimate),
  sprintf("std_err: %.8f", rate16$std.err),
  sprintf("target: %s", rate16$target)
)
writeLines(lines16, "/tmp/grf_stata/tests/fidelity_reports/15_vi_rate_tune/r_test16.txt")
cat(sprintf("  AUTOC (B=500) estimate=%.4f, SE=%.4f\n", rate16$estimate, rate16$std.err))
cat("  [DONE] Test 16 written\n")

# Test 17: RATE with sample.weights (debiasing.weights in Stata)
cat("Test 17: RATE with sample.weights\n")
set.seed(42)
deb_weights <- runif(n, 0.5, 1.5)
set.seed(42)
rate17 <- rank_average_treatment_effect.fit(DR_scores, tau_hat,
                                             target = "AUTOC",
                                             R = 200,
                                             sample.weights = deb_weights)
lines17 <- c(
  "TEST 17: RATE AUTOC with debiasing.weights",
  "STATUS: OK",
  sprintf("estimate: %.8f", rate17$estimate),
  sprintf("std_err: %.8f", rate17$std.err),
  sprintf("target: %s", rate17$target)
)
writeLines(lines17, "/tmp/grf_stata/tests/fidelity_reports/15_vi_rate_tune/r_test17.txt")
cat(sprintf("  AUTOC (debiasing weights) estimate=%.4f, SE=%.4f\n", rate17$estimate, rate17$std.err))
cat("  [DONE] Test 17 written\n")

# ============================================================
# TUNE PARAMETER TESTS (18-24)
# ============================================================
cat("\n--- Tune Parameter Tests ---\n")

# Test 18: Tune mtry + minnodesize (regression_forest)
cat("Test 18: Tune mtry + minnodesize (regression_forest)\n")
set.seed(42)
rf18 <- regression_forest(X, Y,
                           tune.parameters = c("mtry", "min.node.size"),
                           tune.num.trees = 200,
                           tune.num.reps = 50,
                           seed = 42)
tuned18 <- rf18$tuning.output
pred18 <- predict(rf18, X_test)$predictions
cat(sprintf("  Tuned mtry: %d, min.node.size: %d\n",
            rf18$`_tuning.output`$params["mtry"],
            rf18$`_tuning.output`$params["min.node.size"]))
cat(sprintf("  (access via tuning.output)\n"))
# Extract tuned params
tp18 <- tuned18$params
lines18 <- c(
  "TEST 18: Tune mtry + minnodesize regression_forest",
  "STATUS: OK",
  sprintf("tuned_mtry: %s", tp18["mtry"]),
  sprintf("tuned_minnodesize: %s", tp18["min.node.size"]),
  sprintf("pred_mean: %.8f", mean(pred18)),
  sprintf("pred_sd: %.8f", sd(pred18))
)
# Write predictions
pred18_str <- paste(sprintf("%.8f", pred18), collapse = "\n")
writeLines(pred18_str, "/tmp/grf_stata/tests/fidelity_reports/15_vi_rate_tune/r_test18_preds.txt")
writeLines(lines18, "/tmp/grf_stata/tests/fidelity_reports/15_vi_rate_tune/r_test18.txt")
cat("  [DONE] Test 18 written\n")

# Test 19: Tune samplefrac
cat("Test 19: Tune samplefrac\n")
set.seed(42)
rf19 <- regression_forest(X, Y,
                           tune.parameters = c("sample.fraction"),
                           tune.num.trees = 200,
                           tune.num.reps = 50,
                           seed = 42)
tp19 <- rf19$tuning.output$params
pred19 <- predict(rf19, X_test)$predictions
lines19 <- c(
  "TEST 19: Tune samplefrac regression_forest",
  "STATUS: OK",
  sprintf("tuned_samplefrac: %s", tp19["sample.fraction"]),
  sprintf("pred_mean: %.8f", mean(pred19)),
  sprintf("pred_sd: %.8f", sd(pred19))
)
pred19_str <- paste(sprintf("%.8f", pred19), collapse = "\n")
writeLines(pred19_str, "/tmp/grf_stata/tests/fidelity_reports/15_vi_rate_tune/r_test19_preds.txt")
writeLines(lines19, "/tmp/grf_stata/tests/fidelity_reports/15_vi_rate_tune/r_test19.txt")
cat("  [DONE] Test 19 written\n")

# Test 20: Tune all params
cat("Test 20: Tune all params\n")
set.seed(42)
rf20_tune <- regression_forest(X, Y,
                                tune.parameters = "all",
                                tune.num.trees = 200,
                                tune.num.reps = 50,
                                seed = 42)
tp20 <- rf20_tune$tuning.output$params
pred20 <- predict(rf20_tune, X_test)$predictions
lines20 <- c(
  "TEST 20: Tune all params regression_forest",
  "STATUS: OK"
)
for (pnm in names(tp20)) {
  lines20 <- c(lines20, sprintf("tuned_%s: %s", gsub("\\.", "_", pnm), tp20[pnm]))
}
lines20 <- c(lines20,
              sprintf("pred_mean: %.8f", mean(pred20)),
              sprintf("pred_sd: %.8f", sd(pred20)))
pred20_str <- paste(sprintf("%.8f", pred20), collapse = "\n")
writeLines(pred20_str, "/tmp/grf_stata/tests/fidelity_reports/15_vi_rate_tune/r_test20_preds.txt")
writeLines(lines20, "/tmp/grf_stata/tests/fidelity_reports/15_vi_rate_tune/r_test20.txt")
cat("  [DONE] Test 20 written\n")

# Test 21: Tune on causal_forest
cat("Test 21: Tune causal_forest\n")
set.seed(42)
cf21 <- causal_forest(X, Y_causal, W,
                       tune.parameters = c("mtry", "min.node.size"),
                       tune.num.trees = 200,
                       tune.num.reps = 50,
                       seed = 42)
tp21 <- cf21$tuning.output$params
tau21 <- predict(cf21)$predictions
lines21 <- c(
  "TEST 21: Tune mtry + minnodesize causal_forest",
  "STATUS: OK",
  sprintf("tuned_mtry: %s", tp21["mtry"]),
  sprintf("tuned_minnodesize: %s", tp21["min.node.size"]),
  sprintf("tau_mean: %.8f", mean(tau21)),
  sprintf("tau_sd: %.8f", sd(tau21))
)
writeLines(lines21, "/tmp/grf_stata/tests/fidelity_reports/15_vi_rate_tune/r_test21.txt")
cat("  [DONE] Test 21 written\n")

# Test 22: tunenumtrees=100
cat("Test 22: tunenumtrees=100\n")
set.seed(42)
rf22 <- regression_forest(X, Y,
                           tune.parameters = c("mtry", "min.node.size"),
                           tune.num.trees = 100,
                           tune.num.reps = 50,
                           seed = 42)
tp22 <- rf22$tuning.output$params
pred22 <- predict(rf22, X_test)$predictions
lines22 <- c(
  "TEST 22: Tune tunenumtrees=100",
  "STATUS: OK",
  sprintf("tuned_mtry: %s", tp22["mtry"]),
  sprintf("tuned_minnodesize: %s", tp22["min.node.size"]),
  sprintf("pred_mean: %.8f", mean(pred22)),
  sprintf("pred_sd: %.8f", sd(pred22))
)
pred22_str <- paste(sprintf("%.8f", pred22), collapse = "\n")
writeLines(pred22_str, "/tmp/grf_stata/tests/fidelity_reports/15_vi_rate_tune/r_test22_preds.txt")
writeLines(lines22, "/tmp/grf_stata/tests/fidelity_reports/15_vi_rate_tune/r_test22.txt")
cat("  [DONE] Test 22 written\n")

# Test 23: tunenumreps=10
cat("Test 23: tunenumreps=10\n")
set.seed(42)
rf23 <- regression_forest(X, Y,
                           tune.parameters = c("mtry", "min.node.size"),
                           tune.num.trees = 200,
                           tune.num.reps = 10,
                           seed = 42)
tp23 <- rf23$tuning.output$params
pred23 <- predict(rf23, X_test)$predictions
lines23 <- c(
  "TEST 23: Tune tunenumreps=10",
  "STATUS: OK",
  sprintf("tuned_mtry: %s", tp23["mtry"]),
  sprintf("tuned_minnodesize: %s", tp23["min.node.size"]),
  sprintf("pred_mean: %.8f", mean(pred23)),
  sprintf("pred_sd: %.8f", sd(pred23))
)
pred23_str <- paste(sprintf("%.8f", pred23), collapse = "\n")
writeLines(pred23_str, "/tmp/grf_stata/tests/fidelity_reports/15_vi_rate_tune/r_test23_preds.txt")
writeLines(lines23, "/tmp/grf_stata/tests/fidelity_reports/15_vi_rate_tune/r_test23.txt")
cat("  [DONE] Test 23 written\n")

# Test 24: Tune improves MSE
cat("Test 24: Does tuning improve MSE?\n")
set.seed(42)
# Default (no tune)
rf_default <- regression_forest(X, Y, num.trees = 2000, seed = 42)
pred_default <- predict(rf_default, X_test)$predictions
mse_default <- mean((Y_test - pred_default)^2)

# Tuned
rf_tuned <- regression_forest(X, Y,
                               tune.parameters = "all",
                               tune.num.trees = 200,
                               tune.num.reps = 50,
                               num.trees = 2000,
                               seed = 42)
pred_tuned <- predict(rf_tuned, X_test)$predictions
mse_tuned <- mean((Y_test - pred_tuned)^2)

tuned_better <- mse_tuned <= mse_default * 1.05  # tuned within 5% of default or better
cat(sprintf("  MSE default: %.6f, MSE tuned: %.6f\n", mse_default, mse_tuned))
lines24 <- c(
  "TEST 24: Tune improves MSE on held-out data",
  sprintf("STATUS: %s", ifelse(tuned_better, "PASS", "FAIL")),
  sprintf("mse_default: %.8f", mse_default),
  sprintf("mse_tuned: %.8f", mse_tuned),
  sprintf("tuned_better: %d", as.integer(tuned_better))
)
pred24_str <- paste(sprintf("%.8f", pred_tuned), collapse = "\n")
writeLines(pred24_str, "/tmp/grf_stata/tests/fidelity_reports/15_vi_rate_tune/r_test24_preds.txt")
writeLines(lines24, "/tmp/grf_stata/tests/fidelity_reports/15_vi_rate_tune/r_test24.txt")
cat("  [DONE] Test 24 written\n")

# ---- Summary ----
cat("\n=== R Tests Complete ===\n")
cat("All results written to /tmp/grf_stata/tests/fidelity_reports/15_vi_rate_tune/\n")

# Write DGP data for Stata
df <- data.frame(X)
df$y <- Y
df$y_causal <- Y_causal
df$W <- W
df$tau <- tau
df$cluster <- cluster
df$weights_vec <- weights_vec
df$rand_priority <- c(rnorm(n))  # Note: seed 99 in R
set.seed(99)
df$rand_priority <- rnorm(n)

# Test set
df_test <- data.frame(X_test)
colnames(df_test) <- paste0("x", 1:p)
df_test$y_test <- Y_test

# Write CSV
write.csv(df, "/tmp/grf_stata/tests/fidelity_reports/15_vi_rate_tune/data_dgp.csv",
          row.names = FALSE)
write.csv(df_test, "/tmp/grf_stata/tests/fidelity_reports/15_vi_rate_tune/data_test.csv",
          row.names = FALSE)

cat("DGP data written to data_dgp.csv and data_test.csv\n")
