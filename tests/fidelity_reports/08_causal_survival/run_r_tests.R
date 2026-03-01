#!/usr/bin/env Rscript
# R fidelity test for causal_survival_forest (grf 2.5.0)
# Tests 1-16 as specified

suppressPackageStartupMessages(library(grf))

# ---- DGP ----
set.seed(42); n <- 500; p <- 5
X <- matrix(rnorm(n * p), n, p)
W <- rbinom(n, 1, 0.5)
T_control <- rexp(n, rate = exp(0.3 * X[, 1]))
tau_survival <- X[, 1] * 0.5
T_treated <- rexp(n, rate = exp(0.3 * X[, 1] - tau_survival * W))
T_true <- ifelse(W == 1, T_treated, T_control)
C <- rexp(n, rate = 0.2)
Y <- pmin(T_true, C)
D <- as.integer(T_true <= C)
horizon <- median(Y)

cat("DGP summary:\n")
cat(sprintf("  n=%d, p=%d, events=%d (%.1f%%), horizon=%.4f\n",
            n, p, sum(D), mean(D) * 100, horizon))
cat(sprintf("  mean(W)=%.3f\n", mean(W)))

# ---- Helper to save predictions ----
save_preds <- function(tau_hat, filename) {
  df <- data.frame(tau_hat = as.numeric(tau_hat))
  write.csv(df, filename, row.names = FALSE)
  cat(sprintf("  Saved %d predictions -> %s\n", nrow(df), filename))
}

outdir <- "/tmp/grf_stata/tests/fidelity_reports/08_causal_survival"

results <- list()

run_test <- function(test_id, desc, expr, outfile) {
  cat(sprintf("\n=== Test %d: %s ===\n", test_id, desc))
  t0 <- proc.time()
  tryCatch({
    csf <- eval(expr)
    elapsed <- (proc.time() - t0)[["elapsed"]]
    tau_hat <- predict(csf)$predictions
    save_preds(tau_hat, file.path(outdir, outfile))
    cat(sprintf("  OK  (%.1fs)  mean=%.4f sd=%.4f\n",
                elapsed, mean(tau_hat), sd(tau_hat)))
    results[[test_id]] <<- list(status = "OK", desc = desc,
                                 mean = mean(tau_hat), sd = sd(tau_hat),
                                 n = length(tau_hat))
  }, error = function(e) {
    cat(sprintf("  ERROR: %s\n", conditionMessage(e)))
    results[[test_id]] <<- list(status = "ERROR", desc = desc,
                                 error = conditionMessage(e))
  })
}

# Test 1: Default (RMST target, default horizon = median failure time)
run_test(1, "Default RMST, horizon=median(Y)",
  quote(causal_survival_forest(X, Y, W, D,
          num.trees = 500, seed = 42,
          horizon = horizon, target = "RMST")),
  "r_pred_01.csv")

# Test 2: Explicit horizon = Q75 of Y
horizon_q75 <- quantile(Y, 0.75)
run_test(2, sprintf("Explicit horizon=Q75(Y)=%.4f", horizon_q75),
  bquote(causal_survival_forest(X, Y, W, D,
           num.trees = 500, seed = 42,
           horizon = .(horizon_q75), target = "RMST")),
  "r_pred_02.csv")

# Test 3: Survival probability target
run_test(3, "Survival probability target",
  quote(causal_survival_forest(X, Y, W, D,
          num.trees = 500, seed = 42,
          horizon = horizon, target = "survival.probability")),
  "r_pred_03.csv")

# Test 4: nostabilizesplits
run_test(4, "No stabilize splits",
  quote(causal_survival_forest(X, Y, W, D,
          num.trees = 500, seed = 42,
          horizon = horizon, target = "RMST",
          stabilize.splits = FALSE)),
  "r_pred_04.csv")

# Test 5: User-supplied W.hat
W.hat <- rep(0.5, n)  # known propensity
run_test(5, "User-supplied W.hat=0.5",
  quote(causal_survival_forest(X, Y, W, D,
          num.trees = 500, seed = 42,
          horizon = horizon, target = "RMST",
          W.hat = W.hat)),
  "r_pred_05.csv")

# Test 6: With cluster()
set.seed(99)
clusters <- sample(1:50, n, replace = TRUE)
run_test(6, "With cluster()",
  quote(causal_survival_forest(X, Y, W, D,
          num.trees = 500, seed = 42,
          horizon = horizon, target = "RMST",
          clusters = clusters)),
  "r_pred_06.csv")

# Test 7: With sample.weights
wts <- runif(n, 0.5, 1.5)
run_test(7, "With observation weights",
  quote(causal_survival_forest(X, Y, W, D,
          num.trees = 500, seed = 42,
          horizon = horizon, target = "RMST",
          sample.weights = wts)),
  "r_pred_07.csv")

# Test 8: nohonesty
run_test(8, "No honesty",
  quote(causal_survival_forest(X, Y, W, D,
          num.trees = 500, seed = 42,
          horizon = horizon, target = "RMST",
          honesty = FALSE)),
  "r_pred_08.csv")

# Test 9: mtry=2
run_test(9, "mtry=2",
  quote(causal_survival_forest(X, Y, W, D,
          num.trees = 500, seed = 42,
          horizon = horizon, target = "RMST",
          mtry = 2)),
  "r_pred_09.csv")

# Test 10: min.node.size=20
run_test(10, "min.node.size=20",
  quote(causal_survival_forest(X, Y, W, D,
          num.trees = 500, seed = 42,
          horizon = horizon, target = "RMST",
          min.node.size = 20)),
  "r_pred_10.csv")

# Test 11: Heavy censoring (80%)
set.seed(42)
C_heavy <- rexp(n, rate = 0.8)
Y_hc <- pmin(T_true, C_heavy)
D_hc <- as.integer(T_true <= C_heavy)
horizon_hc <- median(Y_hc[D_hc == 1])
cat(sprintf("\nHeavy censoring: events=%d (%.1f%%), horizon=%.4f\n",
            sum(D_hc), mean(D_hc) * 100, horizon_hc))
run_test(11, sprintf("Heavy censoring (~%.0f%% censored)", (1 - mean(D_hc)) * 100),
  bquote(causal_survival_forest(X, Y_hc, W, D_hc,
           num.trees = 500, seed = 42,
           horizon = .(horizon_hc), target = "RMST")),
  "r_pred_11.csv")

# Test 12: Balanced treatment 50/50 (same as default DGP, but explicit)
run_test(12, "Balanced treatment 50/50",
  quote(causal_survival_forest(X, Y, W, D,
          num.trees = 500, seed = 42,
          horizon = horizon, target = "RMST")),
  "r_pred_12.csv")

# Test 13: Unbalanced treatment 70/30
set.seed(123)
W_unbal <- rbinom(n, 1, 0.7)
T_treated_unbal <- rexp(n, rate = exp(0.3 * X[, 1] - tau_survival * W_unbal))
T_true_unbal <- ifelse(W_unbal == 1, T_treated_unbal, T_control)
Y_unbal <- pmin(T_true_unbal, C)
D_unbal <- as.integer(T_true_unbal <= C)
horizon_unbal <- median(Y_unbal)
cat(sprintf("\nUnbalanced: mean(W)=%.3f, events=%d (%.1f%%), horizon=%.4f\n",
            mean(W_unbal), sum(D_unbal), mean(D_unbal) * 100, horizon_unbal))
run_test(13, "Unbalanced treatment 70/30",
  bquote(causal_survival_forest(X, Y_unbal, W_unbal, D_unbal,
           num.trees = 500, seed = 42,
           horizon = .(horizon_unbal), target = "RMST")),
  "r_pred_13.csv")

# Test 14: nuisance forests with fewer trees (num.trees for all steps)
# In grf, there's no direct nuisancetrees param but we use smaller num.trees
run_test(14, "Fewer nuisance trees (num.trees=100)",
  quote(causal_survival_forest(X, Y, W, D,
          num.trees = 100, seed = 42,
          horizon = horizon, target = "RMST")),
  "r_pred_14.csv")

# Test 15: Small horizon (Q25 of failure times)
horizon_q25 <- quantile(Y[D == 1], 0.25)
cat(sprintf("\nSmall horizon Q25 of failure times: %.4f\n", horizon_q25))
run_test(15, sprintf("Small horizon=Q25(failure)=%.4f", horizon_q25),
  bquote(causal_survival_forest(X, Y, W, D,
           num.trees = 500, seed = 42,
           horizon = .(horizon_q25), target = "RMST")),
  "r_pred_15.csv")

# Test 16: Large horizon (Q90 of failure times)
horizon_q90 <- quantile(Y[D == 1], 0.90)
cat(sprintf("\nLarge horizon Q90 of failure times: %.4f\n", horizon_q90))
run_test(16, sprintf("Large horizon=Q90(failure)=%.4f", horizon_q90),
  bquote(causal_survival_forest(X, Y, W, D,
           num.trees = 500, seed = 42,
           horizon = .(horizon_q90), target = "RMST")),
  "r_pred_16.csv")

# ---- Save DGP constants for Stata ----
dgp_constants <- data.frame(
  n = n,
  p = p,
  horizon = horizon,
  horizon_q75 = horizon_q75,
  horizon_q25 = horizon_q25,
  horizon_q90 = horizon_q90,
  horizon_hc = horizon_hc,
  horizon_unbal = horizon_unbal
)
write.csv(dgp_constants, file.path(outdir, "dgp_constants.csv"), row.names = FALSE)

# ---- Save full dataset for Stata ----
df_main <- as.data.frame(X)
colnames(df_main) <- paste0("x", 1:p)
df_main$time  <- Y
df_main$status <- D
df_main$w      <- W
df_main$cluster_var <- clusters
set.seed(99)
df_main$weights_var <- runif(n, 0.5, 1.5)
# unbalanced
df_main$w_unbal  <- W_unbal
df_main$time_unbal <- Y_unbal
df_main$status_unbal <- D_unbal
# heavy censoring
df_main$time_hc <- Y_hc
df_main$status_hc <- D_hc

write.csv(df_main, file.path(outdir, "test_data.csv"), row.names = FALSE)
cat("\nSaved test_data.csv with", nrow(df_main), "rows\n")

# ---- Summary ----
cat("\n\n===== R TEST SUMMARY =====\n")
for (tid in names(results)) {
  r <- results[[tid]]
  if (r$status == "OK") {
    cat(sprintf("  Test %2s: OK  | %s\n", tid, r$desc))
  } else {
    cat(sprintf("  Test %2s: FAIL| %s | %s\n", tid, r$desc, r$error))
  }
}
