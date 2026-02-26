# benchmark_r.R -- R grf speed benchmarks matching test_speed.do
# Run from test-c-plugin-skill-grf/ directory

library(grf)
cat("\n============================================================\n")
cat(" GRF R Package Speed Benchmarks (grf", as.character(packageVersion("grf")), ")\n")
cat("============================================================\n\n")

# ================================================================
# PART 1: Forest type comparison (n=2000, p=10, trees=500)
# ================================================================
cat("--- Part 1: Forest Type Comparison (n=2000, p=10, trees=500) ---\n\n")

set.seed(42)
n <- 2000
p <- 10
X <- matrix(rnorm(n * p), n, p)
Y <- 2*X[,1] + X[,2]^2 + 0.5*X[,3]*X[,4] + rnorm(n)
W <- as.numeric(X[,1] + rnorm(n) > 0)
Z <- as.numeric(X[,2] + rnorm(n) > 0)
Y <- Y + W*(1 + X[,5]) + Z*0.3
Y2 <- X[,1] + X[,3] + rnorm(n)
time_surv <- exp(0.5*X[,1] + rnorm(n))
status_surv <- as.numeric(runif(n) > 0.3)
W_multi <- floor(3*runif(n))

cat(sprintf("  %-40s %8s\n", "Forest Type", "Time (sec)"))
cat(sprintf("  %s\n", paste(rep("-", 55), collapse="")))

# 1. Regression forest
t <- system.time(regression_forest(X, Y, num.trees=500, seed=42))
cat(sprintf("  %-40s %8.2f sec\n", "Regression Forest", t["elapsed"]))

# 2. Causal forest
t <- system.time(causal_forest(X, Y, W, num.trees=500, seed=42))
cat(sprintf("  %-40s %8.2f sec\n", "Causal Forest", t["elapsed"]))

# 3. Quantile forest
t <- system.time(quantile_forest(X, Y, quantiles=c(0.25,0.5,0.75), num.trees=500, seed=42))
cat(sprintf("  %-40s %8.2f sec\n", "Quantile Forest (3 quantiles)", t["elapsed"]))

# 4. Instrumental forest
t <- system.time(instrumental_forest(X, Y, W, Z, num.trees=500, seed=42))
cat(sprintf("  %-40s %8.2f sec\n", "Instrumental Forest", t["elapsed"]))

# 5. Probability forest
t <- system.time(probability_forest(X, as.factor(as.numeric(W > 0)), num.trees=500, seed=42))
cat(sprintf("  %-40s %8.2f sec\n", "Probability Forest", t["elapsed"]))

# 6. Survival forest
t <- system.time(survival_forest(X, time_surv, status_surv, num.trees=500, seed=42))
cat(sprintf("  %-40s %8.2f sec\n", "Survival Forest", t["elapsed"]))

# 7. Causal survival forest
horizon <- median(time_surv[status_surv == 1])
t <- system.time(causal_survival_forest(X, time_surv, W, status_surv, horizon=horizon, num.trees=500, seed=42))
cat(sprintf("  %-40s %8.2f sec\n", "Causal Survival Forest", t["elapsed"]))

# 8. Multi-arm causal forest
t <- system.time(multi_arm_causal_forest(X, Y, as.factor(W_multi), num.trees=500, seed=42))
cat(sprintf("  %-40s %8.2f sec\n", "Multi-arm Causal Forest", t["elapsed"]))

# 9. Multi-regression forest
t <- system.time(multi_regression_forest(X, cbind(Y, Y2), num.trees=500, seed=42))
cat(sprintf("  %-40s %8.2f sec\n", "Multi-regression Forest", t["elapsed"]))

# 10. Local linear regression forest
t <- system.time({
  f <- ll_regression_forest(X, Y, num.trees=500, seed=42)
  predict(f)
})
cat(sprintf("  %-40s %8.2f sec\n", "LL Regression Forest", t["elapsed"]))

# 11. LL regression with LL splits
t <- system.time({
  f <- ll_regression_forest(X, Y, num.trees=500, seed=42, enable.ll.split=TRUE)
  predict(f)
})
cat(sprintf("  %-40s %8.2f sec\n", "LL Regression (LL splits)", t["elapsed"]))

# 12. Boosted regression forest
t <- system.time(boosted_regression_forest(X, Y, num.trees=500, seed=42, boost.steps=3))
cat(sprintf("  %-40s %8.2f sec\n", "Boosted Regression (3 steps)", t["elapsed"]))

# 13. LM forest
t <- system.time(lm_forest(X, Y, W, num.trees=500, seed=42))
cat(sprintf("  %-40s %8.2f sec\n", "LM Forest (1 regressor)", t["elapsed"]))


# ================================================================
# PART 2: Scaling with N
# ================================================================
cat("\n--- Part 2: Scaling with N (Regression Forest, p=10, trees=500) ---\n\n")
cat(sprintf("  %8s %12s %12s\n", "N", "Time (sec)", "Time/N (ms)"))
cat(sprintf("  %s\n", paste(rep("-", 40), collapse="")))

for (nn in c(500, 1000, 2000, 5000, 10000)) {
  set.seed(42)
  X <- matrix(rnorm(nn * 10), nn, 10)
  Y <- 2*X[,1] + X[,2]^2 + rnorm(nn)
  t <- system.time(regression_forest(X, Y, num.trees=500, seed=42))
  cat(sprintf("  %8d %12.3f %12.3f\n", nn, t["elapsed"], t["elapsed"]/nn*1000))
}

cat("\n============================================================\n")
cat(" R benchmarks completed\n")
cat("============================================================\n")
