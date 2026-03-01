#!/usr/bin/env Rscript
# R fidelity tests for causal_forest
# Generates data, runs all 19 tests, saves CSVs for Stata comparison

library(grf)

WORKDIR <- "/tmp/grf_stata/tests/fidelity_reports/02_causal_fit"
dir.create(WORKDIR, recursive = TRUE, showWarnings = FALSE)

results <- list()

# ─────────────────────────────────────────────
# Helper: run one causal forest test, return list of predictions + metadata
# ─────────────────────────────────────────────
run_cf <- function(tag, X, Y, W, ...) {
  cat(sprintf("\n=== TEST %s ===\n", tag))
  cf <- causal_forest(X, Y, W, ...)
  preds <- predict(cf)$predictions
  list(tag = tag, preds = preds, cf = cf, X = X, Y = Y, W = W)
}

# ─────────────────────────────────────────────
# Test 1 – Default (n=500, p=5)
# ─────────────────────────────────────────────
set.seed(42)
n <- 500; p <- 5
X <- matrix(rnorm(n * p), n, p)
colnames(X) <- paste0("x", 1:p)
W <- rbinom(n, 1, 0.5)
tau <- X[,1] + X[,2]
Y <- X[,1] + tau * W + rnorm(n)

res1 <- run_cf("01_default", X, Y, W, num.trees = 500, seed = 42)
df1 <- as.data.frame(X)
df1$y <- Y; df1$w <- W; df1$tau_r <- res1$preds; df1$tau_true <- tau
write.csv(df1, file.path(WORKDIR, "test01_default.csv"), row.names = FALSE)
results[["01_default"]] <- list(tag="01_default", cor=1.0, n=n)

# ─────────────────────────────────────────────
# Test 2 – nostabilizesplits (stabilize.splits=FALSE)
# ─────────────────────────────────────────────
set.seed(42)
res2 <- run_cf("02_nostabilize", X, Y, W,
               num.trees = 500, seed = 42, stabilize.splits = FALSE)
df2 <- as.data.frame(X)
df2$y <- Y; df2$w <- W; df2$tau_r <- res2$preds; df2$tau_true <- tau
write.csv(df2, file.path(WORKDIR, "test02_nostabilize.csv"), row.names = FALSE)

# ─────────────────────────────────────────────
# Test 3 – nuisancetrees=100 (pre-compute nuisance with 100 trees, then supply)
# In grf R, there is no num.trees.for.nuisance; we simulate by fitting
# the nuisance forests with 100 trees externally, then passing Y.hat/W.hat
# ─────────────────────────────────────────────
set.seed(42)
rf_y3 <- regression_forest(X, Y, num.trees = 100, seed = 42)
rf_w3 <- regression_forest(X, W, num.trees = 100, seed = 42)
yhat3 <- predict(rf_y3)$predictions
what3 <- predict(rf_w3)$predictions
res3 <- causal_forest(X, Y, W,
                      num.trees = 500, seed = 42,
                      Y.hat = yhat3, W.hat = what3)
preds3 <- predict(res3)$predictions
df3 <- as.data.frame(X)
df3$y <- Y; df3$w <- W; df3$tau_r <- preds3; df3$tau_true <- tau
df3$yhat <- yhat3; df3$what <- what3
write.csv(df3, file.path(WORKDIR, "test03_nuisancetrees100.csv"), row.names = FALSE)

# ─────────────────────────────────────────────
# Test 4 – User-supplied Y.hat
# ─────────────────────────────────────────────
set.seed(42)
rf_y <- regression_forest(X, Y, num.trees = 500, seed = 42)
yhat_supplied <- predict(rf_y)$predictions

res4 <- causal_forest(X, Y, W,
                      num.trees = 500, seed = 42,
                      Y.hat = yhat_supplied)
preds4 <- predict(res4)$predictions

df4 <- as.data.frame(X)
df4$y <- Y; df4$w <- W; df4$tau_r <- preds4; df4$tau_true <- tau
df4$yhat <- yhat_supplied
write.csv(df4, file.path(WORKDIR, "test04_yhat_supplied.csv"), row.names = FALSE)

# ─────────────────────────────────────────────
# Test 5 – User-supplied W.hat (propensity)
# ─────────────────────────────────────────────
set.seed(42)
rf_w <- regression_forest(X, W, num.trees = 500, seed = 42)
what_supplied <- predict(rf_w)$predictions

res5 <- causal_forest(X, Y, W,
                      num.trees = 500, seed = 42,
                      W.hat = what_supplied)
preds5 <- predict(res5)$predictions

df5 <- as.data.frame(X)
df5$y <- Y; df5$w <- W; df5$tau_r <- preds5; df5$tau_true <- tau
df5$what <- what_supplied
write.csv(df5, file.path(WORKDIR, "test05_what_supplied.csv"), row.names = FALSE)

# ─────────────────────────────────────────────
# Test 6 – Both Y.hat and W.hat supplied
# ─────────────────────────────────────────────
set.seed(42)
res6 <- causal_forest(X, Y, W,
                      num.trees = 500, seed = 42,
                      Y.hat = yhat_supplied,
                      W.hat = what_supplied)
preds6 <- predict(res6)$predictions

df6 <- as.data.frame(X)
df6$y <- Y; df6$w <- W; df6$tau_r <- preds6; df6$tau_true <- tau
df6$yhat <- yhat_supplied; df6$what <- what_supplied
write.csv(df6, file.path(WORKDIR, "test06_both_nuisance.csv"), row.names = FALSE)

# ─────────────────────────────────────────────
# Test 7 – estimate.variance
# ─────────────────────────────────────────────
set.seed(42)
res7 <- causal_forest(X, Y, W,
                      num.trees = 500, seed = 42,
                      compute.oob.predictions = TRUE)
preds7 <- predict(res7, estimate.variance = TRUE)
tau7 <- preds7$predictions
var7 <- preds7$variance.estimates

df7 <- as.data.frame(X)
df7$y <- Y; df7$w <- W; df7$tau_r <- tau7; df7$var_r <- var7; df7$tau_true <- tau
write.csv(df7, file.path(WORKDIR, "test07_variance.csv"), row.names = FALSE)

# ─────────────────────────────────────────────
# Test 8 – cluster()
# ─────────────────────────────────────────────
set.seed(42)
clusters <- rep(1:10, each = 50)
res8 <- causal_forest(X, Y, W,
                      num.trees = 500, seed = 42,
                      clusters = clusters)
preds8 <- predict(res8)$predictions

df8 <- as.data.frame(X)
df8$y <- Y; df8$w <- W; df8$tau_r <- preds8; df8$tau_true <- tau; df8$cluster <- clusters
write.csv(df8, file.path(WORKDIR, "test08_clusters.csv"), row.names = FALSE)

# ─────────────────────────────────────────────
# Test 9 – weights()
# ─────────────────────────────────────────────
set.seed(42)
wts <- runif(n, 0.5, 2.0)
res9 <- causal_forest(X, Y, W,
                      num.trees = 500, seed = 42,
                      sample.weights = wts)
preds9 <- predict(res9)$predictions

df9 <- as.data.frame(X)
df9$y <- Y; df9$w <- W; df9$tau_r <- preds9; df9$tau_true <- tau; df9$wts <- wts
write.csv(df9, file.path(WORKDIR, "test09_weights.csv"), row.names = FALSE)

# ─────────────────────────────────────────────
# Test 10 – equalize.cluster.weights
# ─────────────────────────────────────────────
set.seed(42)
res10 <- causal_forest(X, Y, W,
                       num.trees = 500, seed = 42,
                       clusters = clusters,
                       equalize.cluster.weights = TRUE)
preds10 <- predict(res10)$predictions

df10 <- as.data.frame(X)
df10$y <- Y; df10$w <- W; df10$tau_r <- preds10; df10$tau_true <- tau; df10$cluster <- clusters
write.csv(df10, file.path(WORKDIR, "test10_eq_cluster.csv"), row.names = FALSE)

# ─────────────────────────────────────────────
# Test 11 – nohonesty
# ─────────────────────────────────────────────
set.seed(42)
res11 <- causal_forest(X, Y, W,
                       num.trees = 500, seed = 42,
                       honesty = FALSE)
preds11 <- predict(res11)$predictions

df11 <- as.data.frame(X)
df11$y <- Y; df11$w <- W; df11$tau_r <- preds11; df11$tau_true <- tau
write.csv(df11, file.path(WORKDIR, "test11_nohonesty.csv"), row.names = FALSE)

# ─────────────────────────────────────────────
# Test 12 – mtry=2 + min.node.size=20
# ─────────────────────────────────────────────
set.seed(42)
res12 <- causal_forest(X, Y, W,
                       num.trees = 500, seed = 42,
                       mtry = 2, min.node.size = 20)
preds12 <- predict(res12)$predictions

df12 <- as.data.frame(X)
df12$y <- Y; df12$w <- W; df12$tau_r <- preds12; df12$tau_true <- tau
write.csv(df12, file.path(WORKDIR, "test12_mtry_minnodesize.csv"), row.names = FALSE)

# ─────────────────────────────────────────────
# Test 13 – sample.fraction=0.3
# ─────────────────────────────────────────────
set.seed(42)
res13 <- causal_forest(X, Y, W,
                       num.trees = 500, seed = 42,
                       sample.fraction = 0.3)
preds13 <- predict(res13)$predictions

df13 <- as.data.frame(X)
df13$y <- Y; df13$w <- W; df13$tau_r <- preds13; df13$tau_true <- tau
write.csv(df13, file.path(WORKDIR, "test13_samplefrac.csv"), row.names = FALSE)

# ─────────────────────────────────────────────
# Test 14 – Unbalanced treatment (80/20 split)
# ─────────────────────────────────────────────
set.seed(42)
W14 <- rbinom(n, 1, 0.2)
Y14 <- X[,1] + tau * W14 + rnorm(n)

res14 <- causal_forest(X, Y14, W14,
                       num.trees = 500, seed = 42)
preds14 <- predict(res14)$predictions

df14 <- as.data.frame(X)
df14$y <- Y14; df14$w <- W14; df14$tau_r <- preds14; df14$tau_true <- tau
write.csv(df14, file.path(WORKDIR, "test14_unbalanced.csv"), row.names = FALSE)

# ─────────────────────────────────────────────
# Test 15 – Strong heterogeneity: tau(X) = 5*X1*(X2>0)
# ─────────────────────────────────────────────
set.seed(42)
tau15 <- 5 * X[,1] * (X[,2] > 0)
Y15 <- X[,1] + tau15 * W + rnorm(n)

res15 <- causal_forest(X, Y15, W,
                       num.trees = 500, seed = 42)
preds15 <- predict(res15)$predictions

df15 <- as.data.frame(X)
df15$y <- Y15; df15$w <- W; df15$tau_r <- preds15; df15$tau_true <- tau15
write.csv(df15, file.path(WORKDIR, "test15_strong_heterogeneity.csv"), row.names = FALSE)

# ─────────────────────────────────────────────
# Test 16 – No treatment effect (tau=0)
# ─────────────────────────────────────────────
set.seed(42)
Y16 <- X[,1] + rnorm(n)   # tau = 0

res16 <- causal_forest(X, Y16, W,
                       num.trees = 500, seed = 42)
preds16 <- predict(res16)$predictions

df16 <- as.data.frame(X)
df16$y <- Y16; df16$w <- W; df16$tau_r <- preds16; df16$tau_true <- rep(0, n)
write.csv(df16, file.path(WORKDIR, "test16_no_effect.csv"), row.names = FALSE)

# ─────────────────────────────────────────────
# Test 17 – yhatgen + whatgen (save nuisance estimates)
#   We use the default run from test 1 and check the nuisance saves
# ─────────────────────────────────────────────
# In R: get nuisance from existing forest
yhat17 <- res1$cf$Y.hat
what17 <- res1$cf$W.hat

df17 <- as.data.frame(X)
df17$y <- Y; df17$w <- W; df17$tau_r <- res1$preds
df17$yhat_r <- yhat17; df17$what_r <- what17; df17$tau_true <- tau
write.csv(df17, file.path(WORKDIR, "test17_nuisancegen.csv"), row.names = FALSE)

# ─────────────────────────────────────────────
# Test 18 – Large n=2000
# ─────────────────────────────────────────────
set.seed(42)
n18 <- 2000
X18 <- matrix(rnorm(n18 * p), n18, p)
colnames(X18) <- paste0("x", 1:p)
W18 <- rbinom(n18, 1, 0.5)
tau18 <- X18[,1] + X18[,2]
Y18 <- X18[,1] + tau18 * W18 + rnorm(n18)

res18 <- causal_forest(X18, Y18, W18,
                       num.trees = 500, seed = 42)
preds18 <- predict(res18)$predictions

df18 <- as.data.frame(X18)
df18$y <- Y18; df18$w <- W18; df18$tau_r <- preds18; df18$tau_true <- tau18
write.csv(df18, file.path(WORKDIR, "test18_largen.csv"), row.names = FALSE)

# ─────────────────────────────────────────────
# Test 19 – Combined: cluster + weights + estimate.variance + nostabilizesplits
# ─────────────────────────────────────────────
set.seed(42)
res19 <- causal_forest(X, Y, W,
                       num.trees = 500, seed = 42,
                       clusters = clusters,
                       sample.weights = wts,
                       stabilize.splits = FALSE)
preds19_tau <- predict(res19, estimate.variance = TRUE)
preds19 <- preds19_tau$predictions
var19 <- preds19_tau$variance.estimates

df19 <- as.data.frame(X)
df19$y <- Y; df19$w <- W; df19$tau_r <- preds19; df19$var_r <- var19
df19$tau_true <- tau; df19$cluster <- clusters; df19$wts <- wts
write.csv(df19, file.path(WORKDIR, "test19_combined.csv"), row.names = FALSE)

cat("\n\nAll R tests completed. CSVs written to:\n", WORKDIR, "\n")
cat("\nR session info:\n")
cat(sprintf("  R version: %s\n", R.version.string))
cat(sprintf("  grf version: %s\n", packageVersion("grf")))
