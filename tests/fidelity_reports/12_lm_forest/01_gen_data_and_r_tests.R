## Generate shared DGP and run all R lm_forest tests
library(grf)

set.seed(42); n <- 500; p <- 5
X <- matrix(rnorm(n * p), n, p)
W1 <- rnorm(n)
W2 <- rnorm(n)
W3 <- rnorm(n)  # extra treatment for K=3
beta1 <- X[, 1] + X[, 2]
beta2 <- -X[, 1] + 0.5 * X[, 3]
beta3 <- 0.3 * X[, 2] - 0.4 * X[, 4]  # for K=3
Y <- X[, 1] + beta1 * W1 + beta2 * W2 + rnorm(n)
W <- cbind(W1, W2)

## Save base dataset as CSV for Stata
df <- data.frame(X, W1, W2, W3, Y,
                 cluster_id = rep(1:50, each = 10),
                 obs_weight = runif(n, 0.5, 1.5))
write.csv(df, "/tmp/grf_stata/tests/fidelity_reports/12_lm_forest/data.csv",
          row.names = FALSE)
cat("Data saved\n")

## Helper to save predictions
save_preds <- function(forest, path, prefix = "beta") {
  preds <- predict(forest)$predictions
  # predictions may be 3D array (n, K, 1) or 2D (n, K) or vector (n,)
  if (length(dim(preds)) == 3) preds <- preds[, , 1]
  if (is.null(dim(preds))) preds <- matrix(preds, ncol = 1)
  df_out <- as.data.frame(preds)
  # Ensure columns named beta_1, beta_2, ...
  names(df_out) <- paste0(prefix, "_", seq_len(ncol(df_out)))
  write.csv(df_out, path, row.names = FALSE)
}

## ---- Test 1: Single W (K=1) ----
cat("Test 1: single W\n")
f1 <- lm_forest(X, Y, W1, num.trees = 500, seed = 42)
save_preds(f1, "/tmp/grf_stata/tests/fidelity_reports/12_lm_forest/r_t01.csv")

## ---- Test 2: Two W (K=2) ----
cat("Test 2: two W\n")
f2 <- lm_forest(X, Y, W, num.trees = 500, seed = 42)
save_preds(f2, "/tmp/grf_stata/tests/fidelity_reports/12_lm_forest/r_t02.csv")

## ---- Test 3: Three W (K=3) ----
cat("Test 3: three W\n")
W3mat <- cbind(W1, W2, W3)
Y3 <- X[, 1] + beta1 * W1 + beta2 * W2 + beta3 * W3 + rnorm(n)
f3 <- lm_forest(X, Y3, W3mat, num.trees = 500, seed = 42)
save_preds(f3, "/tmp/grf_stata/tests/fidelity_reports/12_lm_forest/r_t03.csv")

## ---- Test 4: gradient.weights ----
cat("Test 4: gradient.weights\n")
f4 <- lm_forest(X, Y, W, num.trees = 500, seed = 42,
                gradient.weights = c(0.7, 0.3))
save_preds(f4, "/tmp/grf_stata/tests/fidelity_reports/12_lm_forest/r_t04.csv")

## ---- Test 5: stabilize.splits ON ----
cat("Test 5: stabilize.splits\n")
f5 <- lm_forest(X, Y, W, num.trees = 500, seed = 42,
                stabilize.splits = TRUE)
save_preds(f5, "/tmp/grf_stata/tests/fidelity_reports/12_lm_forest/r_t05.csv")

## ---- Test 6: User-supplied Y.hat ----
cat("Test 6: user Y.hat\n")
Yhat <- predict(regression_forest(X, Y, num.trees = 500, seed = 42))$predictions
f6 <- lm_forest(X, Y, W, num.trees = 500, seed = 42, Y.hat = Yhat)
save_preds(f6, "/tmp/grf_stata/tests/fidelity_reports/12_lm_forest/r_t06.csv")
write.csv(data.frame(Yhat = Yhat),
          "/tmp/grf_stata/tests/fidelity_reports/12_lm_forest/r_yhat.csv",
          row.names = FALSE)

## ---- Test 7: User-supplied W.hat ----
cat("Test 7: user W.hat\n")
What1 <- predict(regression_forest(X, W1, num.trees = 500, seed = 42))$predictions
What2 <- predict(regression_forest(X, W2, num.trees = 500, seed = 42))$predictions
f7 <- lm_forest(X, Y, W, num.trees = 500, seed = 42,
                W.hat = cbind(What1, What2))
save_preds(f7, "/tmp/grf_stata/tests/fidelity_reports/12_lm_forest/r_t07.csv")
write.csv(data.frame(What1 = What1, What2 = What2),
          "/tmp/grf_stata/tests/fidelity_reports/12_lm_forest/r_what.csv",
          row.names = FALSE)

## ---- Test 8: Both Y.hat and W.hat ----
cat("Test 8: both Y.hat and W.hat\n")
f8 <- lm_forest(X, Y, W, num.trees = 500, seed = 42,
                Y.hat = Yhat, W.hat = cbind(What1, What2))
save_preds(f8, "/tmp/grf_stata/tests/fidelity_reports/12_lm_forest/r_t08.csv")

## ---- Test 9: cluster ----
cat("Test 9: cluster\n")
cluster_ids <- rep(1:50, each = 10)
f9 <- lm_forest(X, Y, W, num.trees = 500, seed = 42,
                clusters = cluster_ids)
save_preds(f9, "/tmp/grf_stata/tests/fidelity_reports/12_lm_forest/r_t09.csv")

## ---- Test 10: weights ----
cat("Test 10: weights\n")
set.seed(123)
obs_weights <- runif(n, 0.5, 1.5)
f10 <- lm_forest(X, Y, W, num.trees = 500, seed = 42,
                 sample.weights = obs_weights)
save_preds(f10, "/tmp/grf_stata/tests/fidelity_reports/12_lm_forest/r_t10.csv")

## ---- Test 11: nohonesty ----
cat("Test 11: nohonesty\n")
f11 <- lm_forest(X, Y, W, num.trees = 500, seed = 42,
                 honesty = FALSE)
save_preds(f11, "/tmp/grf_stata/tests/fidelity_reports/12_lm_forest/r_t11.csv")

## ---- Test 12: mtry=2 ----
cat("Test 12: mtry=2\n")
f12 <- lm_forest(X, Y, W, num.trees = 500, seed = 42, mtry = 2)
save_preds(f12, "/tmp/grf_stata/tests/fidelity_reports/12_lm_forest/r_t12.csv")

## ---- Test 13: minnodesize=20 ----
cat("Test 13: minnodesize=20\n")
f13 <- lm_forest(X, Y, W, num.trees = 500, seed = 42, min.node.size = 20)
save_preds(f13, "/tmp/grf_stata/tests/fidelity_reports/12_lm_forest/r_t13.csv")

## ---- Test 14: nuisancetrees=100 ----
cat("Test 14: nuisancetrees=100\n")
f14 <- lm_forest(X, Y, W, num.trees = 500, seed = 42,
                 num.threads = 1)
# Note: nuisance.trees argument may not exist in all grf versions
# Try with it, catch if not supported
tryCatch({
  f14 <- lm_forest(X, Y, W, num.trees = 500, seed = 42,
                   nuisance.trees = 100)
  save_preds(f14, "/tmp/grf_stata/tests/fidelity_reports/12_lm_forest/r_t14.csv")
  cat("nuisance.trees=100 succeeded\n")
}, error = function(e) {
  cat("nuisance.trees not supported, using default\n")
  save_preds(f14, "/tmp/grf_stata/tests/fidelity_reports/12_lm_forest/r_t14.csv")
})

## ---- Test 15: estimate.variance ----
cat("Test 15: estimate.variance\n")
f15 <- lm_forest(X, Y, W, num.trees = 500, seed = 42,
                 ci.group.size = 2)
preds15 <- predict(f15, estimate.variance = TRUE)
pmat <- preds15$predictions
# predictions may be 3D array (n, K, 1) or 2D (n, K)
if (length(dim(pmat)) == 3) pmat <- pmat[, , 1]
df15 <- data.frame(
  beta_1 = pmat[, 1],
  beta_2 = pmat[, 2]
)
if (!is.null(preds15$variance.estimates)) {
  vars <- preds15$variance.estimates
  if (length(dim(vars)) == 3) vars <- vars[, , 1]
  df15$var_1 <- vars[, 1]
  df15$var_2 <- vars[, 2]
}
write.csv(df15, "/tmp/grf_stata/tests/fidelity_reports/12_lm_forest/r_t15.csv",
          row.names = FALSE)

## ---- Test 16: Homogeneous beta ----
cat("Test 16: homogeneous beta\n")
set.seed(42)
beta_const <- 2
Y16 <- X[, 1] + beta_const * W1 + rnorm(n)
f16 <- lm_forest(X, Y16, W1, num.trees = 500, seed = 42)
p16 <- predict(f16)$predictions
if (length(dim(p16)) == 3) p16 <- p16[, , 1]
if (is.null(dim(p16))) p16 <- matrix(p16, ncol = 1)
write.csv(data.frame(Y16 = Y16, true_beta = beta_const,
                     pred_beta1 = p16[, 1]),
          "/tmp/grf_stata/tests/fidelity_reports/12_lm_forest/r_t16.csv",
          row.names = FALSE)
cat("Homogeneous beta: mean =", mean(p16[, 1]),
    "sd =", sd(p16[, 1]), "\n")

## ---- Test 17: Strong heterogeneity ----
cat("Test 17: strong heterogeneity\n")
set.seed(42)
beta_strong <- 5 * (X[, 1] > 0)
Y17 <- X[, 1] + beta_strong * W1 + rnorm(n)
f17 <- lm_forest(X, Y17, W1, num.trees = 500, seed = 42)
p17 <- predict(f17)$predictions
if (length(dim(p17)) == 3) p17 <- p17[, , 1]
if (is.null(dim(p17))) p17 <- matrix(p17, ncol = 1)
write.csv(data.frame(Y17 = Y17, true_beta = beta_strong,
                     pred_beta1 = p17[, 1]),
          "/tmp/grf_stata/tests/fidelity_reports/12_lm_forest/r_t17.csv",
          row.names = FALSE)
cat("Strong heterogeneity: cor(true, pred) =",
    cor(beta_strong, p17[, 1]), "\n")

cat("All R tests complete\n")
