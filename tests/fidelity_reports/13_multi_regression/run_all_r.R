## run_all_r.R
## Generate all R reference predictions for multi_regression_forest fidelity tests
## Outputs: CSV files with R predictions + covariates

library(grf)
OUTDIR <- "/tmp/grf_stata/tests/fidelity_reports/13_multi_regression"

## ---- Shared DGP ----
set.seed(42)
n <- 500; p <- 5
X <- matrix(rnorm(n * p), n, p)
colnames(X) <- paste0("x", 1:p)
Y1 <- X[,1] + X[,2] + rnorm(n)
Y2 <- -X[,1] + X[,3] + rnorm(n)
Y3 <- X[,2] * X[,3] + rnorm(n)
Y <- cbind(Y1, Y2, Y3)

## Cluster IDs (50 clusters)
cluster_ids <- rep(1:50, each = 10)
## Weights (uniform positive)
set.seed(99)
weights <- runif(n, 0.5, 2.0)

cat("=== GRF multi_regression_forest fidelity R scripts ===\n")

## ================================================================
## Test 01: 2 outcomes (Y1, Y2)
## ================================================================
cat("\n[01] 2 outcomes ...\n")
rf <- multi_regression_forest(X, Y[,1:2], num.trees = 500, seed = 42)
preds <- predict(rf)$predictions   # n x 2 matrix
df <- as.data.frame(X)
df$y1 <- Y[,1]; df$y2 <- Y[,2]
df$r_pred_y1 <- preds[,1]; df$r_pred_y2 <- preds[,2]
write.csv(df, file.path(OUTDIR, "test01_data.csv"), row.names = FALSE)
cat("  Y1 pred range:", range(preds[,1]), "\n")
cat("  Y2 pred range:", range(preds[,2]), "\n")

## ================================================================
## Test 02: 3 outcomes (Y1, Y2, Y3)
## ================================================================
cat("\n[02] 3 outcomes ...\n")
rf <- multi_regression_forest(X, Y, num.trees = 500, seed = 42)
preds <- predict(rf)$predictions
df <- as.data.frame(X)
df$y1 <- Y[,1]; df$y2 <- Y[,2]; df$y3 <- Y[,3]
df$r_pred_y1 <- preds[,1]; df$r_pred_y2 <- preds[,2]; df$r_pred_y3 <- preds[,3]
write.csv(df, file.path(OUTDIR, "test02_data.csv"), row.names = FALSE)
cat("  Y1/Y2/Y3 pred ranges:", range(preds[,1]), range(preds[,2]), range(preds[,3]), "\n")

## ================================================================
## Test 03: 5 outcomes
## ================================================================
cat("\n[03] 5 outcomes ...\n")
set.seed(42)
Y4 <- X[,4] - X[,5] + rnorm(n)
Y5 <- X[,1]^2 + rnorm(n)
Y5out <- cbind(Y1, Y2, Y3, Y4, Y5)
rf <- multi_regression_forest(X, Y5out, num.trees = 500, seed = 42)
preds <- predict(rf)$predictions
df <- as.data.frame(X)
for (k in 1:5) { df[[paste0("y",k)]] <- Y5out[,k]; df[[paste0("r_pred_y",k)]] <- preds[,k] }
write.csv(df, file.path(OUTDIR, "test03_data.csv"), row.names = FALSE)
cat("  Pred ranges for 5 outcomes done.\n")

## ================================================================
## Test 04: Correlated outcomes (Y2 = Y1 + noise)
## ================================================================
cat("\n[04] Correlated outcomes ...\n")
set.seed(42)
Y_base <- X[,1] + X[,2] + rnorm(n)
Y_corr2 <- Y_base + 0.5 * rnorm(n)   # highly correlated with Y_base
Ycorr <- cbind(Y_base, Y_corr2)
rf <- multi_regression_forest(X, Ycorr, num.trees = 500, seed = 42)
preds <- predict(rf)$predictions
df <- as.data.frame(X)
df$y1 <- Ycorr[,1]; df$y2 <- Ycorr[,2]
df$r_pred_y1 <- preds[,1]; df$r_pred_y2 <- preds[,2]
write.csv(df, file.path(OUTDIR, "test04_data.csv"), row.names = FALSE)
cat("  Correlation of Y1,Y2:", cor(Ycorr[,1], Ycorr[,2]), "\n")

## ================================================================
## Test 05: Independent outcomes
## ================================================================
cat("\n[05] Independent outcomes ...\n")
set.seed(42)
Yind1 <- X[,1] + rnorm(n)
Yind2 <- X[,5] + rnorm(n)
Yind <- cbind(Yind1, Yind2)
rf <- multi_regression_forest(X, Yind, num.trees = 500, seed = 42)
preds <- predict(rf)$predictions
df <- as.data.frame(X)
df$y1 <- Yind[,1]; df$y2 <- Yind[,2]
df$r_pred_y1 <- preds[,1]; df$r_pred_y2 <- preds[,2]
write.csv(df, file.path(OUTDIR, "test05_data.csv"), row.names = FALSE)
cat("  Correlation of Y1,Y2:", cor(Yind[,1], Yind[,2]), "\n")

## ================================================================
## Test 06: cluster()
## ================================================================
cat("\n[06] cluster() ...\n")
rf <- multi_regression_forest(X, Y[,1:2], num.trees = 500, seed = 42,
                               clusters = cluster_ids)
preds <- predict(rf)$predictions
df <- as.data.frame(X)
df$y1 <- Y[,1]; df$y2 <- Y[,2]
df$cluster_id <- cluster_ids
df$r_pred_y1 <- preds[,1]; df$r_pred_y2 <- preds[,2]
write.csv(df, file.path(OUTDIR, "test06_data.csv"), row.names = FALSE)
cat("  Cluster pred ranges done.\n")

## ================================================================
## Test 07: weights()
## ================================================================
cat("\n[07] weights() ...\n")
rf <- multi_regression_forest(X, Y[,1:2], num.trees = 500, seed = 42,
                               sample.weights = weights)
preds <- predict(rf)$predictions
df <- as.data.frame(X)
df$y1 <- Y[,1]; df$y2 <- Y[,2]
df$wt <- weights
df$r_pred_y1 <- preds[,1]; df$r_pred_y2 <- preds[,2]
write.csv(df, file.path(OUTDIR, "test07_data.csv"), row.names = FALSE)
cat("  Weighted pred ranges done.\n")

## ================================================================
## Test 08: nohonesty
## ================================================================
cat("\n[08] nohonesty ...\n")
rf <- multi_regression_forest(X, Y[,1:2], num.trees = 500, seed = 42,
                               honesty = FALSE)
preds <- predict(rf)$predictions
df <- as.data.frame(X)
df$y1 <- Y[,1]; df$y2 <- Y[,2]
df$r_pred_y1 <- preds[,1]; df$r_pred_y2 <- preds[,2]
write.csv(df, file.path(OUTDIR, "test08_data.csv"), row.names = FALSE)
cat("  No-honesty pred ranges done.\n")

## ================================================================
## Test 09: mtry=2
## ================================================================
cat("\n[09] mtry=2 ...\n")
rf <- multi_regression_forest(X, Y[,1:2], num.trees = 500, seed = 42,
                               mtry = 2)
preds <- predict(rf)$predictions
df <- as.data.frame(X)
df$y1 <- Y[,1]; df$y2 <- Y[,2]
df$r_pred_y1 <- preds[,1]; df$r_pred_y2 <- preds[,2]
write.csv(df, file.path(OUTDIR, "test09_data.csv"), row.names = FALSE)
cat("  mtry=2 pred ranges done.\n")

## ================================================================
## Test 10: minnodesize=20
## ================================================================
cat("\n[10] minnodesize=20 ...\n")
rf <- multi_regression_forest(X, Y[,1:2], num.trees = 500, seed = 42,
                               min.node.size = 20)
preds <- predict(rf)$predictions
df <- as.data.frame(X)
df$y1 <- Y[,1]; df$y2 <- Y[,2]
df$r_pred_y1 <- preds[,1]; df$r_pred_y2 <- preds[,2]
write.csv(df, file.path(OUTDIR, "test10_data.csv"), row.names = FALSE)
cat("  minnodesize=20 pred ranges done.\n")

## ================================================================
## Test 11: samplefrac=0.3
## ================================================================
cat("\n[11] samplefrac=0.3 ...\n")
rf <- multi_regression_forest(X, Y[,1:2], num.trees = 500, seed = 42,
                               sample.fraction = 0.3)
preds <- predict(rf)$predictions
df <- as.data.frame(X)
df$y1 <- Y[,1]; df$y2 <- Y[,2]
df$r_pred_y1 <- preds[,1]; df$r_pred_y2 <- preds[,2]
write.csv(df, file.path(OUTDIR, "test11_data.csv"), row.names = FALSE)
cat("  samplefrac=0.3 pred ranges done.\n")

## ================================================================
## Test 12: Combined cluster + weights + nohonesty
## ================================================================
cat("\n[12] Combined: cluster + weights + nohonesty ...\n")
rf <- multi_regression_forest(X, Y[,1:2], num.trees = 500, seed = 42,
                               clusters = cluster_ids,
                               sample.weights = weights,
                               honesty = FALSE)
preds <- predict(rf)$predictions
df <- as.data.frame(X)
df$y1 <- Y[,1]; df$y2 <- Y[,2]
df$cluster_id <- cluster_ids
df$wt <- weights
df$r_pred_y1 <- preds[,1]; df$r_pred_y2 <- preds[,2]
write.csv(df, file.path(OUTDIR, "test12_data.csv"), row.names = FALSE)
cat("  Combined pred ranges done.\n")

## ================================================================
## Test 13: Linear outcomes (Y = X*beta, should predict well)
## ================================================================
cat("\n[13] Linear outcomes ...\n")
set.seed(42)
beta1 <- c(2, -1, 0.5, 0, 0)
beta2 <- c(-1, 0.5, 2, -0.5, 0)
Ylin1 <- X %*% beta1 + 0.1 * rnorm(n)   # very low noise
Ylin2 <- X %*% beta2 + 0.1 * rnorm(n)
Ylin <- cbind(Ylin1, Ylin2)
rf <- multi_regression_forest(X, Ylin, num.trees = 500, seed = 42)
preds <- predict(rf)$predictions
df <- as.data.frame(X)
df$y1 <- Ylin[,1]; df$y2 <- Ylin[,2]
df$r_pred_y1 <- preds[,1]; df$r_pred_y2 <- preds[,2]
write.csv(df, file.path(OUTDIR, "test13_data.csv"), row.names = FALSE)
r_cor_y1 <- cor(preds[,1], Ylin[,1])
r_cor_y2 <- cor(preds[,2], Ylin[,2])
cat("  Linear Y1 R R-vs-truth cor:", r_cor_y1, "\n")
cat("  Linear Y2 R R-vs-truth cor:", r_cor_y2, "\n")

## ================================================================
## Test 14: Nonlinear outcomes (Y = sin(X))
## ================================================================
cat("\n[14] Nonlinear outcomes ...\n")
set.seed(42)
Ynl1 <- sin(X[,1]) + cos(X[,2]) + rnorm(n, sd = 0.5)
Ynl2 <- X[,1] * X[,2] + sin(X[,3]) + rnorm(n, sd = 0.5)
Ynl <- cbind(Ynl1, Ynl2)
rf <- multi_regression_forest(X, Ynl, num.trees = 500, seed = 42)
preds <- predict(rf)$predictions
df <- as.data.frame(X)
df$y1 <- Ynl[,1]; df$y2 <- Ynl[,2]
df$r_pred_y1 <- preds[,1]; df$r_pred_y2 <- preds[,2]
write.csv(df, file.path(OUTDIR, "test14_data.csv"), row.names = FALSE)
r_cor_y1 <- cor(preds[,1], Ynl[,1])
r_cor_y2 <- cor(preds[,2], Ynl[,2])
cat("  Nonlinear Y1 R R-vs-truth cor:", r_cor_y1, "\n")
cat("  Nonlinear Y2 R R-vs-truth cor:", r_cor_y2, "\n")

cat("\n=== All R reference predictions generated ===\n")
