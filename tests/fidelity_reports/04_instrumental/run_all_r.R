## Instrumental Forest Fidelity Tests -- R side
## Generates data CSVs with R predictions for each test scenario
library(grf)
outdir <- "/tmp/grf_stata/tests/fidelity_reports/04_instrumental"

## ============================================================
## DGP helper
## ============================================================
make_dgp <- function(seed = 42, n = 500, p = 5,
                     compliance_scale = 1.0,  # 1.0 = normal, 0.05 = weak, 2.0 = strong
                     tau_const = FALSE) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("x", 1:p)
  Z <- rbinom(n, 1, 0.5)
  compliance <- pnorm(X[, 1]) * compliance_scale
  compliance <- pmin(compliance, 0.99)
  W <- rbinom(n, 1, compliance * Z + (1 - compliance) * 0.1)
  if (tau_const) {
    tau <- rep(1.5, n)
  } else {
    tau <- X[, 1] + X[, 2]
  }
  Y <- X[, 1] + tau * W + X[, 3] * 0.5 + rnorm(n)
  list(X = X, Y = Y, W = W, Z = Z, tau = tau)
}

## ============================================================
## Test 01: Default options
## ============================================================
cat("=== Test 01: Default options ===\n")
d <- make_dgp()
fit <- instrumental_forest(d$X, d$Y, d$W, d$Z,
                           num.trees = 500, seed = 42)
preds <- predict(fit)$predictions

df <- as.data.frame(d$X)
df$y <- d$Y; df$w <- d$W; df$z <- d$Z
df$r_pred <- preds
write.csv(df, file.path(outdir, "test01_data.csv"), row.names = FALSE)
cat("Test 01 done. pred range:", range(preds), "\n\n")

## ============================================================
## Test 02: stabilize.splits=FALSE
## ============================================================
cat("=== Test 02: stabilize.splits=FALSE ===\n")
d <- make_dgp()
fit <- instrumental_forest(d$X, d$Y, d$W, d$Z,
                           num.trees = 500, seed = 42,
                           stabilize.splits = FALSE)
preds <- predict(fit)$predictions

df <- as.data.frame(d$X)
df$y <- d$Y; df$w <- d$W; df$z <- d$Z
df$r_pred <- preds
write.csv(df, file.path(outdir, "test02_data.csv"), row.names = FALSE)
cat("Test 02 done. pred range:", range(preds), "\n\n")

## ============================================================
## Test 03: reduced.form.weight=0.5
## ============================================================
cat("=== Test 03: reduced.form.weight=0.5 ===\n")
d <- make_dgp()
fit <- instrumental_forest(d$X, d$Y, d$W, d$Z,
                           num.trees = 500, seed = 42,
                           reduced.form.weight = 0.5)
preds <- predict(fit)$predictions

df <- as.data.frame(d$X)
df$y <- d$Y; df$w <- d$W; df$z <- d$Z
df$r_pred <- preds
write.csv(df, file.path(outdir, "test03_data.csv"), row.names = FALSE)
cat("Test 03 done. pred range:", range(preds), "\n\n")

## ============================================================
## Test 04: reduced.form.weight=1.0 (pure reduced form)
## ============================================================
cat("=== Test 04: reduced.form.weight=1.0 ===\n")
d <- make_dgp()
fit <- instrumental_forest(d$X, d$Y, d$W, d$Z,
                           num.trees = 500, seed = 42,
                           reduced.form.weight = 1.0)
preds <- predict(fit)$predictions

df <- as.data.frame(d$X)
df$y <- d$Y; df$w <- d$W; df$z <- d$Z
df$r_pred <- preds
write.csv(df, file.path(outdir, "test04_data.csv"), row.names = FALSE)
cat("Test 04 done. pred range:", range(preds), "\n\n")

## ============================================================
## Test 05: User-supplied Y.hat
## ============================================================
cat("=== Test 05: User-supplied Y.hat ===\n")
d <- make_dgp()
# Pre-compute nuisance via separate regression forests
rf_y <- regression_forest(d$X, d$Y, num.trees = 500, seed = 1)
Yhat <- predict(rf_y)$predictions
rf_w <- regression_forest(d$X, d$W, num.trees = 500, seed = 2)
What <- predict(rf_w)$predictions
rf_z <- regression_forest(d$X, d$Z, num.trees = 500, seed = 3)
Zhat <- predict(rf_z)$predictions

# Pass only Y.hat, let W.hat and Z.hat be estimated
# NOTE: R requires all three or none -- so supply all three here
fit <- instrumental_forest(d$X, d$Y, d$W, d$Z,
                           num.trees = 500, seed = 42,
                           Y.hat = Yhat,
                           W.hat = What,
                           Z.hat = Zhat)
preds <- predict(fit)$predictions

df <- as.data.frame(d$X)
df$y <- d$Y; df$w <- d$W; df$z <- d$Z
df$yhat_input <- Yhat
df$what_input <- What
df$zhat_input <- Zhat
df$r_pred <- preds
write.csv(df, file.path(outdir, "test05_data.csv"), row.names = FALSE)
cat("Test 05 done. pred range:", range(preds), "\n\n")

## ============================================================
## Test 06: User-supplied W.hat (same as test05 - all three supplied; here for
##           semantic labeling - W.hat focus)
## ============================================================
cat("=== Test 06: User-supplied W.hat ===\n")
d <- make_dgp()
rf_y <- regression_forest(d$X, d$Y, num.trees = 500, seed = 10)
Yhat <- predict(rf_y)$predictions
rf_w <- regression_forest(d$X, d$W, num.trees = 500, seed = 20)
What <- predict(rf_w)$predictions
rf_z <- regression_forest(d$X, d$Z, num.trees = 500, seed = 30)
Zhat <- predict(rf_z)$predictions

fit <- instrumental_forest(d$X, d$Y, d$W, d$Z,
                           num.trees = 500, seed = 42,
                           Y.hat = Yhat, W.hat = What, Z.hat = Zhat)
preds <- predict(fit)$predictions

df <- as.data.frame(d$X)
df$y <- d$Y; df$w <- d$W; df$z <- d$Z
df$yhat_input <- Yhat
df$what_input <- What
df$zhat_input <- Zhat
df$r_pred <- preds
write.csv(df, file.path(outdir, "test06_data.csv"), row.names = FALSE)
cat("Test 06 done. pred range:", range(preds), "\n\n")

## ============================================================
## Test 07: User-supplied Z.hat
## ============================================================
cat("=== Test 07: User-supplied Z.hat ===\n")
d <- make_dgp()
rf_y <- regression_forest(d$X, d$Y, num.trees = 500, seed = 100)
Yhat <- predict(rf_y)$predictions
rf_w <- regression_forest(d$X, d$W, num.trees = 500, seed = 200)
What <- predict(rf_w)$predictions
rf_z <- regression_forest(d$X, d$Z, num.trees = 500, seed = 300)
Zhat <- predict(rf_z)$predictions

fit <- instrumental_forest(d$X, d$Y, d$W, d$Z,
                           num.trees = 500, seed = 42,
                           Y.hat = Yhat, W.hat = What, Z.hat = Zhat)
preds <- predict(fit)$predictions

df <- as.data.frame(d$X)
df$y <- d$Y; df$w <- d$W; df$z <- d$Z
df$yhat_input <- Yhat
df$what_input <- What
df$zhat_input <- Zhat
df$r_pred <- preds
write.csv(df, file.path(outdir, "test07_data.csv"), row.names = FALSE)
cat("Test 07 done. pred range:", range(preds), "\n\n")

## ============================================================
## Test 08: All three nuisance supplied
## ============================================================
cat("=== Test 08: All three nuisance supplied ===\n")
d <- make_dgp()
# Use lm for simple nuisance
Yhat <- lm(d$Y ~ d$X)$fitted
What <- lm(d$W ~ d$X)$fitted
Zhat <- lm(d$Z ~ d$X)$fitted

fit <- instrumental_forest(d$X, d$Y, d$W, d$Z,
                           num.trees = 500, seed = 42,
                           Y.hat = Yhat, W.hat = What, Z.hat = Zhat)
preds <- predict(fit)$predictions

df <- as.data.frame(d$X)
df$y <- d$Y; df$w <- d$W; df$z <- d$Z
df$yhat_input <- Yhat
df$what_input <- What
df$zhat_input <- Zhat
df$r_pred <- preds
write.csv(df, file.path(outdir, "test08_data.csv"), row.names = FALSE)
cat("Test 08 done. pred range:", range(preds), "\n\n")

## ============================================================
## Test 09: estimate.variance
## ============================================================
cat("=== Test 09: estimate.variance ===\n")
d <- make_dgp()
fit <- instrumental_forest(d$X, d$Y, d$W, d$Z,
                           num.trees = 500, seed = 42)
preds_out <- predict(fit, estimate.variance = TRUE)
preds  <- preds_out$predictions
pvars  <- preds_out$variance.estimates

df <- as.data.frame(d$X)
df$y <- d$Y; df$w <- d$W; df$z <- d$Z
df$r_pred <- preds
df$r_var  <- pvars
write.csv(df, file.path(outdir, "test09_data.csv"), row.names = FALSE)
cat("Test 09 done. var range:", range(pvars), "\n\n")

## ============================================================
## Test 10: cluster()
## ============================================================
cat("=== Test 10: cluster() ===\n")
set.seed(42)
n <- 500; p <- 5
X <- matrix(rnorm(n * p), n, p)
colnames(X) <- paste0("x", 1:p)
clusters <- rep(1:50, each = 10)
Z <- rbinom(n, 1, 0.5)
compliance <- pnorm(X[, 1])
W <- rbinom(n, 1, compliance * Z + (1 - compliance) * 0.1)
tau <- X[, 1] + X[, 2]
Y <- X[, 1] + tau * W + X[, 3] * 0.5 + rnorm(n)

fit <- instrumental_forest(X, Y, W, Z,
                           num.trees = 500, seed = 42,
                           clusters = clusters)
preds <- predict(fit)$predictions

df <- as.data.frame(X)
df$y <- Y; df$w <- W; df$z <- Z; df$cluster_id <- clusters
df$r_pred <- preds
write.csv(df, file.path(outdir, "test10_data.csv"), row.names = FALSE)
cat("Test 10 done. pred range:", range(preds), "\n\n")

## ============================================================
## Test 11: weights()
## ============================================================
cat("=== Test 11: weights() ===\n")
d <- make_dgp()
set.seed(42)
wts <- runif(500, 0.5, 2.0)

fit <- instrumental_forest(d$X, d$Y, d$W, d$Z,
                           num.trees = 500, seed = 42,
                           sample.weights = wts)
preds <- predict(fit)$predictions

df <- as.data.frame(d$X)
df$y <- d$Y; df$w <- d$W; df$z <- d$Z; df$obs_wt <- wts
df$r_pred <- preds
write.csv(df, file.path(outdir, "test11_data.csv"), row.names = FALSE)
cat("Test 11 done. pred range:", range(preds), "\n\n")

## ============================================================
## Test 12: nohonesty
## ============================================================
cat("=== Test 12: nohonesty ===\n")
d <- make_dgp()
fit <- instrumental_forest(d$X, d$Y, d$W, d$Z,
                           num.trees = 500, seed = 42,
                           honesty = FALSE)
preds <- predict(fit)$predictions

df <- as.data.frame(d$X)
df$y <- d$Y; df$w <- d$W; df$z <- d$Z
df$r_pred <- preds
write.csv(df, file.path(outdir, "test12_data.csv"), row.names = FALSE)
cat("Test 12 done. pred range:", range(preds), "\n\n")

## ============================================================
## Test 13: mtry=2 + min.node.size=20
## ============================================================
cat("=== Test 13: mtry=2 + min.node.size=20 ===\n")
d <- make_dgp()
fit <- instrumental_forest(d$X, d$Y, d$W, d$Z,
                           num.trees = 500, seed = 42,
                           mtry = 2, min.node.size = 20)
preds <- predict(fit)$predictions

df <- as.data.frame(d$X)
df$y <- d$Y; df$w <- d$W; df$z <- d$Z
df$r_pred <- preds
write.csv(df, file.path(outdir, "test13_data.csv"), row.names = FALSE)
cat("Test 13 done. pred range:", range(preds), "\n\n")

## ============================================================
## Test 14: Strong instrument (high compliance)
## ============================================================
cat("=== Test 14: Strong instrument ===\n")
set.seed(42)
n <- 500; p <- 5
X <- matrix(rnorm(n * p), n, p)
colnames(X) <- paste0("x", 1:p)
Z <- rbinom(n, 1, 0.5)
# Strong: Z almost perfectly determines W
W <- rbinom(n, 1, ifelse(Z == 1, 0.95, 0.05))
tau <- X[, 1] + X[, 2]
Y <- X[, 1] + tau * W + X[, 3] * 0.5 + rnorm(n)

fit <- instrumental_forest(X, Y, W, Z,
                           num.trees = 500, seed = 42)
preds <- predict(fit)$predictions

df <- as.data.frame(X)
df$y <- Y; df$w <- W; df$z <- Z
df$r_pred <- preds
write.csv(df, file.path(outdir, "test14_data.csv"), row.names = FALSE)
cat("Test 14 done. pred range:", range(preds), "\n\n")

## ============================================================
## Test 15: Weak instrument
## ============================================================
cat("=== Test 15: Weak instrument ===\n")
set.seed(42)
n <- 500; p <- 5
X <- matrix(rnorm(n * p), n, p)
colnames(X) <- paste0("x", 1:p)
Z <- rbinom(n, 1, 0.5)
# Weak: Z barely affects W (low compliance ~5% uplift)
compliance <- 0.05
W <- rbinom(n, 1, compliance * Z + 0.1)
tau <- X[, 1] + X[, 2]
Y <- X[, 1] + tau * W + X[, 3] * 0.5 + rnorm(n)

fit <- instrumental_forest(X, Y, W, Z,
                           num.trees = 500, seed = 42)
preds <- predict(fit)$predictions

df <- as.data.frame(X)
df$y <- Y; df$w <- W; df$z <- Z
df$r_pred <- preds
write.csv(df, file.path(outdir, "test15_data.csv"), row.names = FALSE)
cat("Test 15 done. pred range:", range(preds), "\n\n")

## ============================================================
## Test 16: No heterogeneity (constant tau)
## ============================================================
cat("=== Test 16: No heterogeneity (constant tau) ===\n")
d <- make_dgp(tau_const = TRUE)
fit <- instrumental_forest(d$X, d$Y, d$W, d$Z,
                           num.trees = 500, seed = 42)
preds <- predict(fit)$predictions

df <- as.data.frame(d$X)
df$y <- d$Y; df$w <- d$W; df$z <- d$Z
df$r_pred <- preds
df$true_tau <- d$tau
write.csv(df, file.path(outdir, "test16_data.csv"), row.names = FALSE)
cat("Test 16 done. pred range:", range(preds), "\n\n")

## ============================================================
## Test 17: nuisancetrees=100 (Stata-only option)
## In R: pre-compute nuisance with 100-tree forests, then pass as Y.hat/W.hat/Z.hat
## ============================================================
cat("=== Test 17: nuisancetrees=100 ===\n")
d <- make_dgp()
# Approximate nuisancetrees(100) by pre-computing with 100-tree regression forests
rf_y <- regression_forest(d$X, d$Y, num.trees = 100, seed = 42)
Yhat_100 <- predict(rf_y)$predictions
rf_w <- regression_forest(d$X, d$W, num.trees = 100, seed = 42)
What_100 <- predict(rf_w)$predictions
rf_z <- regression_forest(d$X, d$Z, num.trees = 100, seed = 42)
Zhat_100 <- predict(rf_z)$predictions

fit <- instrumental_forest(d$X, d$Y, d$W, d$Z,
                           num.trees = 500, seed = 42,
                           Y.hat = Yhat_100, W.hat = What_100, Z.hat = Zhat_100)
preds <- predict(fit)$predictions

df <- as.data.frame(d$X)
df$y <- d$Y; df$w <- d$W; df$z <- d$Z
df$yhat_input <- Yhat_100
df$what_input <- What_100
df$zhat_input <- Zhat_100
df$r_pred <- preds
write.csv(df, file.path(outdir, "test17_data.csv"), row.names = FALSE)
cat("Test 17 done. pred range:", range(preds), "\n\n")

## ============================================================
## Test 18: yhatgen + whatgen + zhatgen (save nuisance estimates)
## ============================================================
cat("=== Test 18: Save nuisance estimates ===\n")
d <- make_dgp()
# Run with separate nuisance forests and capture their predictions
rf_y <- regression_forest(d$X, d$Y, num.trees = 500, seed = 42)
Yhat_saved <- predict(rf_y)$predictions
rf_w <- regression_forest(d$X, d$W, num.trees = 500, seed = 42)
What_saved <- predict(rf_w)$predictions
rf_z <- regression_forest(d$X, d$Z, num.trees = 500, seed = 42)
Zhat_saved <- predict(rf_z)$predictions

fit <- instrumental_forest(d$X, d$Y, d$W, d$Z,
                           num.trees = 500, seed = 42)
preds <- predict(fit)$predictions
# Extract nuisance from fit
Yhat_fit <- fit$Y.hat
What_fit <- fit$W.hat
Zhat_fit <- fit$Z.hat

df <- as.data.frame(d$X)
df$y <- d$Y; df$w <- d$W; df$z <- d$Z
df$r_pred <- preds
df$r_yhat <- Yhat_fit
df$r_what <- What_fit
df$r_zhat <- Zhat_fit
write.csv(df, file.path(outdir, "test18_data.csv"), row.names = FALSE)
cat("Test 18 done. Yhat range:", range(Yhat_fit), "\n")
cat("             What range:", range(What_fit), "\n")
cat("             Zhat range:", range(Zhat_fit), "\n\n")

cat("=== All R tests complete ===\n")
