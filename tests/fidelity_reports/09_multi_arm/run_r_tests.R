#!/usr/bin/env Rscript
# ============================================================
# R fidelity tests for multi_arm_causal_forest
# Generates data, runs all 16 tests, saves CSVs for Stata comparison
#
# NOTE: predict(multi_arm_causal_forest(...))$predictions is an array
#       of dimension n x K x 1, so we extract arm k via preds[, k, 1]
# ============================================================

library(grf)

WORKDIR <- "/tmp/grf_stata/tests/fidelity_reports/09_multi_arm"
dir.create(WORKDIR, recursive = TRUE, showWarnings = FALSE)

# Helper: extract arm k from n x K x 1 prediction array
get_arm <- function(preds_arr, k) preds_arr[, k, 1]

# ─────────────────────────────────────────────────────────────────────────────
# Shared DGP function (from spec)
# n=600, p=5, 3 arms: control(0), t1(1), t2(2)
# tau1 = X[,1] + X[,2], tau2 = -X[,1] + 2*X[,3]
# ─────────────────────────────────────────────────────────────────────────────
make_data <- function(seed = 42, n = 600, p = 5,
                      probs = c(1/3, 1/3, 1/3)) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("x", 1:p)
  W <- sample(0:2, n, replace = TRUE, prob = probs)
  W_factor <- factor(W)
  tau1 <- X[, 1] + X[, 2]
  tau2 <- -X[, 1] + 2 * X[, 3]
  Y <- X[, 1] + tau1 * (W == 1) + tau2 * (W == 2) + rnorm(n)
  list(X = X, Y = Y, W = W, W_factor = W_factor, tau1 = tau1, tau2 = tau2)
}

# ─────────────────────────────────────────────────────────────────────────────
# Helper: save a two-arm (K=2) result CSV
# ─────────────────────────────────────────────────────────────────────────────
save_macf <- function(tag, d, cf, extra_cols = list()) {
  preds_arr <- predict(cf)$predictions   # n x K x 1 array
  t1 <- get_arm(preds_arr, 1)
  t2 <- get_arm(preds_arr, 2)
  df <- as.data.frame(d$X)
  df$y  <- d$Y
  df$w  <- d$W
  df$w1 <- as.integer(d$W == 1)
  df$w2 <- as.integer(d$W == 2)
  df$tau_r_t1   <- t1
  df$tau_r_t2   <- t2
  df$tau_true_t1 <- d$tau1
  df$tau_true_t2 <- d$tau2
  for (nm in names(extra_cols)) df[[nm]] <- extra_cols[[nm]]
  outfile <- file.path(WORKDIR, sprintf("test%s.csv", tag))
  write.csv(df, outfile, row.names = FALSE)
  cat(sprintf("  Saved %s  (cor_t1=%.4f  cor_t2=%.4f)\n",
              basename(outfile),
              cor(t1, d$tau1), cor(t2, d$tau2)))
}

# ─────────────────────────────────────────────────────────────────────────────
# TEST 01 — Default 3-arm (2 treatment arms)
# ─────────────────────────────────────────────────────────────────────────────
cat("\n=== TEST 01: 2 treatment arms (default 3-arm) ===\n")
d1 <- make_data(seed = 42)
cf1 <- multi_arm_causal_forest(d1$X, d1$Y, d1$W_factor,
                                num.trees = 500, seed = 42)
save_macf("01_default_2arm", d1, cf1)

# ─────────────────────────────────────────────────────────────────────────────
# TEST 02 — Single treatment arm (binary, 2-arm study)
# W factor: 0=control, 1=treatment1 only
# ─────────────────────────────────────────────────────────────────────────────
cat("\n=== TEST 02: Single treatment arm (binary, 2-arm) ===\n")
set.seed(42)
n2 <- 600; p2 <- 5
X2 <- matrix(rnorm(n2 * p2), n2, p2)
colnames(X2) <- paste0("x", 1:p2)
W2_int <- rbinom(n2, 1, 0.5)
W2_factor <- factor(W2_int)
tau2_true <- X2[, 1] + X2[, 2]
Y2 <- X2[, 1] + tau2_true * W2_int + rnorm(n2)

cf2 <- multi_arm_causal_forest(X2, Y2, W2_factor,
                                num.trees = 500, seed = 42)
preds2_arr <- predict(cf2)$predictions  # n x 1 x 1
t2_1 <- get_arm(preds2_arr, 1)

df2 <- as.data.frame(X2)
df2$y <- Y2; df2$w <- W2_int; df2$w1 <- W2_int
df2$tau_r_t1    <- t2_1
df2$tau_true_t1 <- tau2_true
outfile2 <- file.path(WORKDIR, "test02_single_arm.csv")
write.csv(df2, outfile2, row.names = FALSE)
cat(sprintf("  Saved %s  (cor_t1=%.4f)\n", basename(outfile2),
            cor(t2_1, tau2_true)))

# ─────────────────────────────────────────────────────────────────────────────
# TEST 03 — 4 treatment arms (5-arm total: 0,1,2,3,4)
# ─────────────────────────────────────────────────────────────────────────────
cat("\n=== TEST 03: 4 treatment arms (5-arm total) ===\n")
set.seed(42)
n3 <- 600; p3 <- 5
X3 <- matrix(rnorm(n3 * p3), n3, p3)
colnames(X3) <- paste0("x", 1:p3)
W3_int <- sample(0:4, n3, replace = TRUE)
W3_factor <- factor(W3_int)
tau3_1 <-  X3[, 1] + X3[, 2]
tau3_2 <- -X3[, 1] + 2 * X3[, 3]
tau3_3 <-  X3[, 2] - X3[, 3]
tau3_4 <-  X3[, 1] - X3[, 4]
Y3 <- X3[, 1] +
  tau3_1 * (W3_int == 1) +
  tau3_2 * (W3_int == 2) +
  tau3_3 * (W3_int == 3) +
  tau3_4 * (W3_int == 4) + rnorm(n3)

cf3 <- multi_arm_causal_forest(X3, Y3, W3_factor,
                                num.trees = 500, seed = 42)
preds3_arr <- predict(cf3)$predictions  # n x 4 x 1
t3_1 <- get_arm(preds3_arr, 1); t3_2 <- get_arm(preds3_arr, 2)
t3_3 <- get_arm(preds3_arr, 3); t3_4 <- get_arm(preds3_arr, 4)

df3 <- as.data.frame(X3)
df3$y <- Y3; df3$w <- W3_int
df3$w1 <- as.integer(W3_int == 1); df3$w2 <- as.integer(W3_int == 2)
df3$w3 <- as.integer(W3_int == 3); df3$w4 <- as.integer(W3_int == 4)
df3$tau_r_t1 <- t3_1; df3$tau_r_t2 <- t3_2
df3$tau_r_t3 <- t3_3; df3$tau_r_t4 <- t3_4
df3$tau_true_t1 <- tau3_1; df3$tau_true_t2 <- tau3_2
df3$tau_true_t3 <- tau3_3; df3$tau_true_t4 <- tau3_4
outfile3 <- file.path(WORKDIR, "test03_4arms.csv")
write.csv(df3, outfile3, row.names = FALSE)
cat(sprintf("  Saved %s\n  (cor_t1=%.4f  cor_t2=%.4f  cor_t3=%.4f  cor_t4=%.4f)\n",
            basename(outfile3),
            cor(t3_1, tau3_1), cor(t3_2, tau3_2),
            cor(t3_3, tau3_3), cor(t3_4, tau3_4)))

# ─────────────────────────────────────────────────────────────────────────────
# TEST 04 — Unbalanced arms: 60/20/20 split
# ─────────────────────────────────────────────────────────────────────────────
cat("\n=== TEST 04: Unbalanced arms (60/20/20) ===\n")
d4 <- make_data(seed = 42, probs = c(0.6, 0.2, 0.2))
cf4 <- multi_arm_causal_forest(d4$X, d4$Y, d4$W_factor,
                                num.trees = 500, seed = 42)
save_macf("04_unbalanced", d4, cf4)

# ─────────────────────────────────────────────────────────────────────────────
# TEST 05 — nostabilizesplits
# ─────────────────────────────────────────────────────────────────────────────
cat("\n=== TEST 05: nostabilizesplits ===\n")
d5 <- make_data(seed = 42)
cf5 <- multi_arm_causal_forest(d5$X, d5$Y, d5$W_factor,
                                num.trees = 500, seed = 42,
                                stabilize.splits = FALSE)
save_macf("05_nostabilizesplits", d5, cf5)

# ─────────────────────────────────────────────────────────────────────────────
# TEST 06 — User-supplied Y.hat
# ─────────────────────────────────────────────────────────────────────────────
cat("\n=== TEST 06: User-supplied Y.hat ===\n")
d6 <- make_data(seed = 42)
set.seed(42)
rf_y6 <- regression_forest(d6$X, d6$Y, num.trees = 500, seed = 42)
yhat6 <- predict(rf_y6)$predictions
cf6 <- multi_arm_causal_forest(d6$X, d6$Y, d6$W_factor,
                                num.trees = 500, seed = 42,
                                Y.hat = yhat6)
save_macf("06_yhat_supplied", d6, cf6,
          extra_cols = list(yhat = yhat6))

# ─────────────────────────────────────────────────────────────────────────────
# TEST 07 — User-supplied W.hat (propensity for each arm)
# W.hat must be a matrix of dimension n x K for K treatment arms
# ─────────────────────────────────────────────────────────────────────────────
cat("\n=== TEST 07: User-supplied W.hat ===\n")
d7 <- make_data(seed = 42)
set.seed(42)
# W.hat must be n x (K+1) with colnames = levels(W_factor) = c("0","1","2")
rf_w7_0 <- regression_forest(d7$X, as.numeric(d7$W == 0), num.trees = 500, seed = 42)
rf_w7_1 <- regression_forest(d7$X, as.numeric(d7$W == 1), num.trees = 500, seed = 42)
rf_w7_2 <- regression_forest(d7$X, as.numeric(d7$W == 2), num.trees = 500, seed = 42)
what7_0 <- predict(rf_w7_0)$predictions
what7_1 <- predict(rf_w7_1)$predictions
what7_2 <- predict(rf_w7_2)$predictions
what7_mat <- cbind(what7_0, what7_1, what7_2)
colnames(what7_mat) <- levels(d7$W_factor)   # c("0","1","2")
cf7 <- multi_arm_causal_forest(d7$X, d7$Y, d7$W_factor,
                                num.trees = 500, seed = 42,
                                W.hat = what7_mat)
save_macf("07_what_supplied", d7, cf7,
          extra_cols = list(what0 = what7_0, what1 = what7_1, what2 = what7_2))

# ─────────────────────────────────────────────────────────────────────────────
# TEST 08 — clusters()
# ─────────────────────────────────────────────────────────────────────────────
cat("\n=== TEST 08: clusters() ===\n")
d8 <- make_data(seed = 42)
clusters8 <- rep(1:60, each = 10)  # 60 clusters of 10
cf8 <- multi_arm_causal_forest(d8$X, d8$Y, d8$W_factor,
                                num.trees = 500, seed = 42,
                                clusters = clusters8)
preds8_arr <- predict(cf8)$predictions
t8_1 <- get_arm(preds8_arr, 1); t8_2 <- get_arm(preds8_arr, 2)
df8 <- as.data.frame(d8$X)
df8$y <- d8$Y; df8$w <- d8$W
df8$w1 <- as.integer(d8$W == 1); df8$w2 <- as.integer(d8$W == 2)
df8$cluster_id <- clusters8
df8$tau_r_t1 <- t8_1; df8$tau_r_t2 <- t8_2
df8$tau_true_t1 <- d8$tau1; df8$tau_true_t2 <- d8$tau2
write.csv(df8, file.path(WORKDIR, "test08_cluster.csv"), row.names = FALSE)
cat(sprintf("  Saved test08_cluster.csv  (cor_t1=%.4f  cor_t2=%.4f)\n",
            cor(t8_1, d8$tau1), cor(t8_2, d8$tau2)))

# ─────────────────────────────────────────────────────────────────────────────
# TEST 09 — weights()
# ─────────────────────────────────────────────────────────────────────────────
cat("\n=== TEST 09: weights() ===\n")
d9 <- make_data(seed = 42)
set.seed(99)
wts9 <- runif(600, 0.5, 2.0)
cf9 <- multi_arm_causal_forest(d9$X, d9$Y, d9$W_factor,
                                num.trees = 500, seed = 42,
                                sample.weights = wts9)
preds9_arr <- predict(cf9)$predictions
t9_1 <- get_arm(preds9_arr, 1); t9_2 <- get_arm(preds9_arr, 2)
df9 <- as.data.frame(d9$X)
df9$y <- d9$Y; df9$w <- d9$W
df9$w1 <- as.integer(d9$W == 1); df9$w2 <- as.integer(d9$W == 2)
df9$obs_weight <- wts9
df9$tau_r_t1 <- t9_1; df9$tau_r_t2 <- t9_2
df9$tau_true_t1 <- d9$tau1; df9$tau_true_t2 <- d9$tau2
write.csv(df9, file.path(WORKDIR, "test09_weights.csv"), row.names = FALSE)
cat(sprintf("  Saved test09_weights.csv  (cor_t1=%.4f  cor_t2=%.4f)\n",
            cor(t9_1, d9$tau1), cor(t9_2, d9$tau2)))

# ─────────────────────────────────────────────────────────────────────────────
# TEST 10 — nohonesty (honesty=FALSE)
# ─────────────────────────────────────────────────────────────────────────────
cat("\n=== TEST 10: nohonesty ===\n")
d10 <- make_data(seed = 42)
cf10 <- multi_arm_causal_forest(d10$X, d10$Y, d10$W_factor,
                                 num.trees = 500, seed = 42,
                                 honesty = FALSE)
save_macf("10_nohonesty", d10, cf10)

# ─────────────────────────────────────────────────────────────────────────────
# TEST 11 — mtry=2
# ─────────────────────────────────────────────────────────────────────────────
cat("\n=== TEST 11: mtry=2 ===\n")
d11 <- make_data(seed = 42)
cf11 <- multi_arm_causal_forest(d11$X, d11$Y, d11$W_factor,
                                 num.trees = 500, seed = 42,
                                 mtry = 2)
save_macf("11_mtry2", d11, cf11)

# ─────────────────────────────────────────────────────────────────────────────
# TEST 12 — min.node.size=20
# ─────────────────────────────────────────────────────────────────────────────
cat("\n=== TEST 12: min.node.size=20 ===\n")
d12 <- make_data(seed = 42)
cf12 <- multi_arm_causal_forest(d12$X, d12$Y, d12$W_factor,
                                 num.trees = 500, seed = 42,
                                 min.node.size = 20)
save_macf("12_minnodesize20", d12, cf12)

# ─────────────────────────────────────────────────────────────────────────────
# TEST 13 — Homogeneous effects (tau1=2, tau2=-1)
# ─────────────────────────────────────────────────────────────────────────────
cat("\n=== TEST 13: Homogeneous effects ===\n")
set.seed(42)
n13 <- 600; p13 <- 5
X13 <- matrix(rnorm(n13 * p13), n13, p13)
colnames(X13) <- paste0("x", 1:p13)
W13_int <- sample(0:2, n13, replace = TRUE)
W13_factor <- factor(W13_int)
tau13_1 <- rep(2, n13)
tau13_2 <- rep(-1, n13)
Y13 <- X13[, 1] + 2 * (W13_int == 1) + (-1) * (W13_int == 2) + rnorm(n13)

cf13 <- multi_arm_causal_forest(X13, Y13, W13_factor,
                                 num.trees = 500, seed = 42)
preds13_arr <- predict(cf13)$predictions
t13_1 <- get_arm(preds13_arr, 1); t13_2 <- get_arm(preds13_arr, 2)

df13 <- as.data.frame(X13)
df13$y <- Y13; df13$w <- W13_int
df13$w1 <- as.integer(W13_int == 1); df13$w2 <- as.integer(W13_int == 2)
df13$tau_r_t1 <- t13_1; df13$tau_r_t2 <- t13_2
df13$tau_true_t1 <- tau13_1; df13$tau_true_t2 <- tau13_2
write.csv(df13, file.path(WORKDIR, "test13_homogeneous.csv"), row.names = FALSE)
cat(sprintf("  Saved test13_homogeneous.csv\n  Mean_t1=%.4f (true=2)  Mean_t2=%.4f (true=-1)\n",
            mean(t13_1), mean(t13_2)))

# ─────────────────────────────────────────────────────────────────────────────
# TEST 14 — Strong heterogeneity
# ─────────────────────────────────────────────────────────────────────────────
cat("\n=== TEST 14: Strong heterogeneity ===\n")
set.seed(42)
n14 <- 600; p14 <- 5
X14 <- matrix(rnorm(n14 * p14), n14, p14)
colnames(X14) <- paste0("x", 1:p14)
W14_int <- sample(0:2, n14, replace = TRUE)
W14_factor <- factor(W14_int)
tau14_1 <- 4 * (X14[, 1] + X14[, 2])
tau14_2 <- -4 * X14[, 1] + 6 * X14[, 3]
Y14 <- X14[, 1] + tau14_1 * (W14_int == 1) + tau14_2 * (W14_int == 2) + rnorm(n14)

cf14 <- multi_arm_causal_forest(X14, Y14, W14_factor,
                                 num.trees = 500, seed = 42)
preds14_arr <- predict(cf14)$predictions
t14_1 <- get_arm(preds14_arr, 1); t14_2 <- get_arm(preds14_arr, 2)

df14 <- as.data.frame(X14)
df14$y <- Y14; df14$w <- W14_int
df14$w1 <- as.integer(W14_int == 1); df14$w2 <- as.integer(W14_int == 2)
df14$tau_r_t1 <- t14_1; df14$tau_r_t2 <- t14_2
df14$tau_true_t1 <- tau14_1; df14$tau_true_t2 <- tau14_2
write.csv(df14, file.path(WORKDIR, "test14_strong_het.csv"), row.names = FALSE)
cat(sprintf("  Saved test14_strong_het.csv  (cor_t1=%.4f  cor_t2=%.4f)\n",
            cor(t14_1, tau14_1), cor(t14_2, tau14_2)))

# ─────────────────────────────────────────────────────────────────────────────
# TEST 15 — nuisancetrees=100 (both Y.hat and W.hat from 100-tree forests)
# ─────────────────────────────────────────────────────────────────────────────
cat("\n=== TEST 15: nuisancetrees=100 ===\n")
d15 <- make_data(seed = 42)
set.seed(42)
# W.hat must be n x (K+1) with colnames = levels(W_factor) = c("0","1","2")
rf_y15   <- regression_forest(d15$X, d15$Y, num.trees = 100, seed = 42)
rf_w15_0 <- regression_forest(d15$X, as.numeric(d15$W == 0), num.trees = 100, seed = 42)
rf_w15_1 <- regression_forest(d15$X, as.numeric(d15$W == 1), num.trees = 100, seed = 42)
rf_w15_2 <- regression_forest(d15$X, as.numeric(d15$W == 2), num.trees = 100, seed = 42)
yhat15   <- predict(rf_y15)$predictions
what15_0 <- predict(rf_w15_0)$predictions
what15_1 <- predict(rf_w15_1)$predictions
what15_2 <- predict(rf_w15_2)$predictions
what15_mat <- cbind(what15_0, what15_1, what15_2)
colnames(what15_mat) <- levels(d15$W_factor)

cf15 <- multi_arm_causal_forest(d15$X, d15$Y, d15$W_factor,
                                 num.trees = 500, seed = 42,
                                 Y.hat = yhat15, W.hat = what15_mat)
preds15_arr <- predict(cf15)$predictions
t15_1 <- get_arm(preds15_arr, 1); t15_2 <- get_arm(preds15_arr, 2)

df15 <- as.data.frame(d15$X)
df15$y <- d15$Y; df15$w <- d15$W
df15$w1 <- as.integer(d15$W == 1); df15$w2 <- as.integer(d15$W == 2)
df15$yhat  <- yhat15; df15$what1 <- what15_1; df15$what2 <- what15_2
df15$tau_r_t1 <- t15_1; df15$tau_r_t2 <- t15_2
df15$tau_true_t1 <- d15$tau1; df15$tau_true_t2 <- d15$tau2
write.csv(df15, file.path(WORKDIR, "test15_nuisancetrees100.csv"), row.names = FALSE)
cat(sprintf("  Saved test15_nuisancetrees100.csv  (cor_t1=%.4f  cor_t2=%.4f)\n",
            cor(t15_1, d15$tau1), cor(t15_2, d15$tau2)))

# ─────────────────────────────────────────────────────────────────────────────
# TEST 16 — estimate.variance (variance estimates for each arm)
# ci.group.size >= 2 required for variance estimation
# ─────────────────────────────────────────────────────────────────────────────
cat("\n=== TEST 16: estimate.variance ===\n")
d16 <- make_data(seed = 42)
cf16 <- multi_arm_causal_forest(d16$X, d16$Y, d16$W_factor,
                                 num.trees = 500, seed = 42,
                                 ci.group.size = 2)
preds16 <- predict(cf16, estimate.variance = TRUE)
pred16_arr <- preds16$predictions          # n x 2 x 1
pred16_var <- preds16$variance.estimates   # n x 2 matrix
t16_1 <- get_arm(pred16_arr, 1); t16_2 <- get_arm(pred16_arr, 2)
v16_1 <- pred16_var[, 1]; v16_2 <- pred16_var[, 2]

df16 <- as.data.frame(d16$X)
df16$y <- d16$Y; df16$w <- d16$W
df16$w1 <- as.integer(d16$W == 1); df16$w2 <- as.integer(d16$W == 2)
df16$tau_r_t1  <- t16_1; df16$tau_r_t2  <- t16_2
df16$var_r_t1  <- v16_1; df16$var_r_t2  <- v16_2
df16$tau_true_t1 <- d16$tau1; df16$tau_true_t2 <- d16$tau2
write.csv(df16, file.path(WORKDIR, "test16_variance.csv"), row.names = FALSE)
cat(sprintf("  Saved test16_variance.csv\n  (cor_t1=%.4f  cor_t2=%.4f  mean_var_t1=%.6f  mean_var_t2=%.6f)\n",
            cor(t16_1, d16$tau1), cor(t16_2, d16$tau2),
            mean(v16_1), mean(v16_2)))

cat("\n\n=== All R tests complete. CSVs written to", WORKDIR, "===\n\n")
