## Master R script: Run all 18 ll_regression_forest fidelity tests
## R uses ll_regression_forest() from grf package
## DGP: Y = 3*X1 + 2*X2 - X3 + 0.5*X1*X2 + N(0,1), n=500, p=5
library(grf)

OUTDIR <- "/tmp/grf_stata/tests/fidelity_reports/11_ll_regression"

dgp_standard <- function(seed_val = 42, n = 500, p = 5) {
    set.seed(seed_val)
    X <- matrix(rnorm(n * p), n, p)
    Y <- 3*X[,1] + 2*X[,2] - X[,3] + 0.5*X[,1]*X[,2] + rnorm(n)
    list(X = X, Y = Y)
}

dgp_nonlinear <- function(seed_val = 42, n = 500, p = 5) {
    set.seed(seed_val)
    X <- matrix(rnorm(n * p), n, p)
    Y <- sin(X[,1]) + X[,2]^2 + rnorm(n)
    list(X = X, Y = Y)
}

dgp_linear_pure <- function(seed_val = 42, n = 500, p = 5) {
    set.seed(seed_val)
    X <- matrix(rnorm(n * p), n, p)
    Y <- rowSums(X) + rnorm(n)
    list(X = X, Y = Y)
}

save_csv <- function(X, Y, preds, filename, extra = list()) {
    df <- as.data.frame(X)
    colnames(df) <- paste0("x", 1:ncol(X))
    df$y <- Y
    df$r_pred <- preds
    for (nm in names(extra)) df[[nm]] <- extra[[nm]]
    write.csv(df, filename, row.names = FALSE)
    cat("  Saved:", filename, "\n")
    cat("  Predictions - mean:", round(mean(preds), 4),
        " sd:", round(sd(preds), 4),
        " range:", round(range(preds), 4), "\n")
    invisible(df)
}

## ============================================================
## Test 01: Default LL (no ll splits, default lambda=0.1)
## Stata: grf_ll_regression_forest y x1-x5, gen(pred) ntrees(500) seed(42)
## R: ll_regression_forest then predict with linear.correction.variables=1:5
## ============================================================
cat("\n=== Test 01: Default LL (all vars, lambda=0.1) ===\n")
d <- dgp_standard()
set.seed(42)
rf01 <- ll_regression_forest(d$X, d$Y, num.trees = 500, seed = 42)
preds01 <- predict(rf01, linear.correction.variables = 1:5, ll.lambda = 0.1)$predictions
save_csv(d$X, d$Y, preds01, file.path(OUTDIR, "test01_data.csv"))

## ============================================================
## Test 02: llenable (enable_ll_split=TRUE, all vars)
## Stata: ... llenable
## R: ll_regression_forest(enable.ll.split=TRUE)
## ============================================================
cat("\n=== Test 02: llenable (enable_ll_split=TRUE) ===\n")
d <- dgp_standard()
set.seed(42)
rf02 <- ll_regression_forest(d$X, d$Y, num.trees = 500, seed = 42,
                               enable.ll.split = TRUE)
preds02 <- predict(rf02, linear.correction.variables = 1:5, ll.lambda = 0.1)$predictions
save_csv(d$X, d$Y, preds02, file.path(OUTDIR, "test02_data.csv"))

## ============================================================
## Test 03: llvars(x1 x2) — subset of variables for LL correction
## Stata: ... llvars(x1 x2)
## R: predict(..., linear.correction.variables = c(1,2))
## ============================================================
cat("\n=== Test 03: llvars(x1 x2) — subset LL correction ===\n")
d <- dgp_standard()
set.seed(42)
rf03 <- ll_regression_forest(d$X, d$Y, num.trees = 500, seed = 42)
preds03 <- predict(rf03, linear.correction.variables = c(1, 2), ll.lambda = 0.1)$predictions
save_csv(d$X, d$Y, preds03, file.path(OUTDIR, "test03_data.csv"))

## ============================================================
## Test 04: llsplitvars(x1) — LL splits restricted to x1 only
## Stata: ... llenable llvars(x1)  (llsplitvars maps to ll.split.variables in R)
## R: ll_regression_forest(enable.ll.split=TRUE, ll.split.variables=1)
## Note: In Stata, llvars controls BOTH split vars AND correction vars when enable.ll.split set.
## We compare: Stata llsplitvars(x1) => R ll.split.variables=1 with predict on all vars
## ============================================================
cat("\n=== Test 04: llsplitvars(x1) — restricted split vars ===\n")
d <- dgp_standard()
set.seed(42)
rf04 <- ll_regression_forest(d$X, d$Y, num.trees = 500, seed = 42,
                               enable.ll.split = TRUE,
                               ll.split.variables = 1)
preds04 <- predict(rf04, linear.correction.variables = 1:5, ll.lambda = 0.1)$predictions
save_csv(d$X, d$Y, preds04, file.path(OUTDIR, "test04_data.csv"))

## ============================================================
## Test 05: lllambda=0.01 — small regularization
## Stata: ... lllambda(0.01)
## R: predict(..., ll.lambda = 0.01)
## ============================================================
cat("\n=== Test 05: lllambda=0.01 — small regularization ===\n")
d <- dgp_standard()
set.seed(42)
rf05 <- ll_regression_forest(d$X, d$Y, num.trees = 500, seed = 42)
preds05 <- predict(rf05, linear.correction.variables = 1:5, ll.lambda = 0.01)$predictions
save_csv(d$X, d$Y, preds05, file.path(OUTDIR, "test05_data.csv"))

## ============================================================
## Test 06: lllambda=1.0 — large regularization
## Stata: ... lllambda(1.0)
## R: predict(..., ll.lambda = 1.0)
## ============================================================
cat("\n=== Test 06: lllambda=1.0 — large regularization ===\n")
d <- dgp_standard()
set.seed(42)
rf06 <- ll_regression_forest(d$X, d$Y, num.trees = 500, seed = 42)
preds06 <- predict(rf06, linear.correction.variables = 1:5, ll.lambda = 1.0)$predictions
save_csv(d$X, d$Y, preds06, file.path(OUTDIR, "test06_data.csv"))

## ============================================================
## Test 07: lllambda=10.0 — very large regularization (near standard forest)
## Stata: ... lllambda(10.0)
## R: predict(..., ll.lambda = 10.0)
## ============================================================
cat("\n=== Test 07: lllambda=10.0 — very large (near standard forest) ===\n")
d <- dgp_standard()
set.seed(42)
rf07 <- ll_regression_forest(d$X, d$Y, num.trees = 500, seed = 42)
preds07 <- predict(rf07, linear.correction.variables = 1:5, ll.lambda = 10.0)$predictions
## Also get standard (no LL) for comparison
preds07_noll <- predict(rf07)$predictions
save_csv(d$X, d$Y, preds07, file.path(OUTDIR, "test07_data.csv"),
         extra = list(r_pred_noll = preds07_noll))
cat("  Corr LL(10.0) vs no-LL:", round(cor(preds07, preds07_noll), 4), "\n")

## ============================================================
## Test 08: llweightpenalty — enable weight penalty (ll.weight.penalty=TRUE)
## Stata: ... llweightpenalty
## R: predict(..., ll.weight.penalty = TRUE)
## ============================================================
cat("\n=== Test 08: llweightpenalty ===\n")
d <- dgp_standard()
set.seed(42)
rf08 <- ll_regression_forest(d$X, d$Y, num.trees = 500, seed = 42)
preds08 <- predict(rf08, linear.correction.variables = 1:5,
                   ll.lambda = 0.1, ll.weight.penalty = TRUE)$predictions
save_csv(d$X, d$Y, preds08, file.path(OUTDIR, "test08_data.csv"))

## ============================================================
## Test 09: llcutoff=3 — small cutoff (use leaf betas for small leaves)
## Stata: ... llcutoff(3)
## R: ll_regression_forest(ll.split.cutoff = 3)
## Note: ll_split_cutoff affects the SPLIT phase (tree building), not prediction
## ============================================================
cat("\n=== Test 09: llcutoff=3 ===\n")
d <- dgp_standard()
set.seed(42)
rf09 <- ll_regression_forest(d$X, d$Y, num.trees = 500, seed = 42,
                               enable.ll.split = TRUE,
                               ll.split.cutoff = 3)
preds09 <- predict(rf09, linear.correction.variables = 1:5, ll.lambda = 0.1)$predictions
save_csv(d$X, d$Y, preds09, file.path(OUTDIR, "test09_data.csv"))

## ============================================================
## Test 10: cluster() — with cluster variable
## Stata: ... cluster(clust)
## R: ll_regression_forest(clusters = clust)
## ============================================================
cat("\n=== Test 10: cluster() ===\n")
d <- dgp_standard()
set.seed(42)
n <- nrow(d$X)
clust <- rep(1:50, each = 10)  # 50 clusters of size 10
rf10 <- ll_regression_forest(d$X, d$Y, num.trees = 500, seed = 42,
                              clusters = clust)
preds10 <- predict(rf10, linear.correction.variables = 1:5, ll.lambda = 0.1)$predictions
df10 <- as.data.frame(d$X)
colnames(df10) <- paste0("x", 1:5)
df10$y <- d$Y
df10$r_pred <- preds10
df10$clust <- clust
write.csv(df10, file.path(OUTDIR, "test10_data.csv"), row.names = FALSE)
cat("  Saved: test10_data.csv\n")
cat("  Predictions mean:", round(mean(preds10), 4), " sd:", round(sd(preds10), 4), "\n")

## ============================================================
## Test 11: weights() — with observation weights
## NOTE: grf 2.5.0 ll_regression_forest has NO sample.weights parameter.
## Stata wrapper (grf_ll_regression_forest) accepts weights() via forest-level
## weighting in the C++ plugin. R closest equivalent: ll_regression_forest
## with no weights. We document this API difference; compare Stata weighted
## predictions to R unweighted for qualitative verification, then flag as
## PARTIAL_MATCH (expect lower correlation since weighting differs).
## ============================================================
cat("\n=== Test 11: weights() — NOTE: R ll_regression_forest has no sample.weights ===\n")
d <- dgp_standard()
set.seed(42)
n <- nrow(d$X)
wts <- abs(rnorm(n)) + 0.1
## R: ll_regression_forest without weights (closest available)
rf11 <- ll_regression_forest(d$X, d$Y, num.trees = 500, seed = 42)
preds11 <- predict(rf11, linear.correction.variables = 1:5, ll.lambda = 0.1)$predictions
df11 <- as.data.frame(d$X)
colnames(df11) <- paste0("x", 1:5)
df11$y <- d$Y
df11$r_pred <- preds11   ## R without weights (reference — will differ from Stata weighted)
df11$wt <- wts
write.csv(df11, file.path(OUTDIR, "test11_data.csv"), row.names = FALSE)
cat("  Saved: test11_data.csv\n")
cat("  R predictions (unweighted): mean =", round(mean(preds11), 4),
    " sd =", round(sd(preds11), 4), "\n")
cat("  NOTE: Stata will use weights — expect lower correlation vs Test 01\n")

## ============================================================
## Test 12: nohonesty — without honesty
## Stata: ... nohonesty
## R: ll_regression_forest(honesty = FALSE)
## ============================================================
cat("\n=== Test 12: nohonesty ===\n")
d <- dgp_standard()
set.seed(42)
rf12 <- ll_regression_forest(d$X, d$Y, num.trees = 500, seed = 42,
                               honesty = FALSE)
preds12 <- predict(rf12, linear.correction.variables = 1:5, ll.lambda = 0.1)$predictions
save_csv(d$X, d$Y, preds12, file.path(OUTDIR, "test12_data.csv"))

## ============================================================
## Test 13: mtry=2 — restricted splitting variable count
## Stata: ... mtry(2)
## R: ll_regression_forest(mtry = 2)
## ============================================================
cat("\n=== Test 13: mtry=2 ===\n")
d <- dgp_standard()
set.seed(42)
rf13 <- ll_regression_forest(d$X, d$Y, num.trees = 500, seed = 42,
                               mtry = 2)
preds13 <- predict(rf13, linear.correction.variables = 1:5, ll.lambda = 0.1)$predictions
save_csv(d$X, d$Y, preds13, file.path(OUTDIR, "test13_data.csv"))

## ============================================================
## Test 14: minnodesize=20 — larger leaves
## Stata: ... minnodesize(20)
## R: ll_regression_forest(min.node.size = 20)
## ============================================================
cat("\n=== Test 14: minnodesize=20 ===\n")
d <- dgp_standard()
set.seed(42)
rf14 <- ll_regression_forest(d$X, d$Y, num.trees = 500, seed = 42,
                               min.node.size = 20)
preds14 <- predict(rf14, linear.correction.variables = 1:5, ll.lambda = 0.1)$predictions
save_csv(d$X, d$Y, preds14, file.path(OUTDIR, "test14_data.csv"))

## ============================================================
## Test 15: Combined — llvars(x1 x2) + lllambda=0.5 + llweightpenalty
## Stata: ... llvars(x1 x2) lllambda(0.5) llweightpenalty
## R: predict(..., linear.correction.variables = c(1,2), ll.lambda = 0.5, ll.weight.penalty = TRUE)
## ============================================================
cat("\n=== Test 15: Combined llvars(x1 x2) + lllambda=0.5 + llweightpenalty ===\n")
d <- dgp_standard()
set.seed(42)
rf15 <- ll_regression_forest(d$X, d$Y, num.trees = 500, seed = 42)
preds15 <- predict(rf15, linear.correction.variables = c(1, 2),
                   ll.lambda = 0.5, ll.weight.penalty = TRUE)$predictions
save_csv(d$X, d$Y, preds15, file.path(OUTDIR, "test15_data.csv"))

## ============================================================
## Test 16: Compare LL vs non-LL (standard regression forest)
## Both on linear DGP to show LL improvement
## Stata: grf_ll_regression_forest vs grf_regression_forest
## R: predict with/without linear.correction.variables
## ============================================================
cat("\n=== Test 16: LL vs non-LL (linear DGP) ===\n")
d <- dgp_standard()
set.seed(42)
rf16 <- ll_regression_forest(d$X, d$Y, num.trees = 500, seed = 42)
preds16_ll   <- predict(rf16, linear.correction.variables = 1:5, ll.lambda = 0.1)$predictions
preds16_noll <- predict(rf16)$predictions
## Standard regression forest for comparison
set.seed(42)
rf16_std <- regression_forest(d$X, d$Y, num.trees = 500, seed = 42)
preds16_std <- predict(rf16_std)$predictions
true_mu16 <- 3*d$X[,1] + 2*d$X[,2] - d$X[,3] + 0.5*d$X[,1]*d$X[,2]
mse_ll   <- mean((preds16_ll   - d$Y)^2)
mse_noll <- mean((preds16_noll - d$Y)^2)
mse_std  <- mean((preds16_std  - d$Y)^2)
cat("  MSE LL:       ", round(mse_ll, 4), "\n")
cat("  MSE no-LL:    ", round(mse_noll, 4), "\n")
cat("  MSE std RF:   ", round(mse_std, 4), "\n")
cat("  LL lower MSE than no-LL:", mse_ll < mse_noll, "\n")
df16 <- as.data.frame(d$X); colnames(df16) <- paste0("x", 1:5)
df16$y <- d$Y; df16$r_pred <- preds16_ll; df16$r_pred_noll <- preds16_noll
df16$r_pred_std <- preds16_std; df16$true_mu <- true_mu16
write.csv(df16, file.path(OUTDIR, "test16_data.csv"), row.names = FALSE)
cat("  Saved: test16_data.csv\n")

## ============================================================
## Test 17: Nonlinear DGP — Y = sin(X1) + X2^2
## Stata: grf_ll_regression_forest y x1-x5, gen(pred) ntrees(500) seed(42)
## R: ll_regression_forest + predict with LL
## ============================================================
cat("\n=== Test 17: Nonlinear DGP (sin+quadratic) ===\n")
d <- dgp_nonlinear()
set.seed(42)
rf17 <- ll_regression_forest(d$X, d$Y, num.trees = 500, seed = 42)
preds17 <- predict(rf17, linear.correction.variables = 1:5, ll.lambda = 0.1)$predictions
save_csv(d$X, d$Y, preds17, file.path(OUTDIR, "test17_data.csv"))

## ============================================================
## Test 18: Pure Linear DGP — Y = sum(X) — LL should excel
## Stata: grf_ll_regression_forest y x1-x5, gen(pred) ntrees(500) seed(42)
## R: ll_regression_forest + predict with LL
## ============================================================
cat("\n=== Test 18: Pure Linear DGP (Y = sum(X)) — LL excels ===\n")
d <- dgp_linear_pure()
set.seed(42)
rf18 <- ll_regression_forest(d$X, d$Y, num.trees = 500, seed = 42)
preds18_ll <- predict(rf18, linear.correction.variables = 1:5, ll.lambda = 0.1)$predictions
## Standard forest for comparison
set.seed(42)
rf18_std <- regression_forest(d$X, d$Y, num.trees = 500, seed = 42)
preds18_std <- predict(rf18_std)$predictions
mse_ll18  <- mean((preds18_ll  - d$Y)^2)
mse_std18 <- mean((preds18_std - d$Y)^2)
cat("  MSE LL (pure linear):", round(mse_ll18, 4), "\n")
cat("  MSE std RF:          ", round(mse_std18, 4), "\n")
cat("  LL lower MSE:         ", mse_ll18 < mse_std18, "\n")
df18 <- as.data.frame(d$X); colnames(df18) <- paste0("x", 1:5)
df18$y <- d$Y; df18$r_pred <- preds18_ll; df18$r_pred_std <- preds18_std
write.csv(df18, file.path(OUTDIR, "test18_data.csv"), row.names = FALSE)
cat("  Saved: test18_data.csv\n")

cat("\n\n=== All 18 R tests complete ===\n")
cat("Output directory:", OUTDIR, "\n")
