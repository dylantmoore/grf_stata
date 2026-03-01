## ATE Fidelity Tests - R Reference
## All 19 tests for average_treatment_effect from grf

library(grf)

outdir <- "/tmp/grf_stata/tests/fidelity_reports/03_ate"
results <- list()

save_result <- function(test_id, test_name, ate, se,
                        expected_error=FALSE, r_errored=FALSE,
                        note="", true_ate=NA_real_) {
  results[[length(results)+1]] <<- list(
    test_id=test_id, test_name=test_name,
    ate=ate, se=se,
    expected_error=expected_error, r_errored=r_errored,
    note=note, true_ate=true_ate
  )
  if (is.na(ate)) {
    cat(sprintf("  R Test %02d [%s]: errored/NA (expected=%s)\n",
                test_id, test_name, expected_error))
  } else {
    cat(sprintf("  R Test %02d [%s]: ATE=%.6f SE=%.6f\n",
                test_id, test_name, ate, se))
  }
}

## ============================================================
## Helper: generate standard DGP
## ============================================================
gen_data <- function(n, p=5, seed=42, tau_fn=NULL, prob_treat=0.5) {
  set.seed(seed)
  X <- matrix(runif(n * p), n, p)
  colnames(X) <- paste0("x", 1:p)
  W <- rbinom(n, 1, prob_treat)
  if (is.null(tau_fn)) {
    tau <- rep(2, n)
  } else {
    tau <- tau_fn(X)
  }
  Y <- X[,1] + 2 * X[,2] + tau * W + rnorm(n, 0, 0.5)
  list(X=X, Y=Y, W=W, tau=tau)
}

## ============================================================
## SHARED DATA: n=500 default DGP (used in tests 1-4,7-10)
## ============================================================
d <- gen_data(500)
cat("Fitting main causal forest (n=500)...\n")
cf <- causal_forest(d$X, d$Y, d$W, num.trees=2000, seed=42)
df_main <- as.data.frame(d$X)
df_main$y <- d$Y
df_main$w <- d$W

## ============================================================
## TEST 1: ATE (all) — default AIPW
## ============================================================
cat("Running Test 01: ATE (all)...\n")
res <- average_treatment_effect(cf)
save_result(1, "ATE_all", res["estimate"], res["std.err"])
write.csv(df_main, file.path(outdir, "test01_data.csv"), row.names=FALSE)

## ============================================================
## TEST 2: ATT (treated)
## ============================================================
cat("Running Test 02: ATT (treated)...\n")
res <- average_treatment_effect(cf, target.sample="treated")
save_result(2, "ATT_treated", res["estimate"], res["std.err"])

## ============================================================
## TEST 3: ATC (control)
## ============================================================
cat("Running Test 03: ATC (control)...\n")
res <- average_treatment_effect(cf, target.sample="control")
save_result(3, "ATC_control", res["estimate"], res["std.err"])

## ============================================================
## TEST 4: Overlap weights
## ============================================================
cat("Running Test 04: Overlap...\n")
res <- average_treatment_effect(cf, target.sample="overlap")
save_result(4, "ATE_overlap", res["estimate"], res["std.err"])

## ============================================================
## TEST 5: Debiasing weights
## ============================================================
cat("Running Test 05: Debiasing weights...\n")
dbw <- abs(d$X[,1]) + 0.5
res <- average_treatment_effect(cf, debiasing.weights=dbw)
save_result(5, "ATE_debias", res["estimate"], res["std.err"])
df5 <- df_main
df5$dbw <- dbw
write.csv(df5, file.path(outdir, "test05_data.csv"), row.names=FALSE)

## ============================================================
## TEST 6: Large sample n=2000
## ============================================================
cat("Running Test 06: Large n=2000...\n")
d6 <- gen_data(2000, seed=43)
cat("Fitting causal forest n=2000...\n")
cf6 <- causal_forest(d6$X, d6$Y, d6$W, num.trees=2000, seed=42)
res <- average_treatment_effect(cf6)
save_result(6, "ATE_large_n", res["estimate"], res["std.err"])
df6 <- as.data.frame(d6$X)
df6$y <- d6$Y
df6$w <- d6$W
write.csv(df6, file.path(outdir, "test06_data.csv"), row.names=FALSE)

## ============================================================
## TEST 7: TMLE ATE
## ============================================================
cat("Running Test 07: TMLE ATE...\n")
res <- average_treatment_effect(cf, method="TMLE", target.sample="all")
save_result(7, "TMLE_ATE", res["estimate"], res["std.err"],
            note="TMLE method")

## ============================================================
## TEST 8: TMLE ATT
## ============================================================
cat("Running Test 08: TMLE ATT...\n")
res <- average_treatment_effect(cf, method="TMLE", target.sample="treated")
save_result(8, "TMLE_ATT", res["estimate"], res["std.err"],
            note="TMLE method")

## ============================================================
## TEST 9: TMLE ATC
## ============================================================
cat("Running Test 09: TMLE ATC...\n")
res <- average_treatment_effect(cf, method="TMLE", target.sample="control")
save_result(9, "TMLE_ATC", res["estimate"], res["std.err"],
            note="TMLE method")

## ============================================================
## TEST 10: TMLE + overlap → silently redirects to AIPW in grf 2.5
## ============================================================
cat("Running Test 10: TMLE + overlap redirect...\n")
note10 <- ""
res <- withCallingHandlers(
  average_treatment_effect(cf, method="TMLE", target.sample="overlap"),
  message=function(m) {
    note10 <<- conditionMessage(m)
    invokeRestart("muffleMessage")
  }
)
if (note10 == "") note10 <- "silent_redirect_to_AIPW_in_grf2.5"
save_result(10, "TMLE_overlap_redirect", res["estimate"], res["std.err"],
            note=note10)

## ============================================================
## TEST 11: TMLE + debiasing.weights
## In grf 2.5, R does NOT error — debiasing.weights is silently IGNORED
## when method=TMLE. Stata should error with rc=198.
## ============================================================
cat("Running Test 11: TMLE + debiasing.weights...\n")
res11 <- tryCatch(
  average_treatment_effect(cf, method="TMLE", debiasing.weights=dbw),
  error=function(e) { cat("  R errored:", conditionMessage(e), "\n"); NULL }
)
if (!is.null(res11)) {
  save_result(11, "TMLE_debias_error", res11["estimate"], res11["std.err"],
              expected_error=FALSE, r_errored=FALSE,
              note="R_ignores_debias_with_TMLE_grf2.5")
} else {
  save_result(11, "TMLE_debias_error", NA_real_, NA_real_,
              expected_error=TRUE, r_errored=TRUE,
              note="R_errored")
}

## ============================================================
## TEST 12: Constant treatment effect tau=2
## ============================================================
cat("Running Test 12: Constant tau=2...\n")
d12 <- gen_data(500, seed=100, tau_fn=function(X) rep(2, nrow(X)))
cat("Fitting causal forest tau=2...\n")
cf12 <- causal_forest(d12$X, d12$Y, d12$W, num.trees=2000, seed=42)
res <- average_treatment_effect(cf12)
save_result(12, "const_tau2", res["estimate"], res["std.err"], true_ate=2.0)
df12 <- as.data.frame(d12$X)
df12$y <- d12$Y
df12$w <- d12$W
write.csv(df12, file.path(outdir, "test12_data.csv"), row.names=FALSE)

## ============================================================
## TEST 13: Heterogeneous tau(X) = 3*X1
## ============================================================
cat("Running Test 13: Heterogeneous tau=3*X1...\n")
d13 <- gen_data(500, seed=101, tau_fn=function(X) 3 * X[,1])
cat("Fitting causal forest hetero tau...\n")
cf13 <- causal_forest(d13$X, d13$Y, d13$W, num.trees=2000, seed=42)
res <- average_treatment_effect(cf13)
true_ate13 <- mean(3 * d13$X[,1])
save_result(13, "hetero_tau", res["estimate"], res["std.err"],
            true_ate=true_ate13)
df13 <- as.data.frame(d13$X)
df13$y <- d13$Y
df13$w <- d13$W
write.csv(df13, file.path(outdir, "test13_data.csv"), row.names=FALSE)

## ============================================================
## TEST 14: Unbalanced treatment (80% treated)
## ============================================================
cat("Running Test 14: Unbalanced 80% treated...\n")
d14 <- gen_data(500, seed=102, prob_treat=0.8)
cat("Fitting causal forest unbalanced...\n")
cf14 <- causal_forest(d14$X, d14$Y, d14$W, num.trees=2000, seed=42)
res <- average_treatment_effect(cf14)
save_result(14, "unbalanced_80", res["estimate"], res["std.err"])
df14 <- as.data.frame(d14$X)
df14$y <- d14$Y
df14$w <- d14$W
write.csv(df14, file.path(outdir, "test14_data.csv"), row.names=FALSE)

## ============================================================
## TEST 15: With clusters (20 clusters)
## ============================================================
cat("Running Test 15: With clusters...\n")
set.seed(103)
n15 <- 500; p15 <- 5
X15 <- matrix(runif(n15 * p15), n15, p15)
colnames(X15) <- paste0("x", 1:p15)
clusters15 <- sample(1:20, n15, replace=TRUE)
W15 <- rbinom(n15, 1, 0.5)
Y15 <- X15[,1] + 2*X15[,2] + 2*W15 + rnorm(n15, 0, 0.5)
cat("Fitting causal forest with clusters...\n")
cf15 <- causal_forest(X15, Y15, W15, clusters=clusters15, num.trees=2000, seed=42)
res <- average_treatment_effect(cf15)
save_result(15, "cluster_ATE", res["estimate"], res["std.err"])
df15 <- as.data.frame(X15)
df15$y <- Y15
df15$w <- W15
df15$cluster_id <- clusters15
write.csv(df15, file.path(outdir, "test15_data.csv"), row.names=FALSE)

## ============================================================
## TEST 16: With clusters + TMLE → R errors "not implemented"
## ============================================================
cat("Running Test 16: Clusters + TMLE...\n")
res16 <- tryCatch(
  average_treatment_effect(cf15, method="TMLE"),
  error=function(e) { cat("  R errored:", conditionMessage(e), "\n"); NULL }
)
if (!is.null(res16)) {
  save_result(16, "cluster_TMLE", res16["estimate"], res16["std.err"])
} else {
  save_result(16, "cluster_TMLE", NA_real_, NA_real_,
              expected_error=TRUE, r_errored=TRUE,
              note="TMLE_not_implemented_with_clusters")
}

## ============================================================
## TEST 17: Clustered + equalize.cluster.weights → R errors "unused argument"
## ============================================================
cat("Running Test 17: Clusters + equalize weights...\n")
res17 <- tryCatch(
  average_treatment_effect(cf15, target.sample="overlap", equalize.cluster.weights=TRUE),
  error=function(e) { cat("  R errored:", conditionMessage(e), "\n"); NULL }
)
if (!is.null(res17)) {
  save_result(17, "cluster_equalize", res17["estimate"], res17["std.err"])
} else {
  # Fall back to plain overlap ATE with cluster forest
  res17b <- average_treatment_effect(cf15, target.sample="overlap")
  save_result(17, "cluster_equalize", res17b["estimate"], res17b["std.err"],
              note="equalize.cluster.weights_not_in_grf2.5_used_overlap")
}

## ============================================================
## TEST 18: Large treatment effect tau=10
## ============================================================
cat("Running Test 18: Large tau=10...\n")
d18 <- gen_data(500, seed=104, tau_fn=function(X) rep(10, nrow(X)))
cat("Fitting causal forest large tau...\n")
cf18 <- causal_forest(d18$X, d18$Y, d18$W, num.trees=2000, seed=42)
res <- average_treatment_effect(cf18)
save_result(18, "large_tau10", res["estimate"], res["std.err"], true_ate=10.0)
df18 <- as.data.frame(d18$X)
df18$y <- d18$Y
df18$w <- d18$W
write.csv(df18, file.path(outdir, "test18_data.csv"), row.names=FALSE)

## ============================================================
## TEST 19: Near-zero treatment effect tau=0.01
## ============================================================
cat("Running Test 19: Near-zero tau=0.01...\n")
d19 <- gen_data(500, seed=105, tau_fn=function(X) rep(0.01, nrow(X)))
cat("Fitting causal forest near-zero tau...\n")
cf19 <- causal_forest(d19$X, d19$Y, d19$W, num.trees=2000, seed=42)
res <- average_treatment_effect(cf19)
save_result(19, "nearzero_tau", res["estimate"], res["std.err"], true_ate=0.01)
df19 <- as.data.frame(d19$X)
df19$y <- d19$Y
df19$w <- d19$W
write.csv(df19, file.path(outdir, "test19_data.csv"), row.names=FALSE)

## ============================================================
## Save all R results as CSV
## ============================================================
cat("\nSaving R results...\n")
df_results <- do.call(rbind, lapply(results, function(r) {
  data.frame(
    test_id        = r$test_id,
    test_name      = r$test_name,
    ate            = ifelse(is.null(r$ate) || is.na(r$ate),
                            NA_real_, as.numeric(r$ate)),
    se             = ifelse(is.null(r$se) || is.na(r$se),
                            NA_real_, as.numeric(r$se)),
    true_ate       = ifelse(is.null(r$true_ate) || is.na(r$true_ate),
                            NA_real_, as.numeric(r$true_ate)),
    expected_error = ifelse(is.null(r$expected_error), FALSE, r$expected_error),
    r_errored      = ifelse(is.null(r$r_errored), FALSE, r$r_errored),
    note           = ifelse(is.null(r$note), "", r$note),
    stringsAsFactors=FALSE
  )
}))

write.csv(df_results,
          file.path(outdir, "r_results.csv"),
          row.names=FALSE)
cat("R results saved to r_results.csv\n")
cat("Done!\n")
