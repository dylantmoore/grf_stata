## Compute correlations between R and Stata predictions for all 18 tests
## Returns a data.frame of results for inclusion in the report
library(grf)

OUTDIR <- "/tmp/grf_stata/tests/fidelity_reports/11_ll_regression"

tests <- list(
    list(id =  1, label = "Default LL (all vars, lambda=0.1)",
         file = "test01_stata.csv", r_col = "r_pred", s_col = "stata_pred",
         threshold = 0.90, note = ""),
    list(id =  2, label = "llenable (enable_ll_split=TRUE)",
         file = "test02_stata.csv", r_col = "r_pred", s_col = "stata_pred",
         threshold = 0.90, note = ""),
    list(id =  3, label = "llvars(x1 x2) — subset LL correction",
         file = "test03_stata.csv", r_col = "r_pred", s_col = "stata_pred",
         threshold = 0.90, note = ""),
    list(id =  4, label = "llsplitvars(x1) — restricted split vars",
         file = "test04_stata.csv", r_col = "r_pred", s_col = "stata_pred",
         threshold = 0.90, note = ""),
    list(id =  5, label = "lllambda=0.01 — small regularization",
         file = "test05_stata.csv", r_col = "r_pred", s_col = "stata_pred",
         threshold = 0.90, note = ""),
    list(id =  6, label = "lllambda=1.0 — large regularization",
         file = "test06_stata.csv", r_col = "r_pred", s_col = "stata_pred",
         threshold = 0.90, note = ""),
    list(id =  7, label = "lllambda=10.0 — very large regularization",
         file = "test07_stata.csv", r_col = "r_pred", s_col = "stata_pred",
         threshold = 0.90, note = ""),
    list(id =  8, label = "llweightpenalty",
         file = "test08_stata.csv", r_col = "r_pred", s_col = "stata_pred",
         threshold = 0.90, note = ""),
    list(id =  9, label = "llcutoff=3",
         file = "test09_stata.csv", r_col = "r_pred", s_col = "stata_pred",
         threshold = 0.90, note = ""),
    list(id = 10, label = "cluster()",
         file = "test10_stata.csv", r_col = "r_pred", s_col = "stata_pred",
         threshold = 0.90, note = ""),
    list(id = 11, label = "weights() — R lacks sample.weights in ll_regression_forest",
         file = "test11_stata.csv", r_col = "r_pred", s_col = "stata_pred",
         threshold = 0.80,
         note = "PARTIAL_MATCH: R ll_regression_forest has no sample.weights param; R uses unweighted predictions"),
    list(id = 12, label = "nohonesty",
         file = "test12_stata.csv", r_col = "r_pred", s_col = "stata_pred",
         threshold = 0.90, note = ""),
    list(id = 13, label = "mtry=2",
         file = "test13_stata.csv", r_col = "r_pred", s_col = "stata_pred",
         threshold = 0.90, note = ""),
    list(id = 14, label = "minnodesize=20",
         file = "test14_stata.csv", r_col = "r_pred", s_col = "stata_pred",
         threshold = 0.90, note = ""),
    list(id = 15, label = "Combined: llvars(x1 x2) + lllambda=0.5 + llweightpenalty",
         file = "test15_stata.csv", r_col = "r_pred", s_col = "stata_pred",
         threshold = 0.90, note = ""),
    list(id = 16, label = "LL vs non-LL comparison (linear DGP)",
         file = "test16_stata.csv", r_col = "r_pred", s_col = "stata_pred_ll",
         threshold = 0.90, note = "Extra: MSE comparison LL vs std RF"),
    list(id = 17, label = "Nonlinear DGP (Y = sin(X1) + X2^2)",
         file = "test17_stata.csv", r_col = "r_pred", s_col = "stata_pred",
         threshold = 0.90, note = ""),
    list(id = 18, label = "Pure Linear DGP (Y = sum(X)) — LL should excel",
         file = "test18_stata.csv", r_col = "r_pred", s_col = "stata_pred",
         threshold = 0.90, note = "Extra: LL MSE vs std RF MSE")
)

results <- data.frame(
    test_id   = integer(),
    label     = character(),
    corr      = numeric(),
    threshold = numeric(),
    pass      = logical(),
    r_mean    = numeric(),
    r_sd      = numeric(),
    s_mean    = numeric(),
    s_sd      = numeric(),
    note      = character(),
    stringsAsFactors = FALSE
)

cat("=" , strrep("=", 60), "\n")
cat("   R-vs-Stata Fidelity: ll_regression_forest (18 tests)\n")
cat("=" , strrep("=", 60), "\n\n")

for (t in tests) {
    fp <- file.path(OUTDIR, t$file)
    if (!file.exists(fp)) {
        cat(sprintf("  Test %02d [MISSING]: %s\n", t$id, t$label))
        next
    }
    df <- read.csv(fp)
    r_v <- df[[t$r_col]]
    s_v <- df[[t$s_col]]
    corr <- cor(r_v, s_v, use = "complete.obs")
    pass <- corr >= t$threshold

    results <- rbind(results, data.frame(
        test_id   = t$id,
        label     = t$label,
        corr      = round(corr, 4),
        threshold = t$threshold,
        pass      = pass,
        r_mean    = round(mean(r_v, na.rm = TRUE), 4),
        r_sd      = round(sd(r_v, na.rm = TRUE), 4),
        s_mean    = round(mean(s_v, na.rm = TRUE), 4),
        s_sd      = round(sd(s_v, na.rm = TRUE), 4),
        note      = t$note,
        stringsAsFactors = FALSE
    ))
    status <- if (pass) "PASS" else "FAIL"
    cat(sprintf("  Test %02d [%s] corr=%.4f (threshold=%.2f): %s\n",
                t$id, status, corr, t$threshold, t$label))
    if (t$note != "") cat(sprintf("         NOTE: %s\n", t$note))
}

## Extra comparisons for Tests 16 and 18
cat("\n--- Extra MSE Comparisons ---\n")

## Test 16: LL vs non-LL for linear DGP
df16 <- read.csv(file.path(OUTDIR, "test16_stata.csv"))
r_mse_ll   <- mean((df16$r_pred - df16$y)^2, na.rm = TRUE)
r_mse_noll <- mean((df16$r_pred_noll - df16$y)^2, na.rm = TRUE)
r_mse_std  <- mean((df16$r_pred_std - df16$y)^2, na.rm = TRUE)
s_mse_ll   <- mean((df16$stata_pred_ll - df16$y)^2, na.rm = TRUE)
s_mse_noll <- mean((df16$stata_pred_noll - df16$y)^2, na.rm = TRUE)
cat("\nTest 16: LL vs non-LL (linear+interaction DGP):\n")
cat(sprintf("  R   MSE LL  (linear.correction.variables=1:5, lambda=0.1): %.4f\n", r_mse_ll))
cat(sprintf("  R   MSE noLL (plain regression forest predict):              %.4f\n", r_mse_noll))
cat(sprintf("  R   MSE std forest:                                          %.4f\n", r_mse_std))
cat(sprintf("  Stata MSE LL:                                                %.4f\n", s_mse_ll))
cat(sprintf("  Stata MSE noLL:                                              %.4f\n", s_mse_noll))
cat(sprintf("  R LL lower than R noLL:     %s\n", r_mse_ll < r_mse_noll))
cat(sprintf("  Stata LL lower than no-LL:  %s\n", s_mse_ll < s_mse_noll))

## Test 18: Pure linear DGP
df18 <- read.csv(file.path(OUTDIR, "test18_stata.csv"))
r_mse_ll18  <- mean((df18$r_pred - df18$y)^2, na.rm = TRUE)
r_mse_std18 <- mean((df18$r_pred_std - df18$y)^2, na.rm = TRUE)
s_mse_ll18  <- mean((df18$stata_pred - df18$y)^2, na.rm = TRUE)
s_mse_std18 <- mean((df18$stata_pred_std - df18$y)^2, na.rm = TRUE)
cat("\nTest 18: Pure Linear DGP (Y = sum(X)):\n")
cat(sprintf("  R   MSE LL:  %.4f  Std RF: %.4f  LL better: %s\n",
            r_mse_ll18, r_mse_std18, r_mse_ll18 < r_mse_std18))
cat(sprintf("  Stata MSE LL: %.4f  Std RF: %.4f  LL better: %s\n",
            s_mse_ll18, s_mse_std18, s_mse_ll18 < s_mse_std18))

## Summary
n_pass <- sum(results$pass)
n_total <- nrow(results)
cat(sprintf("\n\n=== SUMMARY: %d / %d tests PASSED (threshold >= corr required) ===\n",
            n_pass, n_total))

write.csv(results, file.path(OUTDIR, "correlation_results.csv"), row.names = FALSE)
cat("Results saved to: correlation_results.csv\n")

## Return results invisibly for use in report
invisible(results)
