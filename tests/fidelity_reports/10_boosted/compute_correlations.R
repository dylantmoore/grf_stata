## compute_correlations.R -- Compare R vs Stata predictions for all 16 tests
library(grf)

results <- list()

compute_stats <- function(test_id, test_name, stata_file, r_col="r_pred", stata_col="stata_pred") {
  if (!file.exists(stata_file)) {
    cat("WARNING: file not found:", stata_file, "\n")
    return(NULL)
  }
  df <- read.csv(stata_file)

  if (!(r_col %in% names(df)) || !(stata_col %in% names(df))) {
    cat("WARNING: columns missing in", stata_file, "\n")
    cat("  Available:", paste(names(df), collapse=", "), "\n")
    return(NULL)
  }

  r_pred    <- df[[r_col]]
  s_pred    <- df[[stata_col]]
  true_mu   <- if ("true_mu" %in% names(df)) df$true_mu else NULL

  # Filter to non-missing
  ok <- !is.na(r_pred) & !is.na(s_pred)
  r_pred <- r_pred[ok]; s_pred <- s_pred[ok]

  corr_rs   <- cor(r_pred, s_pred)
  rmse_rs   <- sqrt(mean((r_pred - s_pred)^2))

  r_vs_true <- if (!is.null(true_mu)) cor(r_pred, true_mu[ok]) else NA
  s_vs_true <- if (!is.null(true_mu)) cor(s_pred, true_mu[ok]) else NA

  cat(sprintf("Test %s (%s):\n", test_id, test_name))
  cat(sprintf("  R-Stata corr:   %.6f  (PASS: %s)\n", corr_rs, ifelse(corr_rs > 0.90, "YES", "NO")))
  cat(sprintf("  R-Stata RMSE:   %.6f\n", rmse_rs))
  if (!is.na(r_vs_true)) cat(sprintf("  R-vs-true:      %.6f\n", r_vs_true))
  if (!is.na(s_vs_true)) cat(sprintf("  Stata-vs-true:  %.6f\n", s_vs_true))
  cat("\n")

  list(test_id=test_id, test_name=test_name, corr=corr_rs, rmse=rmse_rs,
       r_vs_true=r_vs_true, s_vs_true=s_vs_true, n=sum(ok),
       pass=corr_rs > 0.90)
}

basedir <- "/tmp/grf_stata/tests/fidelity_reports/10_boosted"

cat("=================================================================\n")
cat("Boosted Regression Forest: R vs Stata Fidelity Comparison\n")
cat("=================================================================\n\n")

tests_config <- list(
  list("01", "Default (auto-tune steps)", "test01_stata.csv"),
  list("02", "boost.steps=1", "test02_stata.csv"),
  list("03", "boost.steps=3", "test03_stata.csv"),
  list("04", "boost.steps=5", "test04_stata.csv"),
  list("05", "boost.max.steps=10", "test05_stata.csv"),
  list("06", "boost.error.reduction=0.90", "test06_stata.csv"),
  list("07", "boost.error.reduction=0.99", "test07_stata.csv"),
  list("08", "boost.trees.tune=50", "test08_stata.csv"),
  list("09", "nostabilizesplits", "test09_stata.csv"),
  list("10", "cluster()", "test10_stata.csv"),
  list("11", "weights()", "test11_stata.csv"),
  list("12", "nohonesty", "test12_stata.csv"),
  list("13", "mtry=2", "test13_stata.csv"),
  list("15", "Linear data", "test15_stata.csv"),
  list("16", "Highly nonlinear data", "test16_stata.csv")
)

all_results <- list()
for (tc in tests_config) {
  res <- compute_stats(tc[[1]], tc[[2]], file.path(basedir, tc[[3]]))
  if (!is.null(res)) all_results[[length(all_results)+1]] <- res
}

# Test 14 special: compare BRF vs RF
cat("Test 14 (Boosted vs Plain RF):\n")
df14 <- read.csv(file.path(basedir, "test14_stata.csv"))
cat(sprintf("  BRF: Stata-vs-R corr: %.6f  (PASS: %s)\n",
    cor(df14$brf_pred, df14$r_pred_brf), ifelse(cor(df14$brf_pred, df14$r_pred_brf) > 0.90, "YES", "NO")))
cat(sprintf("  RF:  Stata-vs-R corr: %.6f  (PASS: %s)\n",
    cor(df14$rf_pred, df14$r_pred_rf), ifelse(cor(df14$rf_pred, df14$r_pred_rf) > 0.90, "YES", "NO")))
cat(sprintf("  Boosted MSE (vs true): %.6f\n", mean(df14$brf_sq_err)))
cat(sprintf("  Plain RF MSE (vs true): %.6f\n", mean(df14$rf_sq_err)))
improvement <- (mean(df14$rf_sq_err) - mean(df14$brf_sq_err)) / mean(df14$rf_sq_err) * 100
cat(sprintf("  Boosted improvement over RF: %.2f%%\n\n", improvement))

brf14_corr <- cor(df14$brf_pred, df14$r_pred_brf)
rf14_corr  <- cor(df14$rf_pred,  df14$r_pred_rf)
all_results[[length(all_results)+1]] <- list(test_id="14a", test_name="BRF vs R (test14)", corr=brf14_corr, pass=brf14_corr>0.90)
all_results[[length(all_results)+1]] <- list(test_id="14b", test_name="RF vs R (test14)",  corr=rf14_corr,  pass=rf14_corr>0.90)

# Summary
cat("=================================================================\n")
cat("SUMMARY\n")
cat("=================================================================\n")
n_pass <- sum(sapply(all_results, function(r) r$pass))
n_total <- length(all_results)
cat(sprintf("Tests: %d PASS, %d FAIL, %d total\n\n", n_pass, n_total - n_pass, n_total))

for (r in all_results) {
  status <- if (r$pass) "PASS" else "FAIL"
  cat(sprintf("  [%s] Test %s: %s (corr=%.4f)\n", status, r$test_id, r$test_name, r$corr))
}

# Save to RDS for report generation
saveRDS(all_results, file.path(basedir, "correlation_results.rds"))
cat("\nResults saved to correlation_results.rds\n")
