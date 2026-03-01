## Master comparison script for regression_forest fidelity tests
## Computes Pearson correlations between R and Stata predictions

results <- list()

safe_cor <- function(x, y) {
  valid <- is.finite(x) & is.finite(y)
  if (sum(valid) < 2) return(NA)
  cor(x[valid], y[valid])
}

compare_test <- function(test_id, test_name, stata_file,
                         r_col = "r_pred", s_col = "stata_pred",
                         var_r_col = NULL, var_s_col = NULL,
                         extra_notes = "") {
  cat(sprintf("\n=== Test %s: %s ===\n", test_id, test_name))
  if (!file.exists(stata_file)) {
    cat("  MISSING file:", stata_file, "\n")
    return(list(id = test_id, name = test_name, corr = NA, result = "MISSING", notes = extra_notes))
  }
  df <- read.csv(stata_file)
  cat("  Rows:", nrow(df), "| Cols:", paste(names(df), collapse=", "), "\n")

  if (!r_col %in% names(df) || !s_col %in% names(df)) {
    cat("  ERROR: Missing expected columns\n")
    return(list(id = test_id, name = test_name, corr = NA, result = "ERROR", notes = extra_notes))
  }

  cc <- safe_cor(df[[r_col]], df[[s_col]])
  cat(sprintf("  Prediction correlation (R vs Stata): %.6f\n", cc))
  cat(sprintf("  R preds: mean=%.4f sd=%.4f\n", mean(df[[r_col]], na.rm=TRUE), sd(df[[r_col]], na.rm=TRUE)))
  cat(sprintf("  Stata preds: mean=%.4f sd=%.4f\n", mean(df[[s_col]], na.rm=TRUE), sd(df[[s_col]], na.rm=TRUE)))

  result_str <- ifelse(is.na(cc), "FAIL",
                 ifelse(cc > 0.95, "STRONG PASS",
                   ifelse(cc > 0.90, "PASS", "FAIL")))
  cat("  Result:", result_str, "\n")

  var_corr <- NA
  var_result <- NA
  if (!is.null(var_r_col) && !is.null(var_s_col) && var_r_col %in% names(df) && var_s_col %in% names(df)) {
    var_corr <- safe_cor(df[[var_r_col]], df[[var_s_col]])
    cat(sprintf("  Variance correlation (R vs Stata): %.6f\n", var_corr))
    var_result <- ifelse(is.na(var_corr), "FAIL", ifelse(var_corr > 0.85, "PASS", "FAIL"))
    cat("  Variance result:", var_result, "\n")
  }

  if (!is.null(extra_notes) && nchar(extra_notes) > 0) cat("  Note:", extra_notes, "\n")

  list(id = test_id, name = test_name, corr = cc, result = result_str,
       var_corr = var_corr, var_result = var_result, notes = extra_notes)
}

base_dir <- "/tmp/grf_stata/tests/fidelity_reports/01_regression"

# Test 01: Default
results[[1]] <- compare_test("01", "Default options (n=500, p=5, ntrees=500, seed=42)",
  file.path(base_dir, "test01_stata.csv"))

# Test 02: nohonesty
results[[2]] <- compare_test("02", "nohonesty (honesty=FALSE)",
  file.path(base_dir, "test02_stata.csv"))

# Test 03: mtry=2
results[[3]] <- compare_test("03", "mtry=2",
  file.path(base_dir, "test03_stata.csv"))

# Test 04: minnodesize=20
results[[4]] <- compare_test("04", "min.node.size=20",
  file.path(base_dir, "test04_stata.csv"))

# Test 05: samplefrac=0.3
results[[5]] <- compare_test("05", "sample.fraction=0.3",
  file.path(base_dir, "test05_stata.csv"))

# Test 06: honestyfrac=0.7
results[[6]] <- compare_test("06", "honesty.fraction=0.7",
  file.path(base_dir, "test06_stata.csv"))

# Test 07: alpha=0.15
results[[7]] <- compare_test("07", "alpha=0.15",
  file.path(base_dir, "test07_stata.csv"))

# Test 08: imbalancepenalty=1.0
results[[8]] <- compare_test("08", "imbalance.penalty=1.0",
  file.path(base_dir, "test08_stata.csv"))

# Test 09: estimatevariance
results[[9]] <- compare_test("09", "estimate.variance (ci.group.size=2)",
  file.path(base_dir, "test09_stata.csv"),
  var_r_col = "r_var", var_s_col = "stata_var")

# Test 10: cluster
results[[10]] <- compare_test("10", "clusters (10 clusters)",
  file.path(base_dir, "test10_stata.csv"))

# Test 11: weights
results[[11]] <- compare_test("11", "sample.weights (random positive weights)",
  file.path(base_dir, "test11_stata.csv"))

# Test 12: equalizeclusterweights
results[[12]] <- compare_test("12", "equalize.cluster.weights=TRUE",
  file.path(base_dir, "test12_stata.csv"))

# Test 13: nomia
results[[13]] <- compare_test("13", "nomia (no MIA) — R default MIA always on; tested on complete data",
  file.path(base_dir, "test13_stata.csv"),
  extra_notes = "grf 2.5.0 has no enable.missing.indicator param; test uses complete data (no NAs). R uses MIA=on (default); Stata uses nomia. Results should still match on complete data.")

# Test 14: cigroupsize=2
results[[14]] <- compare_test("14", "ci.group.size=2",
  file.path(base_dir, "test14_stata.csv"))

# Test 15: combined nohonesty+mtry=3+minnodesize=15
results[[15]] <- compare_test("15", "Combined: nohonesty + mtry=3 + min.node.size=15",
  file.path(base_dir, "test15_stata.csv"))

# Test 16: cluster+weights+estimatevariance
results[[16]] <- compare_test("16", "Combined: cluster + weights + estimatevariance",
  file.path(base_dir, "test16_stata.csv"),
  var_r_col = "r_var", var_s_col = "stata_var")

# Test 17: large p=20
results[[17]] <- compare_test("17", "Large p (p=20)",
  file.path(base_dir, "test17_stata.csv"))

# Test 18a: ntrees=100
results[[18]] <- compare_test("18a", "ntrees=100",
  file.path(base_dir, "test18a_stata.csv"))

# Test 18b: ntrees=2000
results[[19]] <- compare_test("18b", "ntrees=2000",
  file.path(base_dir, "test18b_stata.csv"))

# Test 19: missing data
results[[20]] <- compare_test("19", "Missing data (10% NAs in x1, x3) with MIA on",
  file.path(base_dir, "test19_stata.csv"))

# Test 20: seed reproducibility
cat("\n=== Test 20: Seed reproducibility ===\n")
df20 <- read.csv(file.path(base_dir, "test20_stata.csv"))
cat("  Stata internal reproducibility (same seed, corr):",
    safe_cor(df20$stata_pred1, df20$stata_pred2), "\n")
cat("  Stata exact match:", all(df20$stata_pred1 == df20$stata_pred2, na.rm=TRUE), "\n")
stata_r_corr <- safe_cor(df20$r_pred, df20$stata_pred1)
cat(sprintf("  R vs Stata corr: %.6f\n", stata_r_corr))
r_check <- read.csv(file.path(base_dir, "test20_r_seed_check.csv"))
cat("  R internal reproducibility corr:", r_check$r_corr, "| exact:", r_check$exact_match == 1, "\n")
result_str20 <- ifelse(is.na(stata_r_corr), "FAIL",
                  ifelse(stata_r_corr > 0.95, "STRONG PASS",
                    ifelse(stata_r_corr > 0.90, "PASS", "FAIL")))
cat("  Result:", result_str20, "\n")
results[[21]] <- list(id = "20", name = "Seed reproducibility (seed=123)",
                      corr = stata_r_corr, result = result_str20,
                      var_corr = NA, var_result = NA,
                      notes = paste0("Stata internal corr=", safe_cor(df20$stata_pred1, df20$stata_pred2),
                                     "; exact match=", all(df20$stata_pred1 == df20$stata_pred2)))

# Summary
cat("\n\n========== SUMMARY ==========\n")
cat(sprintf("%-5s %-55s %10s %10s\n", "Test", "Name", "Pred-corr", "Result"))
cat(strrep("-", 85), "\n")
for (r in results) {
  varstr <- if (!is.null(r$var_corr) && !is.na(r$var_corr))
              sprintf(" | var-corr=%.4f [%s]", r$var_corr, r$var_result) else ""
  cat(sprintf("%-5s %-55s %10.6f %10s%s\n",
    r$id,
    substr(r$name, 1, 55),
    ifelse(is.na(r$corr), -999, r$corr),
    r$result,
    varstr))
}
n_pass <- sum(sapply(results, function(r) grepl("PASS", r$result)))
n_total <- length(results)
cat(sprintf("\nTotal PASS: %d / %d\n", n_pass, n_total))

# Save results to CSV for the report
out <- do.call(rbind, lapply(results, function(r) {
  data.frame(test = r$id, name = r$name,
             corr = ifelse(is.null(r$corr) || is.na(r$corr), NA, r$corr),
             result = r$result,
             var_corr = ifelse(is.null(r$var_corr) || is.na(r$var_corr), NA, r$var_corr),
             var_result = ifelse(is.null(r$var_result) || is.na(r$var_result), "", as.character(r$var_result)),
             notes = ifelse(is.null(r$notes) || is.na(r$notes), "", r$notes),
             stringsAsFactors = FALSE)
}))
write.csv(out, file.path(base_dir, "results_summary.csv"), row.names = FALSE)
cat("\nResults saved to results_summary.csv\n")
