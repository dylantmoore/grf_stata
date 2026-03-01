## Comparison script for all probability_forest fidelity tests
## Computes Pearson correlations, class agreement rates, and probability sum checks

library(grf)

WORKDIR <- "/tmp/grf_stata/tests/fidelity_reports/06_probability"

cat("\n", strrep("=", 70), "\n")
cat("PROBABILITY FOREST FIDELITY: R vs STATA COMPARISON\n")
cat(strrep("=", 70), "\n\n")

results <- list()

# ============================================================
# Helper: compute correlation for each class
# ============================================================
compare_class_probs <- function(df, r_cols, s_cols, test_name) {
  nclasses <- length(r_cols)
  cors <- numeric(nclasses)
  for (c in seq_along(r_cols)) {
    r_vals <- df[[r_cols[c]]]
    s_vals <- df[[s_cols[c]]]
    # Check for constant columns (all predictions identical)
    if (sd(r_vals, na.rm=TRUE) < 1e-10 || sd(s_vals, na.rm=TRUE) < 1e-10) {
      cors[c] <- NA
    } else {
      cors[c] <- cor(r_vals, s_vals, use = "complete.obs")
    }
  }
  names(cors) <- paste0("class_", 0:(nclasses-1))

  # Compute argmax class agreement
  r_mat <- as.matrix(df[, r_cols])
  s_mat <- as.matrix(df[, s_cols])
  r_class <- apply(r_mat, 1, which.max) - 1
  s_class <- apply(s_mat, 1, which.max) - 1
  agree_rate <- mean(r_class == s_class, na.rm = TRUE)

  # Sum check
  r_sum <- rowSums(r_mat)
  s_sum <- rowSums(s_mat)
  r_max_dev <- max(abs(r_sum - 1))
  s_max_dev <- max(abs(s_sum - 1))

  list(
    test = test_name,
    nclasses = nclasses,
    cors = cors,
    min_cor = min(cors, na.rm = TRUE),
    mean_cor = mean(cors, na.rm = TRUE),
    agree_rate = agree_rate,
    r_sum_dev = r_max_dev,
    s_sum_dev = s_max_dev,
    pass_cor = all(cors >= 0.90, na.rm = TRUE),
    pass_agree = agree_rate >= 0.80,
    pass_rsum = r_max_dev < 0.01,
    pass_ssum = s_max_dev < 0.01
  )
}

print_result <- function(res) {
  status <- if (res$pass_cor && res$pass_agree) "[PASS]" else "[FAIL]"
  cat(sprintf("\n%s %s\n", status, res$test))
  cat(sprintf("  Classes: %d\n", res$nclasses))
  cat("  Per-class Pearson correlations:\n")
  for (nm in names(res$cors)) {
    v <- res$cors[[nm]]
    pass <- !is.na(v) && v >= 0.90
    cat(sprintf("    %-12s r = %6.4f  %s\n", nm, ifelse(is.na(v), NA, v),
                ifelse(is.na(v), "(NA)", ifelse(pass, "PASS", "FAIL"))))
  }
  cat(sprintf("  Min correlation:     %.4f  %s\n", res$min_cor,
              ifelse(res$pass_cor, "PASS", "FAIL")))
  cat(sprintf("  Class agreement:     %.4f  %s\n", res$agree_rate,
              ifelse(res$pass_agree, "PASS", "FAIL")))
  cat(sprintf("  R prob sum max dev:  %.2e  %s\n", res$r_sum_dev,
              ifelse(res$pass_rsum, "PASS", "FAIL")))
  cat(sprintf("  Stata prob sum dev:  %.2e  %s\n", res$s_sum_dev,
              ifelse(res$pass_ssum, "PASS", "FAIL")))
}

# ============================================================
# Test 01: Binary classification
# ============================================================
df01 <- read.csv(file.path(WORKDIR, "test01_stata.csv"))
res01 <- compare_class_probs(df01,
  r_cols = c("r_pred_c0", "r_pred_c1"),
  s_cols = c("stata_pred_c0", "stata_pred_c1"),
  test_name = "Test 01: Binary Classification (n=500)")
print_result(res01)
results[["test01"]] <- res01

# ============================================================
# Test 02: 3-class
# ============================================================
df02 <- read.csv(file.path(WORKDIR, "test02_stata.csv"))
res02 <- compare_class_probs(df02,
  r_cols = c("r_pred_c0", "r_pred_c1", "r_pred_c2"),
  s_cols = c("stata_pred_c0", "stata_pred_c1", "stata_pred_c2"),
  test_name = "Test 02: 3-Class Multinomial (n=500)")
print_result(res02)
results[["test02"]] <- res02

# ============================================================
# Test 03: 5-class
# ============================================================
df03 <- read.csv(file.path(WORKDIR, "test03_stata.csv"))
res03 <- compare_class_probs(df03,
  r_cols = paste0("r_pred_c", 0:4),
  s_cols = paste0("stata_pred_c", 0:4),
  test_name = "Test 03: 5-Class Multinomial (n=500)")
print_result(res03)
results[["test03"]] <- res03

# ============================================================
# Test 04: nclasses specified vs auto
# ============================================================
df04 <- read.csv(file.path(WORKDIR, "test04_stata.csv"))
cat("\n[INFO] Test 04: nclasses specified vs auto-detect\n")
# Compare spec vs auto (should be identical)
cor_spec_auto_c0 <- cor(df04$stata_spec_c0, df04$stata_auto_c0)
cor_spec_auto_c1 <- cor(df04$stata_spec_c1, df04$stata_auto_c1)
max_diff_c0 <- max(abs(df04$stata_spec_c0 - df04$stata_auto_c0))
max_diff_c1 <- max(abs(df04$stata_spec_c1 - df04$stata_auto_c1))
cat(sprintf("  spec vs auto c0: cor=%.6f, max_diff=%.2e\n", cor_spec_auto_c0, max_diff_c0))
cat(sprintf("  spec vs auto c1: cor=%.6f, max_diff=%.2e\n", cor_spec_auto_c1, max_diff_c1))
identical_results <- max_diff_c0 < 1e-10 && max_diff_c1 < 1e-10
cat(sprintf("  Identical predictions: %s\n", ifelse(identical_results, "YES", "NO")))
# Also compare spec vs R
res04 <- compare_class_probs(
  data.frame(r_pred_c0 = df04$r_pred_c0, r_pred_c1 = df04$r_pred_c1,
             stata_pred_c0 = df04$stata_spec_c0, stata_pred_c1 = df04$stata_spec_c1),
  r_cols = c("r_pred_c0", "r_pred_c1"),
  s_cols = c("stata_pred_c0", "stata_pred_c1"),
  test_name = "Test 04: nclasses(2) specified vs auto (binary, n=500)")
print_result(res04)
results[["test04"]] <- res04

# ============================================================
# Test 05: Unbalanced classes
# ============================================================
df05 <- read.csv(file.path(WORKDIR, "test05_stata.csv"))
res05 <- compare_class_probs(df05,
  r_cols = c("r_pred_c0", "r_pred_c1"),
  s_cols = c("stata_pred_c0", "stata_pred_c1"),
  test_name = "Test 05: Unbalanced Classes ~84/16 (n=500)")
print_result(res05)
results[["test05"]] <- res05

# ============================================================
# Test 06: Cluster
# ============================================================
df06 <- read.csv(file.path(WORKDIR, "test06_stata.csv"))
res06 <- compare_class_probs(df06,
  r_cols = c("r_pred_c0", "r_pred_c1"),
  s_cols = c("stata_pred_c0", "stata_pred_c1"),
  test_name = "Test 06: With cluster() - 50 clusters (n=500)")
print_result(res06)
results[["test06"]] <- res06

# ============================================================
# Test 07: Weights
# ============================================================
df07 <- read.csv(file.path(WORKDIR, "test07_stata.csv"))
res07 <- compare_class_probs(df07,
  r_cols = c("r_pred_c0", "r_pred_c1"),
  s_cols = c("stata_pred_c0", "stata_pred_c1"),
  test_name = "Test 07: With sample weights (n=500)")
print_result(res07)
results[["test07"]] <- res07

# ============================================================
# Test 08: nohonesty
# ============================================================
df08 <- read.csv(file.path(WORKDIR, "test08_stata.csv"))
res08 <- compare_class_probs(df08,
  r_cols = c("r_pred_c0", "r_pred_c1"),
  s_cols = c("stata_pred_c0", "stata_pred_c1"),
  test_name = "Test 08: nohonesty (honesty=FALSE, n=500)")
print_result(res08)
results[["test08"]] <- res08

# ============================================================
# Test 09: mtry=2
# ============================================================
df09 <- read.csv(file.path(WORKDIR, "test09_stata.csv"))
res09 <- compare_class_probs(df09,
  r_cols = c("r_pred_c0", "r_pred_c1"),
  s_cols = c("stata_pred_c0", "stata_pred_c1"),
  test_name = "Test 09: mtry=2 (restricted splitting, n=500)")
print_result(res09)
results[["test09"]] <- res09

# ============================================================
# Test 10: minnodesize=20
# ============================================================
df10 <- read.csv(file.path(WORKDIR, "test10_stata.csv"))
res10 <- compare_class_probs(df10,
  r_cols = c("r_pred_c0", "r_pred_c1"),
  s_cols = c("stata_pred_c0", "stata_pred_c1"),
  test_name = "Test 10: minnodesize=20 (larger leaves, n=500)")
print_result(res10)
results[["test10"]] <- res10

# ============================================================
# Test 11: Combined cluster + weights + nohonesty
# ============================================================
df11 <- read.csv(file.path(WORKDIR, "test11_stata.csv"))
res11 <- compare_class_probs(df11,
  r_cols = c("r_pred_c0", "r_pred_c1"),
  s_cols = c("stata_pred_c0", "stata_pred_c1"),
  test_name = "Test 11: Combined cluster+weights+nohonesty (n=500)")
print_result(res11)
results[["test11"]] <- res11

# ============================================================
# Test 12: Probability sum check (3-class)
# ============================================================
df12 <- read.csv(file.path(WORKDIR, "test12_stata.csv"))
cat("\n[INFO] Test 12: Probability sum check (3-class)\n")
# R sums
r_sums <- df12$r_sum
cat(sprintf("  R: max|sum-1| = %.2e, all < 0.01: %s\n",
            max(abs(r_sums - 1)), all(abs(r_sums - 1) < 0.01)))
# Stata sums
s_sums <- df12$stata_sum
cat(sprintf("  Stata: max|sum-1| = %.2e, all < 0.01: %s\n",
            max(abs(s_sums - 1)), all(abs(s_sums - 1) < 0.01)))
# Also compute correlations
res12 <- compare_class_probs(df12,
  r_cols = c("r_pred_c0", "r_pred_c1", "r_pred_c2"),
  s_cols = c("stata_pred_c0", "stata_pred_c1", "stata_pred_c2"),
  test_name = "Test 12: Probability sum = 1 (3-class, n=500)")
print_result(res12)
results[["test12"]] <- res12

# ============================================================
# Test 13: Large sample n=2000
# ============================================================
df13 <- read.csv(file.path(WORKDIR, "test13_stata.csv"))
res13 <- compare_class_probs(df13,
  r_cols = c("r_pred_c0", "r_pred_c1"),
  s_cols = c("stata_pred_c0", "stata_pred_c1"),
  test_name = "Test 13: Large sample n=2000 (binary)")
print_result(res13)
results[["test13"]] <- res13

# ============================================================
# Summary table
# ============================================================
cat("\n\n", strrep("=", 70), "\n")
cat("SUMMARY TABLE\n")
cat(strrep("=", 70), "\n\n")
cat(sprintf("%-8s %-42s %-8s %-8s %-8s %-8s\n",
            "Test", "Description", "MinCor", "Agree", "SumOK", "Result"))
cat(strrep("-", 90), "\n")

for (nm in names(results)) {
  res <- results[[nm]]
  ssum_ok <- res$pass_rsum && res$pass_ssum
  overall <- res$pass_cor && res$pass_agree && ssum_ok
  cat(sprintf("%-8s %-42s %-8.4f %-8.4f %-8s %-8s\n",
              nm,
              substr(res$test, 1, 42),
              ifelse(is.finite(res$min_cor), res$min_cor, NA),
              res$agree_rate,
              ifelse(ssum_ok, "PASS", "FAIL"),
              ifelse(overall, "PASS", "FAIL")))
}

n_pass <- sum(sapply(results, function(r) r$pass_cor && r$pass_agree))
n_total <- length(results)
cat(sprintf("\nOverall: %d / %d tests passed\n", n_pass, n_total))

# Save summary as CSV for report
summary_df <- do.call(rbind, lapply(names(results), function(nm) {
  r <- results[[nm]]
  data.frame(
    test = nm,
    description = r$test,
    nclasses = r$nclasses,
    min_cor = round(r$min_cor, 4),
    mean_cor = round(r$mean_cor, 4),
    agree_rate = round(r$agree_rate, 4),
    r_sum_dev = format(r$r_sum_dev, scientific = TRUE, digits = 3),
    s_sum_dev = format(r$s_sum_dev, scientific = TRUE, digits = 3),
    pass_cor = r$pass_cor,
    pass_agree = r$pass_agree,
    pass_sum = r$pass_rsum && r$pass_ssum,
    overall = r$pass_cor && r$pass_agree && r$pass_rsum && r$pass_ssum,
    stringsAsFactors = FALSE
  )
}))
write.csv(summary_df, file.path(WORKDIR, "summary.csv"), row.names = FALSE)
cat("\nSummary saved to summary.csv\n")
