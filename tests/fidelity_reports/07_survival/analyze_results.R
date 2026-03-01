## Analyze R vs Stata fidelity results for survival_forest tests
## Computes per-column Pearson correlations and summary statistics

OUTDIR <- "/tmp/grf_stata/tests/fidelity_reports/07_survival"

## Helper: read stata CSV and compute per-column correlations
analyze_test <- function(testname, ncols, prefix_r = "r_s", prefix_s = "pred_s",
                         corr_threshold = 0.85) {
  fpath <- file.path(OUTDIR, paste0(testname, "_stata.csv"))
  if (!file.exists(fpath)) {
    cat("MISSING:", fpath, "\n")
    return(NULL)
  }
  df <- read.csv(fpath)

  r_cols  <- paste0(prefix_r, 1:ncols)
  s_cols  <- paste0(prefix_s, 1:ncols)

  # Check all expected columns exist
  missing_r <- setdiff(r_cols, names(df))
  missing_s <- setdiff(s_cols, names(df))
  if (length(missing_r) > 0 || length(missing_s) > 0) {
    cat("MISSING COLUMNS in", testname, ":\n  R:", missing_r, "\n  S:", missing_s, "\n")
    # Reduce to available columns
    available <- intersect(1:ncols, which(r_cols %in% names(df) & s_cols %in% names(df)))
    if (length(available) == 0) return(NULL)
    ncols <- max(available)
    r_cols <- r_cols[available]
    s_cols <- s_cols[available]
  }

  corrs <- numeric(ncols)
  for (j in seq_along(r_cols)) {
    rv <- df[[r_cols[j]]]
    sv <- df[[s_cols[j]]]
    ok <- !is.na(rv) & !is.na(sv)
    if (sum(ok) < 3) { corrs[j] <- NA; next }
    corrs[j] <- cor(rv[ok], sv[ok])
  }

  min_c <- min(corrs, na.rm = TRUE)
  mean_c <- mean(corrs, na.rm = TRUE)
  n_pass <- sum(corrs >= corr_threshold, na.rm = TRUE)
  n_fail <- sum(corrs < corr_threshold, na.rm = TRUE)

  list(
    test      = testname,
    ncols     = ncols,
    corrs     = corrs,
    min_corr  = min_c,
    mean_corr = mean_c,
    n_pass    = n_pass,
    n_fail    = n_fail,
    pass      = (n_fail == 0),
    threshold = corr_threshold
  )
}

## Helper: check monotonically decreasing survival curves
check_monotone <- function(testname, ncols, prefix_s = "pred_s") {
  fpath <- file.path(OUTDIR, paste0(testname, "_stata.csv"))
  if (!file.exists(fpath)) return(NULL)
  df <- read.csv(fpath)
  s_cols <- paste0(prefix_s, 1:ncols)
  avail <- intersect(s_cols, names(df))
  if (length(avail) < 2) return(TRUE)

  mat <- as.matrix(df[, avail])
  # For each row, check non-increasing
  n_nonmono <- 0
  for (i in 1:nrow(mat)) {
    row <- mat[i, ]
    row <- row[!is.na(row)]
    if (length(row) < 2) next
    if (any(diff(row) > 1e-10)) n_nonmono <- n_nonmono + 1
  }
  list(n_nonmono = n_nonmono, n_total = nrow(mat),
       pct_mono = 100 * (1 - n_nonmono / nrow(mat)))
}

## ---- Run analyses ----
results <- list()

# Tests 01-17: survival curve comparisons (predtype pairs)
results[["test01"]] <- analyze_test("test01", 20)
results[["test02"]] <- analyze_test("test02", 50)
results[["test03"]] <- analyze_test("test03", 100)
results[["test04"]] <- analyze_test("test04", 20)
results[["test05"]] <- analyze_test("test05", 20)
results[["test06"]] <- analyze_test("test06", 20)
results[["test07"]] <- analyze_test("test07", 20)
results[["test08"]] <- analyze_test("test08", 20)
results[["test09"]] <- analyze_test("test09", 20)
results[["test10"]] <- analyze_test("test10", 20)
results[["test11"]] <- analyze_test("test11", 20)
results[["test12"]] <- analyze_test("test12", 20)
results[["test13"]] <- analyze_test("test13", 20)
results[["test16"]] <- analyze_test("test16", 20)
results[["test17"]] <- analyze_test("test17", 50)

## Tests 14-15: expected survival
analyze_esurv <- function(testname, corr_threshold = 0.90) {
  fpath <- file.path(OUTDIR, paste0(testname, "_stata.csv"))
  if (!file.exists(fpath)) { cat("MISSING:", fpath, "\n"); return(NULL) }
  df <- read.csv(fpath)
  if (!all(c("stata_esurv", "r_esurv") %in% names(df))) {
    cat("Missing esurv columns in", testname, "\n")
    return(NULL)
  }
  ok <- !is.na(df$stata_esurv) & !is.na(df$r_esurv)
  corr <- cor(df$stata_esurv[ok], df$r_esurv[ok])
  list(test = testname, corr = corr,
       pass = corr >= corr_threshold, threshold = corr_threshold,
       stata_mean = mean(df$stata_esurv[ok]),
       r_mean     = mean(df$r_esurv[ok]),
       stata_sd   = sd(df$stata_esurv[ok]),
       r_sd       = sd(df$r_esurv[ok]))
}

# Expected survival (for test14 we need special column names from stata CSV)
esurv14 <- analyze_esurv("test14", 0.90)
esurv15 <- analyze_esurv("test15", 0.90)

## Test 15 consistency: check stata_esurv negatively correlated with x1
test15_df <- read.csv(file.path(OUTDIR, "test15_stata.csv"))
cor_x1_stata <- NA
cor_x1_r     <- NA
if ("stata_esurv" %in% names(test15_df) && "x1" %in% names(test15_df)) {
  ok15 <- !is.na(test15_df$stata_esurv) & !is.na(test15_df$x1)
  cor_x1_stata <- cor(test15_df$stata_esurv[ok15], test15_df$x1[ok15])
}
if ("r_esurv" %in% names(test15_df) && "x1" %in% names(test15_df)) {
  ok15r <- !is.na(test15_df$r_esurv) & !is.na(test15_df$x1)
  cor_x1_r <- cor(test15_df$r_esurv[ok15r], test15_df$x1[ok15r])
}

## Monotone checks (on Stata output)
mono <- list()
for (t in c("test01","test02","test03","test04","test05","test06","test07",
            "test08","test09","test10","test11","test12","test13","test16")) {
  ncl <- switch(t, "test02"=50, "test03"=100, 20)
  mono[[t]] <- check_monotone(t, ncl)
}
mono[["test17"]] <- check_monotone("test17", 50)

## ---- Print results summary ----
cat("\n")
cat("=============================================================\n")
cat("  R vs Stata Fidelity: survival_forest / grf_expected_survival\n")
cat("=============================================================\n\n")

cat(sprintf("%-12s %-6s %-8s %-8s %-8s %-8s %-8s\n",
            "Test", "NCols", "MinCorr", "MeanCorr", "NPass", "NFail", "Status"))
cat(strrep("-", 70), "\n")

test_labels <- list(
  test01 = "01 Default (noutput=20, KM)",
  test02 = "02 noutput=50",
  test03 = "03 noutput=100",
  test04 = "04 Nelson-Aalen",
  test05 = "05 KM explicit",
  test06 = "06 nofastlogrank",
  test07 = "07 cluster()",
  test08 = "08 weights()",
  test09 = "09 nohonesty",
  test10 = "10 mtry=2",
  test11 = "11 minnodesize=20",
  test12 = "12 High event rate",
  test13 = "13 Low event rate",
  test16 = "16 numfailures(50)",
  test17 = "17 Combined"
)

for (tn in names(test_labels)) {
  res <- results[[tn]]
  if (is.null(res)) {
    cat(sprintf("%-30s  ERROR: no results\n", test_labels[[tn]]))
    next
  }
  status <- if (res$pass) "PASS" else "FAIL"
  cat(sprintf("%-30s  cols=%3d  min=%.4f  mean=%.4f  pass=%d  fail=%d  [%s]\n",
              test_labels[[tn]], res$ncols, res$min_corr, res$mean_corr,
              res$n_pass, res$n_fail, status))
}

cat("\n")
cat("--- Expected Survival ---\n")
for (es in list(esurv14, esurv15)) {
  if (is.null(es)) { cat("MISSING\n"); next }
  lab <- if (es$test == "test14") "14 E[T|X] base" else "15 E[T|X] consistency"
  status <- if (es$pass) "PASS" else "FAIL"
  cat(sprintf("%-30s  corr=%.4f  (threshold=%.2f)  [%s]\n",
              lab, es$corr, es$threshold, status))
  cat(sprintf("  R mean=%.4f  Stata mean=%.4f\n", es$r_mean, es$stata_mean))
}

cat("\n")
cat("--- E[T|X] Directional Consistency (Test 15) ---\n")
cat(sprintf("  cor(r_esurv,    x1) = %.4f  (should be < 0)\n", cor_x1_r))
cat(sprintf("  cor(stata_esurv,x1) = %.4f  (should be < 0)\n", cor_x1_stata))

cat("\n")
cat("--- Monotone Survival Curves (Stata output) ---\n")
for (tn in names(mono)) {
  m <- mono[[tn]]
  if (is.null(m)) next
  status <- if (m$n_nonmono == 0) "PASS" else "FAIL"
  cat(sprintf("  %-20s  mono=%.1f%%  nonmono=%d/%d  [%s]\n",
              tn, m$pct_mono, m$n_nonmono, m$n_total, status))
}

## ---- Save results as RData for report generation ----
save(results, esurv14, esurv15, cor_x1_r, cor_x1_stata,
     mono, test_labels,
     file = file.path(OUTDIR, "analysis_results.RData"))
cat("\nResults saved to analysis_results.RData\n")

## ---- Per-column correlation detail for key tests ----
cat("\n=== Per-column correlation detail ===\n")
for (tn in c("test01","test02","test03","test04","test05","test06","test07","test08")) {
  res <- results[[tn]]
  if (is.null(res)) next
  cat(sprintf("\n%s (n=%d cols):\n", test_labels[[tn]], res$ncols))
  # Show first 10 and last 5 columns
  show_idx <- unique(c(1:min(10, res$ncols),
                       max(1, res$ncols-4):res$ncols))
  for (j in show_idx) {
    c_str <- if (!is.na(res$corrs[j])) sprintf("%.4f", res$corrs[j]) else "  NA  "
    status <- if (!is.na(res$corrs[j]) && res$corrs[j] >= res$threshold) "ok" else "!!"
    cat(sprintf("  col %3d: r=%.4f [%s]\n", j, res$corrs[j], status))
  }
}
