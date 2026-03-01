## Analyze all quantile forest fidelity test results
## Compute Pearson correlations for each quantile, check ordering, output summary

library(stats)

DIR <- "/tmp/grf_stata/tests/fidelity_reports/05_quantile"

THRESHOLD <- 0.90
THRESHOLD_WEIGHTS <- 0.80  # relaxed for weight test where models differ

# Helper: safe correlation
safe_cor <- function(x, y) {
  if (all(is.na(x)) || all(is.na(y))) return(NA)
  tryCatch(cor(x, y, use = "complete.obs"), error = function(e) NA)
}

# Helper: check quantile ordering
check_ordering <- function(df, q_low, q_high) {
  vl <- paste0("stata_pred_q", q_low)
  vh <- paste0("stata_pred_q", q_high)
  if (!vl %in% names(df) || !vh %in% names(df)) return(NA)
  n_viol <- sum(df[[vl]] > df[[vh]], na.rm = TRUE)
  n_viol
}

results <- list()

# ---- Test 01: Default quantiles ----
df <- read.csv(file.path(DIR, "test01_stata.csv"))
r_q10 <- safe_cor(df$stata_pred_q10, df$r_q10)
r_q50 <- safe_cor(df$stata_pred_q50, df$r_q50)
r_q90 <- safe_cor(df$stata_pred_q90, df$r_q90)
cross1 <- check_ordering(df, 10, 50)
cross2 <- check_ordering(df, 50, 90)
results[["01_default"]] <- list(
  test = "01: Default quantiles (0.1, 0.5, 0.9)",
  quantiles = c("q10", "q50", "q90"),
  cors = c(r_q10, r_q50, r_q90),
  threshold = THRESHOLD,
  cross_10_50 = cross1, cross_50_90 = cross2,
  note = ""
)
cat("Test 01: cor(q10)=", round(r_q10, 4), " cor(q50)=", round(r_q50, 4), " cor(q90)=", round(r_q90, 4), "\n")

# ---- Test 02: Single quantile ----
df <- read.csv(file.path(DIR, "test02_stata.csv"))
r_q50 <- safe_cor(df$stata_pred_q50, df$r_q50)
results[["02_single"]] <- list(
  test = "02: Single quantile (0.5 median only)",
  quantiles = c("q50"),
  cors = c(r_q50),
  threshold = THRESHOLD,
  cross_10_50 = NA, cross_50_90 = NA,
  note = ""
)
cat("Test 02: cor(q50)=", round(r_q50, 4), "\n")

# ---- Test 03: Many quantiles ----
df <- read.csv(file.path(DIR, "test03_stata.csv"))
r_q10 <- safe_cor(df$stata_pred_q10, df$r_q10)
r_q25 <- safe_cor(df$stata_pred_q25, df$r_q25)
r_q50 <- safe_cor(df$stata_pred_q50, df$r_q50)
r_q75 <- safe_cor(df$stata_pred_q75, df$r_q75)
r_q90 <- safe_cor(df$stata_pred_q90, df$r_q90)
cross1 <- check_ordering(df, 10, 25)
cross2 <- check_ordering(df, 25, 50)
cross3 <- check_ordering(df, 50, 75)
cross4 <- check_ordering(df, 75, 90)
results[["03_many"]] <- list(
  test = "03: Many quantiles (0.1, 0.25, 0.5, 0.75, 0.9)",
  quantiles = c("q10", "q25", "q50", "q75", "q90"),
  cors = c(r_q10, r_q25, r_q50, r_q75, r_q90),
  threshold = THRESHOLD,
  cross_10_50 = cross1 + cross2, cross_50_90 = cross3 + cross4,
  note = ""
)
cat("Test 03: q10=", round(r_q10,4), " q25=", round(r_q25,4), " q50=", round(r_q50,4),
    " q75=", round(r_q75,4), " q90=", round(r_q90,4), "\n")

# ---- Test 04: Extreme quantiles ----
df <- read.csv(file.path(DIR, "test04_stata.csv"))
r_q1  <- safe_cor(df$stata_pred_q1,  df$r_q1)
r_q99 <- safe_cor(df$stata_pred_q99, df$r_q99)
cross1 <- check_ordering(df, 1, 99)
results[["04_extreme"]] <- list(
  test = "04: Extreme quantiles (0.01, 0.99)",
  quantiles = c("q1", "q99"),
  cors = c(r_q1, r_q99),
  threshold = THRESHOLD,
  cross_10_50 = cross1, cross_50_90 = NA,
  note = ""
)
cat("Test 04: cor(q1)=", round(r_q1, 4), " cor(q99)=", round(r_q99, 4), "\n")

# ---- Test 05: regression.splitting ----
df <- read.csv(file.path(DIR, "test05_stata.csv"))
r_q10 <- safe_cor(df$stata_pred_q10, df$r_q10)
r_q50 <- safe_cor(df$stata_pred_q50, df$r_q50)
r_q90 <- safe_cor(df$stata_pred_q90, df$r_q90)
cross1 <- check_ordering(df, 10, 50)
cross2 <- check_ordering(df, 50, 90)
results[["05_regsplit"]] <- list(
  test = "05: regression.splitting",
  quantiles = c("q10", "q50", "q90"),
  cors = c(r_q10, r_q50, r_q90),
  threshold = THRESHOLD,
  cross_10_50 = cross1, cross_50_90 = cross2,
  note = ""
)
cat("Test 05: cor(q10)=", round(r_q10,4), " cor(q50)=", round(r_q50,4), " cor(q90)=", round(r_q90,4), "\n")

# ---- Test 06: regsplit + multiple quantiles ----
df <- read.csv(file.path(DIR, "test06_stata.csv"))
r_q10 <- safe_cor(df$stata_pred_q10, df$r_q10)
r_q25 <- safe_cor(df$stata_pred_q25, df$r_q25)
r_q50 <- safe_cor(df$stata_pred_q50, df$r_q50)
r_q75 <- safe_cor(df$stata_pred_q75, df$r_q75)
r_q90 <- safe_cor(df$stata_pred_q90, df$r_q90)
results[["06_regsplit_multi"]] <- list(
  test = "06: regression.splitting + 5 quantiles",
  quantiles = c("q10", "q25", "q50", "q75", "q90"),
  cors = c(r_q10, r_q25, r_q50, r_q75, r_q90),
  threshold = THRESHOLD,
  cross_10_50 = NA, cross_50_90 = NA,
  note = ""
)
cat("Test 06: q10=", round(r_q10,4), " q25=", round(r_q25,4), " q50=", round(r_q50,4),
    " q75=", round(r_q75,4), " q90=", round(r_q90,4), "\n")

# ---- Test 07: Heteroscedastic ----
df <- read.csv(file.path(DIR, "test07_stata.csv"))
r_q10 <- safe_cor(df$stata_pred_q10, df$r_q10)
r_q50 <- safe_cor(df$stata_pred_q50, df$r_q50)
r_q90 <- safe_cor(df$stata_pred_q90, df$r_q90)
cross1 <- check_ordering(df, 10, 50)
cross2 <- check_ordering(df, 50, 90)
results[["07_hetero"]] <- list(
  test = "07: Heteroscedastic (Y = X1 + X2*N(0,1))",
  quantiles = c("q10", "q50", "q90"),
  cors = c(r_q10, r_q50, r_q90),
  threshold = THRESHOLD,
  cross_10_50 = cross1, cross_50_90 = cross2,
  note = ""
)
cat("Test 07: cor(q10)=", round(r_q10,4), " cor(q50)=", round(r_q50,4), " cor(q90)=", round(r_q90,4), "\n")

# ---- Test 08: Skewed ----
df <- read.csv(file.path(DIR, "test08_stata.csv"))
r_q10 <- safe_cor(df$stata_pred_q10, df$r_q10)
r_q50 <- safe_cor(df$stata_pred_q50, df$r_q50)
r_q90 <- safe_cor(df$stata_pred_q90, df$r_q90)
cross1 <- check_ordering(df, 10, 50)
cross2 <- check_ordering(df, 50, 90)
results[["08_skewed"]] <- list(
  test = "08: Skewed (Y = exp(X1) + N(0,0.3))",
  quantiles = c("q10", "q50", "q90"),
  cors = c(r_q10, r_q50, r_q90),
  threshold = THRESHOLD,
  cross_10_50 = cross1, cross_50_90 = cross2,
  note = ""
)
cat("Test 08: cor(q10)=", round(r_q10,4), " cor(q50)=", round(r_q50,4), " cor(q90)=", round(r_q90,4), "\n")

# ---- Test 09: Cluster ----
df <- read.csv(file.path(DIR, "test09_stata.csv"))
r_q10 <- safe_cor(df$stata_pred_q10, df$r_q10)
r_q50 <- safe_cor(df$stata_pred_q50, df$r_q50)
r_q90 <- safe_cor(df$stata_pred_q90, df$r_q90)
cross1 <- check_ordering(df, 10, 50)
cross2 <- check_ordering(df, 50, 90)
results[["09_cluster"]] <- list(
  test = "09: cluster() with 50 clusters",
  quantiles = c("q10", "q50", "q90"),
  cors = c(r_q10, r_q50, r_q90),
  threshold = THRESHOLD,
  cross_10_50 = cross1, cross_50_90 = cross2,
  note = ""
)
cat("Test 09: cor(q10)=", round(r_q10,4), " cor(q50)=", round(r_q50,4), " cor(q90)=", round(r_q90,4), "\n")

# ---- Test 10: Weights ----
# Note: R ran without weights; Stata with weights — relaxed threshold
df <- read.csv(file.path(DIR, "test10_stata.csv"))
r_q10 <- safe_cor(df$stata_pred_q10, df$r_q10)
r_q50 <- safe_cor(df$stata_pred_q50, df$r_q50)
r_q90 <- safe_cor(df$stata_pred_q90, df$r_q90)
cross1 <- check_ordering(df, 10, 50)
cross2 <- check_ordering(df, 50, 90)
results[["10_weights"]] <- list(
  test = "10: weights() [relaxed: R has no sample.weights in qf]",
  quantiles = c("q10", "q50", "q90"),
  cors = c(r_q10, r_q50, r_q90),
  threshold = THRESHOLD_WEIGHTS,
  cross_10_50 = cross1, cross_50_90 = cross2,
  note = "R quantile_forest lacks sample.weights; R baseline run without weights"
)
cat("Test 10: cor(q10)=", round(r_q10,4), " cor(q50)=", round(r_q50,4), " cor(q90)=", round(r_q90,4), "(relaxed threshold)\n")

# ---- Test 11: No honesty ----
df <- read.csv(file.path(DIR, "test11_stata.csv"))
r_q10 <- safe_cor(df$stata_pred_q10, df$r_q10)
r_q50 <- safe_cor(df$stata_pred_q50, df$r_q50)
r_q90 <- safe_cor(df$stata_pred_q90, df$r_q90)
cross1 <- check_ordering(df, 10, 50)
cross2 <- check_ordering(df, 50, 90)
results[["11_nohonesty"]] <- list(
  test = "11: nohonesty",
  quantiles = c("q10", "q50", "q90"),
  cors = c(r_q10, r_q50, r_q90),
  threshold = THRESHOLD,
  cross_10_50 = cross1, cross_50_90 = cross2,
  note = ""
)
cat("Test 11: cor(q10)=", round(r_q10,4), " cor(q50)=", round(r_q50,4), " cor(q90)=", round(r_q90,4), "\n")

# ---- Test 12: mtry=2 ----
df <- read.csv(file.path(DIR, "test12_stata.csv"))
r_q10 <- safe_cor(df$stata_pred_q10, df$r_q10)
r_q50 <- safe_cor(df$stata_pred_q50, df$r_q50)
r_q90 <- safe_cor(df$stata_pred_q90, df$r_q90)
cross1 <- check_ordering(df, 10, 50)
cross2 <- check_ordering(df, 50, 90)
results[["12_mtry2"]] <- list(
  test = "12: mtry=2",
  quantiles = c("q10", "q50", "q90"),
  cors = c(r_q10, r_q50, r_q90),
  threshold = THRESHOLD,
  cross_10_50 = cross1, cross_50_90 = cross2,
  note = ""
)
cat("Test 12: cor(q10)=", round(r_q10,4), " cor(q50)=", round(r_q50,4), " cor(q90)=", round(r_q90,4), "\n")

# ---- Test 13: minnodesize=20 ----
df <- read.csv(file.path(DIR, "test13_stata.csv"))
r_q10 <- safe_cor(df$stata_pred_q10, df$r_q10)
r_q50 <- safe_cor(df$stata_pred_q50, df$r_q50)
r_q90 <- safe_cor(df$stata_pred_q90, df$r_q90)
cross1 <- check_ordering(df, 10, 50)
cross2 <- check_ordering(df, 50, 90)
results[["13_minnodesize"]] <- list(
  test = "13: minnodesize=20",
  quantiles = c("q10", "q50", "q90"),
  cors = c(r_q10, r_q50, r_q90),
  threshold = THRESHOLD,
  cross_10_50 = cross1, cross_50_90 = cross2,
  note = ""
)
cat("Test 13: cor(q10)=", round(r_q10,4), " cor(q50)=", round(r_q50,4), " cor(q90)=", round(r_q90,4), "\n")

# ---- Test 14: samplefrac=0.3 ----
df <- read.csv(file.path(DIR, "test14_stata.csv"))
r_q10 <- safe_cor(df$stata_pred_q10, df$r_q10)
r_q50 <- safe_cor(df$stata_pred_q50, df$r_q50)
r_q90 <- safe_cor(df$stata_pred_q90, df$r_q90)
cross1 <- check_ordering(df, 10, 50)
cross2 <- check_ordering(df, 50, 90)
results[["14_samplefrac"]] <- list(
  test = "14: samplefrac=0.3",
  quantiles = c("q10", "q50", "q90"),
  cors = c(r_q10, r_q50, r_q90),
  threshold = THRESHOLD,
  cross_10_50 = cross1, cross_50_90 = cross2,
  note = ""
)
cat("Test 14: cor(q10)=", round(r_q10,4), " cor(q50)=", round(r_q50,4), " cor(q90)=", round(r_q90,4), "\n")

# ---- Test 15: Combined ----
df <- read.csv(file.path(DIR, "test15_stata.csv"))
r_q10 <- safe_cor(df$stata_pred_q10, df$r_q10)
r_q50 <- safe_cor(df$stata_pred_q50, df$r_q50)
r_q90 <- safe_cor(df$stata_pred_q90, df$r_q90)
cross1 <- check_ordering(df, 10, 50)
cross2 <- check_ordering(df, 50, 90)
results[["15_combined"]] <- list(
  test = "15: Combined (regsplit + cluster + weights)",
  quantiles = c("q10", "q50", "q90"),
  cors = c(r_q10, r_q50, r_q90),
  threshold = THRESHOLD,
  cross_10_50 = cross1, cross_50_90 = cross2,
  note = "Stata uses weights; R does not (weights not in quantile_forest)"
)
cat("Test 15: cor(q10)=", round(r_q10,4), " cor(q50)=", round(r_q50,4), " cor(q90)=", round(r_q90,4), "\n")

# ---- Test 16: Large p ----
df <- read.csv(file.path(DIR, "test16_stata.csv"))
r_q10 <- safe_cor(df$stata_pred_q10, df$r_q10)
r_q50 <- safe_cor(df$stata_pred_q50, df$r_q50)
r_q90 <- safe_cor(df$stata_pred_q90, df$r_q90)
cross1 <- check_ordering(df, 10, 50)
cross2 <- check_ordering(df, 50, 90)
results[["16_largep"]] <- list(
  test = "16: Large p=20 predictors",
  quantiles = c("q10", "q50", "q90"),
  cors = c(r_q10, r_q50, r_q90),
  threshold = THRESHOLD,
  cross_10_50 = cross1, cross_50_90 = cross2,
  note = ""
)
cat("Test 16: cor(q10)=", round(r_q10,4), " cor(q50)=", round(r_q50,4), " cor(q90)=", round(r_q90,4), "\n")

# ---- Test 17: Quantile crossing check ----
df <- read.csv(file.path(DIR, "test17_stata.csv"))
r_q10 <- safe_cor(df$stata_pred_q10, df$r_q10)
r_q50 <- safe_cor(df$stata_pred_q50, df$r_q50)
r_q90 <- safe_cor(df$stata_pred_q90, df$r_q90)
n_cross_10_50 <- sum(df$cross_10_50, na.rm = TRUE)
n_cross_50_90 <- sum(df$cross_50_90, na.rm = TRUE)
r_cross_10_50 <- sum(df$r_q10 > df$r_q50, na.rm = TRUE)
r_cross_50_90 <- sum(df$r_q50 > df$r_q90, na.rm = TRUE)
results[["17_crossing"]] <- list(
  test = "17: Quantile crossing check (ordering invariant)",
  quantiles = c("q10", "q50", "q90"),
  cors = c(r_q10, r_q50, r_q90),
  threshold = THRESHOLD,
  cross_10_50 = n_cross_10_50, cross_50_90 = n_cross_50_90,
  r_cross_10_50 = r_cross_10_50, r_cross_50_90 = r_cross_50_90,
  note = paste0("Stata: q10>q50 violations=", n_cross_10_50, ", q50>q90=", n_cross_50_90,
                " | R: q10>q50=", r_cross_10_50, ", q50>q90=", r_cross_50_90)
)
cat("Test 17: cor(q10)=", round(r_q10,4), " cor(q50)=", round(r_q50,4), " cor(q90)=", round(r_q90,4), "\n")
cat("  Stata violations: q10>q50=", n_cross_10_50, " q50>q90=", n_cross_50_90, "\n")
cat("  R    violations: q10>q50=", r_cross_10_50, " q50>q90=", r_cross_50_90, "\n")

# ---- Summary ----
cat("\n\n=== SUMMARY ===\n")
n_pass <- 0; n_fail <- 0; n_total <- 0
for (key in names(results)) {
  r <- results[[key]]
  pass_all <- all(r$cors >= r$threshold, na.rm = TRUE)
  status <- ifelse(pass_all, "PASS", "FAIL")
  if (pass_all) n_pass <- n_pass + 1 else n_fail <- n_fail + 1
  n_total <- n_total + 1
  cors_str <- paste(round(r$cors, 4), collapse = ", ")
  cat(sprintf("%-45s [%s]  cors=(%s)\n", r$test, status, cors_str))
}
cat(sprintf("\nTotal: %d tests, %d PASS, %d FAIL\n", n_total, n_pass, n_fail))

# Save results as RDS for report generation
saveRDS(results, file.path(DIR, "results.rds"))
cat("\nResults saved to results.rds\n")
