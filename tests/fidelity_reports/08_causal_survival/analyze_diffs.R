#!/usr/bin/env Rscript
# Deep analysis of R vs Stata prediction differences

outdir <- "/tmp/grf_stata/tests/fidelity_reports/08_causal_survival"

cat("=== Deep analysis of R vs Stata causal_survival_forest predictions ===\n\n")

# Load results
results <- readRDS(file.path(outdir, "correlation_results.rds"))

# Focus on tests with lower correlation
low_cor_tests <- which(sapply(results, function(r) !is.null(r$cor) && r$cor < 0.85))

for (tid in 1:16) {
  r_file <- file.path(outdir, sprintf("r_pred_%02d.csv", tid))
  s_file <- file.path(outdir, sprintf("s_pred_%02d.csv", tid))

  if (!file.exists(r_file) || !file.exists(s_file)) next

  r_pred <- read.csv(r_file)$tau_hat
  s_pred <- read.csv(s_file)[[1]]

  cat(sprintf("Test %2d: ", tid))
  cat(sprintf("cor=%.4f  R[%+.4f±%.4f, %+.4f to %+.4f]  Stata[%+.4f±%.4f, %+.4f to %+.4f]\n",
              results[[tid]]$cor,
              mean(r_pred), sd(r_pred), min(r_pred), max(r_pred),
              mean(s_pred), sd(s_pred), min(s_pred), max(s_pred)))

  # Check Spearman correlation (rank-based - more robust to scale differences)
  spear <- cor(r_pred, s_pred, method="spearman")
  cat(sprintf("         Spearman=%.4f  ", spear))

  # Check if scale difference: if we rescale Stata predictions
  # lm(r ~ s)
  fit <- lm(r_pred ~ s_pred)
  cat(sprintf("  lm(R~Stata): b0=%.4f, b1=%.4f, R2=%.4f\n",
              coef(fit)[1], coef(fit)[2], summary(fit)$r.squared))
}

# Overall rank correlation
cat("\n\n=== Testing if sign agreement is consistent ===\n")
for (tid in 1:16) {
  r_file <- file.path(outdir, sprintf("r_pred_%02d.csv", tid))
  s_file <- file.path(outdir, sprintf("s_pred_%02d.csv", tid))
  if (!file.exists(r_file) || !file.exists(s_file)) next

  r_pred <- read.csv(r_file)$tau_hat
  s_pred <- read.csv(s_file)[[1]]

  sign_agree <- mean(sign(r_pred) == sign(s_pred))
  cat(sprintf("Test %2d: sign agreement = %.1f%%\n", tid, sign_agree * 100))
}

# Compare scale: ratio of SDs
cat("\n=== SD ratios (R vs Stata) ===\n")
for (tid in 1:16) {
  r_file <- file.path(outdir, sprintf("r_pred_%02d.csv", tid))
  s_file <- file.path(outdir, sprintf("s_pred_%02d.csv", tid))
  if (!file.exists(r_file) || !file.exists(s_file)) next

  r_pred <- read.csv(r_file)$tau_hat
  s_pred <- read.csv(s_file)[[1]]

  sd_ratio <- sd(r_pred) / sd(s_pred)
  cat(sprintf("Test %2d: SD(R)/SD(Stata) = %.3f  (R spread is %.1fx Stata spread)\n",
              tid, sd_ratio, sd_ratio))
}
