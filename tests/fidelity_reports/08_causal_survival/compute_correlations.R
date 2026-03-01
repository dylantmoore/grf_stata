#!/usr/bin/env Rscript
# Compute Pearson correlations between R and Stata predictions for all 16 tests

outdir <- "/tmp/grf_stata/tests/fidelity_reports/08_causal_survival"

# Test metadata
tests <- list(
  list(id=1,  desc="Default RMST, horizon=median(Y)",             target="RMST",    threshold=0.85),
  list(id=2,  desc="Explicit horizon=Q75(Y)=1.1489",              target="RMST",    threshold=0.85),
  list(id=3,  desc="Survival probability target",                  target="SP",      threshold=0.75),
  list(id=4,  desc="No stabilize splits",                          target="RMST",    threshold=0.85),
  list(id=5,  desc="User-supplied W.hat=0.5",                      target="RMST",    threshold=0.85),
  list(id=6,  desc="With cluster(50 groups)",                      target="RMST",    threshold=0.85),
  list(id=7,  desc="With observation weights",                     target="RMST",    threshold=0.85),
  list(id=8,  desc="No honesty",                                   target="RMST",    threshold=0.85),
  list(id=9,  desc="mtry=2",                                       target="RMST",    threshold=0.85),
  list(id=10, desc="min.node.size=20",                             target="RMST",    threshold=0.85),
  list(id=11, desc="Heavy censoring (~40% censored at rate=0.8)",  target="RMST",    threshold=0.85),
  list(id=12, desc="Balanced treatment 50/50 (same as Test 1)",    target="RMST",    threshold=0.85),
  list(id=13, desc="Unbalanced treatment 70/30",                   target="RMST",    threshold=0.85),
  list(id=14, desc="Fewer trees (num.trees=100)",                  target="RMST",    threshold=0.85),
  list(id=15, desc="Small horizon=Q25(failure times)=0.2173",      target="RMST",    threshold=0.85),
  list(id=16, desc="Large horizon=Q90(failure times)=2.0006",      target="RMST",    threshold=0.85)
)

results <- list()

cat("Computing R vs Stata correlations for causal_survival_forest...\n\n")
cat(sprintf("%-3s  %-50s  %6s  %6s  %6s  %6s  %6s  %-6s\n",
            "ID", "Description", "Cor", "R.mean", "S.mean", "R.sd", "S.sd", "Status"))
cat(strrep("-", 110), "\n")

for (t in tests) {
  tid <- sprintf("%02d", t$id)
  r_file <- file.path(outdir, sprintf("r_pred_%s.csv", tid))
  s_file <- file.path(outdir, sprintf("s_pred_%s.csv", tid))

  if (!file.exists(r_file) || !file.exists(s_file)) {
    cat(sprintf("%2d   %-50s  MISSING FILES\n", t$id, t$desc))
    results[[t$id]] <- list(id=t$id, desc=t$desc, status="MISSING", cor=NA)
    next
  }

  r_pred <- read.csv(r_file)$tau_hat
  s_pred <- read.csv(s_file)[[1]]

  if (length(r_pred) != length(s_pred)) {
    cat(sprintf("%2d   %-50s  LENGTH MISMATCH: R=%d Stata=%d\n",
                t$id, t$desc, length(r_pred), length(s_pred)))
    results[[t$id]] <- list(id=t$id, desc=t$desc, status="MISMATCH", cor=NA)
    next
  }

  # Handle any NA
  ok <- !is.na(r_pred) & !is.na(s_pred)
  n_ok <- sum(ok)

  if (n_ok < 10) {
    cat(sprintf("%2d   %-50s  TOO FEW NON-MISSING: %d\n", t$id, t$desc, n_ok))
    results[[t$id]] <- list(id=t$id, desc=t$desc, status="FEW_OBS", cor=NA)
    next
  }

  cor_val <- cor(r_pred[ok], s_pred[ok], method = "pearson")
  threshold <- t$threshold
  passed <- !is.na(cor_val) && cor_val > threshold

  status <- if (passed) "PASS" else "FAIL"
  cat(sprintf("%2d   %-50s  %6.4f  %6.4f  %6.4f  %6.4f  %6.4f  %-6s\n",
              t$id, substr(t$desc, 1, 50),
              cor_val,
              mean(r_pred[ok]), mean(s_pred[ok]),
              sd(r_pred[ok]), sd(s_pred[ok]),
              status))

  results[[t$id]] <- list(
    id = t$id, desc = t$desc, target = t$target,
    cor = cor_val, threshold = threshold,
    r_mean = mean(r_pred[ok]), s_mean = mean(s_pred[ok]),
    r_sd = sd(r_pred[ok]), s_sd = sd(s_pred[ok]),
    n = n_ok, n_ok = n_ok,
    r_min = min(r_pred[ok]), r_max = max(r_pred[ok]),
    s_min = min(s_pred[ok]), s_max = max(s_pred[ok]),
    status = status, passed = passed
  )
}

cat(strrep("-", 110), "\n")

n_pass <- sum(sapply(results, function(r) !is.null(r$passed) && r$passed))
n_fail <- sum(sapply(results, function(r) !is.null(r$passed) && !r$passed))
n_total <- length(results)

cat(sprintf("\nSummary: %d/%d PASS, %d FAIL\n", n_pass, n_total, n_fail))

# Save results
saveRDS(results, file.path(outdir, "correlation_results.rds"))
cat("Saved correlation_results.rds\n")
