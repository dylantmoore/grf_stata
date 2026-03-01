## compute_correlations.R
## Read paired R/Stata CSV outputs and compute Pearson correlations per outcome

OUTDIR <- "/tmp/grf_stata/tests/fidelity_reports/13_multi_regression"

results <- list()

safe_cor <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3) return(NA_real_)
  cor(x[ok], y[ok])
}

## ----------------------------------------------------------------
## Helper: compute cors for a given test CSV
## ----------------------------------------------------------------
check_test <- function(test_id, csv_file, outcomes) {
  path <- file.path(OUTDIR, csv_file)
  if (!file.exists(path)) {
    cat(sprintf("[%s] MISSING: %s\n", test_id, csv_file))
    return(NULL)
  }
  df <- read.csv(path)
  cors <- sapply(outcomes, function(k) {
    rc <- paste0("r_pred_y", k)
    sc <- paste0("stata_pred_y", k)
    if (!all(c(rc, sc) %in% names(df))) {
      cat(sprintf("  [%s] Missing column for outcome %d\n", test_id, k))
      return(NA_real_)
    }
    safe_cor(df[[rc]], df[[sc]])
  })
  names(cors) <- paste0("Y", outcomes)
  cors
}

## ----------------------------------------------------------------
## Run all 14 tests
## ----------------------------------------------------------------

tests <- list(
  list(id="01", label="2 outcomes (default)", csv="test01_stata.csv", K=1:2),
  list(id="02", label="3 outcomes",            csv="test02_stata.csv", K=1:3),
  list(id="03", label="5 outcomes",            csv="test03_stata.csv", K=1:5),
  list(id="04", label="Correlated outcomes",   csv="test04_stata.csv", K=1:2),
  list(id="05", label="Independent outcomes",  csv="test05_stata.csv", K=1:2),
  list(id="06", label="cluster()",             csv="test06_stata.csv", K=1:2),
  list(id="07", label="weights()",             csv="test07_stata.csv", K=1:2),
  list(id="08", label="nohonesty",             csv="test08_stata.csv", K=1:2),
  list(id="09", label="mtry=2",                csv="test09_stata.csv", K=1:2),
  list(id="10", label="minnodesize=20",        csv="test10_stata.csv", K=1:2),
  list(id="11", label="samplefrac=0.3",        csv="test11_stata.csv", K=1:2),
  list(id="12", label="Combined: cluster+weights+nohonesty", csv="test12_stata.csv", K=1:2),
  list(id="13", label="Linear outcomes",       csv="test13_stata.csv", K=1:2),
  list(id="14", label="Nonlinear outcomes",    csv="test14_stata.csv", K=1:2)
)

all_results <- list()
for (t in tests) {
  cors <- check_test(t$id, t$csv, t$K)
  if (is.null(cors)) next
  pass <- all(cors > 0.90, na.rm = TRUE)
  status <- ifelse(pass, "PASS", "FAIL")
  cat(sprintf("[%s] %s -> %s\n", t$id, t$label, status))
  for (k in seq_along(cors)) {
    cat(sprintf("  Y%d cor = %.4f %s\n", t$K[k], cors[k],
                ifelse(cors[k] > 0.90, "ok", "BELOW 0.90")))
  }
  all_results[[t$id]] <- list(
    id = t$id,
    label = t$label,
    cors = cors,
    status = status
  )
}

## ----------------------------------------------------------------
## Summary table
## ----------------------------------------------------------------
cat("\n=== SUMMARY ===\n")
cat(sprintf("%-4s  %-45s  %-16s  %-6s\n", "ID", "Test", "Correlations", "Status"))
cat(strrep("-", 80), "\n")
for (r in all_results) {
  cor_str <- paste(sprintf("%.4f", r$cors), collapse=" | ")
  cat(sprintf("%-4s  %-45s  %-16s  %-6s\n",
      r$id, r$label, cor_str, r$status))
}

n_pass <- sum(sapply(all_results, function(r) r$status == "PASS"))
n_fail <- sum(sapply(all_results, function(r) r$status == "FAIL"))
cat(sprintf("\nTotal: %d PASS, %d FAIL out of %d tests\n", n_pass, n_fail, length(all_results)))

## Save machine-readable results
saveRDS(all_results, file.path(OUTDIR, "correlation_results.rds"))
cat("Saved to correlation_results.rds\n")
