## Compare R and Stata ATE results

r   <- read.csv("/tmp/grf_stata/tests/fidelity_reports/03_ate/r_results.csv",
               stringsAsFactors=FALSE)
st  <- read.csv("/tmp/grf_stata/tests/fidelity_reports/03_ate/stata_results.csv",
               stringsAsFactors=FALSE)

# Stata uses "." for missing — convert to NA
st$ate <- suppressWarnings(as.numeric(st$ate))
st$se  <- suppressWarnings(as.numeric(st$se))
r$ate  <- suppressWarnings(as.numeric(r$ate))
r$se   <- suppressWarnings(as.numeric(r$se))

## Merge
m <- merge(r, st, by=c("test_id","test_name"), suffixes=c("_r","_stata"))
m <- m[order(m$test_id),]

## Build comparison table
cat("\n=== ATE Fidelity Comparison: R vs Stata ===\n\n")
cat(sprintf("%-4s %-25s %10s %10s %10s %8s %8s %8s %7s  %-8s\n",
    "ID", "Name", "ATE_R", "ATE_Stata", "Diff", "|z|",
    "SE_R", "SE_Stata", "SE_reld", "Result"))
cat(paste(rep("-",100), collapse=""), "\n")

outcomes <- character(nrow(m))

for (i in seq_len(nrow(m))) {
  row  <- m[i,]
  tid  <- row$test_id
  tn   <- row$test_name
  ar   <- as.numeric(row$ate_r)
  ast  <- as.numeric(row$ate_stata)
  sr   <- as.numeric(row$se_r)
  sst  <- as.numeric(row$se_stata)
  err_r  <- row$r_errored
  err_st <- row$errored

  # Case 1: Both errored
  if (is.na(ar) && (is.na(ast) || err_st == 1)) {
    cat(sprintf("%-4d %-25s  Both errored/NA — behavioral match\n", tid, tn))
    outcomes[i] <- "PASS (both error)"
    next
  }
  # Case 2: Only R errored, Stata gave a value
  if (is.na(ar) && !is.na(ast)) {
    cat(sprintf("%-4d %-25s  R_errored=NA  Stata=%.6f — behavioral diff\n",
                tid, tn, ast))
    outcomes[i] <- "BEHAVIORAL_DIFF"
    next
  }
  # Case 3: Only Stata errored, R gave value
  if (!is.na(ar) && is.na(ast)) {
    cat(sprintf("%-4d %-25s  R=%.6f  Stata=NA (errored=%d) — behavioral diff\n",
                tid, tn, ar, err_st))
    outcomes[i] <- "BEHAVIORAL_DIFF"
    next
  }

  diff   <- ar - ast
  maxse  <- max(sr, sst, na.rm=TRUE)
  z      <- abs(diff) / maxse
  serdiff<- abs(sr - sst) / max(sr, 1e-12)
  zpass  <- z < 3
  sepass <- serdiff < 0.5
  result <- if(zpass && sepass) "PASS" else if(zpass) "PASS(z)|FAIL(se)" else "FAIL"
  outcomes[i] <- result
  cat(sprintf("%-4d %-25s %10.5f %10.5f %10.5f %8.3f %8.5f %8.5f %7.3f  %s\n",
              tid, tn, ar, ast, diff, z, sr, sst, serdiff, result))
}

cat("\n\nSummary:\n")
cat("PASS:            ", sum(outcomes == "PASS"), "\n")
cat("PASS(z)|FAIL(se):", sum(outcomes == "PASS(z)|FAIL(se)"), "\n")
cat("FAIL:            ", sum(outcomes == "FAIL"), "\n")
cat("BEHAVIORAL_DIFF: ", sum(outcomes == "BEHAVIORAL_DIFF"), "\n")
cat("PASS (both err): ", sum(outcomes == "PASS (both error)"), "\n")

# Save for report
write.csv(data.frame(test_id=m$test_id, test_name=m$test_name,
                     ate_r=m$ate_r, se_r=m$se_r,
                     ate_stata=m$ate_stata, se_stata=m$se_stata,
                     outcome=outcomes),
          "/tmp/grf_stata/tests/fidelity_reports/03_ate/comparison.csv",
          row.names=FALSE)
