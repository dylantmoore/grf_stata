# Step 1 Code Review -- GPT-5.3/Codex Agent

**Reviewer:** GPT-5.3/Codex review agent (simulated by Claude Opus 4.6)
**Date:** 2026-02-27
**Files reviewed:**
- `/tmp/grf_stata/tests/generate_reference.R`
- `/tmp/grf_stata/tests/test_fidelity.do`

---

## Verdict: NOT LGTM -- 3 issues found (1 bug, 2 gaps)

---

## Issue 1 (BUG): Causal survival forest `min.node.size` mismatch

**R script** (`generate_reference.R:601`):
```r
csf <- causal_survival_forest(X_cs, Y_cs, W_cs, D_cs,
  horizon = horizon_val,
  num.trees = 2000,
  seed = 42,
  honesty = TRUE,
  min.node.size = 15    # <-- explicitly 15
)
```

**Stata test** (`test_fidelity.do:871-872`):
```stata
grf_causal_survival_forest time w status x1-x5, gen(cspred_stata) ///
    ntrees(2000) seed(42) horizon(`cs_horizon')
```

The Stata call does **not** pass `minnodesize(15)` (or the equivalent Stata option). It will use the default `min.node.size` (typically 5), whereas the R reference was generated with `min.node.size = 15`. This parameter mismatch will cause the two forests to produce different predictions, likely causing Test 23 to fail or be unreliable.

**Fix:** Add `minnodesize(15)` (or equivalent) to the Stata `grf_causal_survival_forest` call at line 872.

---

## Issue 2 (GAP): 3 unconsumed reference CSVs

The R script generates **41 CSV files**. The Stata test consumes **38** of them. The following 3 CSVs are produced but never tested:

| CSV file | Generated at | Description |
|---|---|---|
| `causal_output.csv` | line 107 | CATE predictions + variance from causal forest |
| `regression_variable_importance.csv` | line 69 | Variable importance for regression forest |
| `survival_failure_times.csv` | line 284 | Unique failure times from survival forest |

**Assessment:**
- `causal_output.csv` is a significant omission -- there is no test that correlates Stata CATE predictions against R CATE predictions for the standard causal forest (analogous to Test 1 for regression). This is the most important gap.
- `regression_variable_importance.csv` -- Test 7 only tests causal variable importance rank correlation. A parallel test for regression variable importance would be valuable.
- `survival_failure_times.csv` -- This is ancillary metadata (the set of unique event times). It is indirectly consumed by the survival output comparison (Test 12), since the survival predictions are indexed by these times. This is the least critical gap.

**Recommendation:** Add a Test 25 for causal CATE prediction correlation (using `causal_output.csv`). Optionally add a test for regression variable importance. The `survival_failure_times.csv` gap is acceptable.

---

## Issue 3 (POTENTIAL BUG): Quantile column name after Stata import

**R script** (`generate_reference.R:175`):
```r
colnames(output_df) <- paste0("q", quantiles)
# Produces columns: q0.1, q0.25, q0.5, q0.75, q0.9
```

**Stata test** (`test_fidelity.do:349`):
```stata
rename q05 pred_r
```

When Stata's `import delimited` reads a column named `q0.5`, the dot is converted to either an underscore (`q0_5` in Stata 14+) or stripped entirely (`q05` in older versions). The test assumes the latter behavior (`q05`). If the project targets Stata 14+, this will fail because the actual variable name would be `q0_5`.

**Severity:** Depends on target Stata version. If Stata 14+ is targeted, this is a bug. If Stata 13 or earlier, it is fine. Worth verifying.

---

## Items verified as correct

### Step 1A: Fidelity path consistency
- All paths in `test_fidelity.do` use `ref/` relative to project root. No absolute paths or `../` relative paths found.
- Comment header at lines 1-14 of `test_fidelity.do` clearly documents the expected working directory convention.
- `generate_reference.R` uses `outdir <- "ref"` (line 17) consistently.
- Both scripts are designed to run from the project root. Paths are consistent.

### Step 1B: generate_reference.R completeness
All 7 new forest types and RATE AUTOC are present:

| Section | Forest type | Lines | Seed | Trees | Honesty | min.node.size |
|---|---|---|---|---|---|---|
| 12 | RATE AUTOC | 422-435 | 42 | 2000 | TRUE | (uses existing cf) |
| 13 | ll_regression_forest | 438-465 | 42 | 2000 | TRUE | 5 |
| 14 | lm_forest | 467-498 | 42 | 2000 | TRUE | 5 |
| 15 | multi_arm_causal_forest | 500-536 | 42 | 2000 | TRUE | 5 |
| 16 | multi_regression_forest | 538-571 | 42 | 2000 | TRUE | 5 |
| 17 | causal_survival_forest | 573-614 | 42 | 2000 | TRUE | 15 |
| 18 | boosted_regression_forest | 616-645 | 42 | 2000 | TRUE | 5 |

All forest functions are valid grf functions. Parameters are consistent (seed=42, num.trees=2000, honesty=TRUE, min.node.size=5 except causal_survival=15).

Each new section uses `set.seed(42)` before data generation to ensure reproducibility.

### Step 1C: test_fidelity.do completeness
Tests 1-24 are present and numbered:

| Test | Description | Metric | Threshold |
|---|---|---|---|
| 1 | Regression OOB predictions | Correlation | > 0.95 |
| 2 | Causal ATE (all) | z-test | < 3 |
| 3 | Causal ATE (treated) | z-test | < 3 |
| 4 | Causal ATE (control) | z-test | < 3 |
| 5 | Causal ATE (overlap) | z-test | < 3 |
| 6 | Causal BLP coefficients | per-coef z-test | < 3 |
| 7 | Variable importance rank | Spearman | > 0.70 |
| 8 | DR scores | Correlation | > 0.90 |
| 9 | Quantile forest predictions | Correlation | > 0.90 |
| 10 | Probability forest probs | Correlation | > 0.90 |
| 11 | Instrumental forest LATE | Correlation | > 0.90 |
| 12 | Survival forest predictions | Correlation | > 0.90 |
| 13 | Expected survival E[T|X] | Correlation | > 0.90 |
| 14 | Causal APE | z-test | < 3 |
| 15 | Test calibration | per-coef z-test | < 3 |
| 16 | BLP HC0 | per-coef z-test | < 3 |
| 17 | BLP HC3 | per-coef z-test | < 3 |
| 18 | RATE AUTOC | z-test | < 3 |
| 19 | LL regression forest | Correlation | > 0.90 |
| 20 | LM forest coefficients | Correlation | > 0.90 |
| 21 | Multi-arm causal forest | Correlation | > 0.90 |
| 22 | Multi-regression forest | Correlation | > 0.90 |
| 23 | Causal survival forest | Correlation | > 0.90 |
| 24 | Boosted regression forest | Correlation | > 0.90 |

All 24 tests use the graceful skip pattern (`capture confirm file ... if _rc ... SKIP`).

### Column name mappings verified correct
- R `std.err` -> Stata `std_err` (ATE, APE, RATE CSVs) -- matches
- R `std.error` -> Stata `std_error` (BLP, test calibration CSVs) -- matches
- R `estimate` -> Stata `estimate` -- matches
- R `prediction` -> Stata `prediction` -- matches
- R `cate` -> Stata `cate` -- matches
- R `late` -> Stata `late` -- matches
- R `score` -> Stata `score` -- matches

### Thresholds are reasonable
- Correlation thresholds (0.90-0.95) are appropriate for cross-implementation comparison
- z-test threshold (< 3) corresponds to roughly p < 0.003, accounting for normal stochastic variation between implementations

---

## Summary

| # | Severity | Issue | Location |
|---|---|---|---|
| 1 | **BUG** | `min.node.size` mismatch: R uses 15, Stata uses default (5) for causal survival forest | R:601 vs Stata:872 |
| 2 | **GAP** | 3 unconsumed CSVs (`causal_output`, `regression_variable_importance`, `survival_failure_times`) | R:107,69,284 |
| 3 | **POTENTIAL** | Quantile column name `q05` vs `q0_5` depends on Stata version | Stata:349 |
