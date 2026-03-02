# Step 1 Code Review -- Opus 4.6

**Files reviewed:**
- `/tmp/grf_stata/tests/generate_reference.R`
- `/tmp/grf_stata/tests/test_fidelity.do`

**Date:** 2026-02-27

---

## 1. Step 1A: Path Consistency and Header

**Status: PASS**

- `generate_reference.R` line 17: `outdir <- "ref"` -- paths are relative to CWD.
- `test_fidelity.do` lines 1-14: Comment header clearly documents expected working directory (project root), notes that both scripts use `ref/` relative to CWD, and states all tests gracefully skip if files are missing.
- Every CSV reference in test_fidelity.do uses `"ref/..."` paths (confirmed across all 24 tests). No stale absolute paths or `tests/ref/` paths found.
- `generate_reference.R` uses `file.path(outdir, ...)` consistently via the `save_csv` helper.

No issues.

---

## 2. Step 1B: generate_reference.R -- New Forest Types and RATE AUTOC

**Status: PASS with 1 issue (see Issue #1 below)**

### CSV Inventory (all CSVs produced by generate_reference.R)

| # | Section | CSV filename | Column names |
|---|---------|-------------|--------------|
| 1 | Regression | regression_input.csv | x1..x5, y |
| 2 | Regression | regression_output.csv | prediction, variance |
| 3 | Regression | regression_variable_importance.csv | variable, importance |
| 4 | Causal | causal_input.csv | x1..x5, y, w |
| 5 | Causal | causal_output.csv | cate, variance |
| 6 | Causal | causal_ate.csv | estimate, std.err |
| 7 | Causal | causal_variable_importance.csv | variable, importance |
| 8 | Causal | causal_blp.csv | variable, estimate, std.error, t.value, p.value |
| 9 | Causal | causal_test_calibration.csv | term, estimate, std.error, t.value, p.value |
| 10 | Quantile | quantile_input.csv | x1..x5, y |
| 11 | Quantile | quantile_output.csv | q0.1, q0.25, q0.5, q0.75, q0.9 |
| 12 | Instrumental | instrumental_input.csv | x1..x5, y, w, z |
| 13 | Instrumental | instrumental_output.csv | late, variance |
| 14 | Probability | probability_input.csv | x1..x5, y |
| 15 | Probability | probability_output.csv | class0, class1, class2 |
| 16 | Survival | survival_input.csv | x1..x5, time, status |
| 17 | Survival | survival_failure_times.csv | failure_time |
| 18 | Survival | survival_output.csv | t1..t20 |
| 19 | Causal ATE | causal_ate_treated.csv | estimate, std.err |
| 20 | Causal ATE | causal_ate_control.csv | estimate, std.err |
| 21 | Causal ATE | causal_ate_overlap.csv | estimate, std.err |
| 22 | APE | causal_ape_input.csv | x1..x5, y, w |
| 23 | APE | causal_ape.csv | estimate, std.err |
| 24 | BLP HC0 | causal_blp_hc0.csv | variable, estimate, std.error, t.value, p.value |
| 25 | BLP HC3 | causal_blp_hc3.csv | variable, estimate, std.error, t.value, p.value |
| 26 | DR Scores | causal_scores.csv | score |
| 27 | Expected Surv | survival_expected.csv | expected_time |
| 28 | RATE AUTOC | causal_rate_autoc.csv | estimate, std.err |
| 29 | LL Regression | ll_regression_input.csv | x1..x5, y |
| 30 | LL Regression | ll_regression_output.csv | prediction |
| 31 | LM Forest | lm_forest_input.csv | x1..x5, y |
| 32 | LM Forest | lm_forest_output.csv | coef1..coefN |
| 33 | Multi-arm | multi_arm_input.csv | x1..x5, y, w |
| 34 | Multi-arm | multi_arm_output.csv | contrast1..contrastN |
| 35 | Multi-regression | multi_regression_input.csv | x1..x5, y1, y2 |
| 36 | Multi-regression | multi_regression_output.csv | pred_y1, pred_y2 |
| 37 | Causal Survival | causal_survival_input.csv | x1..x5, time, status, w |
| 38 | Causal Survival | causal_survival_output.csv | cate |
| 39 | Causal Survival | causal_survival_horizon.csv | horizon |
| 40 | Boosted Reg | boosted_regression_input.csv | x1..x5, y |
| 41 | Boosted Reg | boosted_regression_output.csv | prediction |

**Total: 41 CSVs generated.**

### R function calls and parameters -- correctness check

All forest fits use consistent parameters:
- `num.trees = 2000` -- correct
- `seed = 42` -- correct
- `honesty = TRUE` -- correct
- `min.node.size = 5` -- correct for all except causal_survival which correctly uses `min.node.size = 15` (line 601)

Specific checks:
- **ll_regression_forest** (line 453): Correct function name and signature. `ll_regression_forest(X, Y, ...)` is the grf API.
- **lm_forest** (line 483): `lm_forest(X_lm, Y_lm, X_lm, ...)` -- the third argument is the gradient covariates `W`. Passing `X_lm` as `W` means the model estimates local linear coefficients for all X variables. This is a valid usage.
- **multi_arm_causal_forest** (line 522): `multi_arm_causal_forest(X, Y, W, ...)` with `W` as a factor. Correct API.
- **multi_regression_forest** (line 558): `multi_regression_forest(X, Y, ...)` with Y as a matrix. Correct API.
- **causal_survival_forest** (line 596): `causal_survival_forest(X, Y, W, D, horizon=..., ...)`. Correct API. `horizon` is computed as `median(Y_cs[D_cs == 1])` which is a reasonable choice.
- **boosted_regression_forest** (line 632): Correct function name and signature.
- **rank_average_treatment_effect** (line 429): Called as `rank_average_treatment_effect(cf, pred_cf_cate)`. Correct -- the second argument is the prioritization scores (using CATE predictions gives AUTOC).
- **best_linear_projection** with vcov.type (lines 356, 368): Correct parameter name `vcov.type`.
- **average_partial_effect** (line 345): Called on a continuous-treatment causal forest. Correct usage.

### Seed management

The R script calls `set.seed(42)` at line 16 (global). Sections 1-12 do NOT re-seed, so the RNG state flows sequentially. However, sections 13-18 (the new forest types) each call `set.seed(42)` again (lines 442, 472, 505, 543, 578, 620). This means:
- Sections 13-18 each start with a fresh seed(42), making them individually reproducible.
- Sections 1-12 depend on cumulative RNG state from seed(42) at line 16.

This is intentional and correct for reproducibility of each new forest type independently.

---

## 3. Step 1C: test_fidelity.do -- Fidelity Comparisons

**Status: PASS with 2 issues (see Issues #1, #2 below)**

### Test inventory and CSV consumption mapping

| Test # | Description | Input CSV | Reference CSV | Metric | Threshold |
|--------|------------|-----------|--------------|--------|-----------|
| 1 | Regression OOB | regression_input | regression_output | correlation | > 0.95 |
| 2 | Causal ATE (all) | causal_input | causal_ate | z-test | < 3 |
| 3 | Causal ATE (treated) | causal_input | causal_ate_treated | z-test | < 3 |
| 4 | Causal ATE (control) | causal_input | causal_ate_control | z-test | < 3 |
| 5 | Causal ATE (overlap) | causal_input | causal_ate_overlap | z-test | < 3 |
| 6 | Causal BLP | causal_input | causal_blp | z-test | < 3 |
| 7 | Variable importance | causal_input | causal_variable_importance | Spearman | > 0.70 |
| 8 | DR scores | causal_input | causal_scores | correlation | > 0.90 |
| 9 | Quantile | quantile_input | quantile_output | correlation | > 0.90 |
| 10 | Probability | probability_input | probability_output | correlation | > 0.90 |
| 11 | Instrumental | instrumental_input | instrumental_output | correlation | > 0.90 |
| 12 | Survival | survival_input | survival_output | correlation | > 0.90 |
| 13 | Expected survival | survival_input | survival_expected | correlation | > 0.90 |
| 14 | Causal APE | causal_ape_input | causal_ape | z-test | < 3 |
| 15 | Test calibration | causal_input | causal_test_calibration | z-test | < 3 |
| 16 | BLP HC0 | causal_input | causal_blp_hc0 | z-test | < 3 |
| 17 | BLP HC3 | causal_input | causal_blp_hc3 | z-test | < 3 |
| 18 | RATE AUTOC | causal_input | causal_rate_autoc | z-test | < 3 |
| 19 | LL regression | ll_regression_input | ll_regression_output | correlation | > 0.90 |
| 20 | LM forest | lm_forest_input | lm_forest_output | correlation | > 0.90 |
| 21 | Multi-arm | multi_arm_input | multi_arm_output | correlation | > 0.90 |
| 22 | Multi-regression | multi_regression_input | multi_regression_output | correlation | > 0.90 |
| 23 | Causal survival | causal_survival_input | causal_survival_output + causal_survival_horizon | correlation | > 0.90 |
| 24 | Boosted regression | boosted_regression_input | boosted_regression_output | correlation | > 0.90 |

**Total: 24 tests (tests 9-24 are the 16 new ones). Meets the 16+ requirement.**

### Graceful skip check

All 24 tests use the `capture confirm file "ref/..."` pattern followed by `if _rc { display as text "SKIP: ..." }`. Confirmed correct for every test block.

### Column name cross-references (R CSV columns vs. Stata `rename` / access)

All verified correct:
- Test 9 (quantile): R produces `q0.1, q0.25, q0.5, q0.75, q0.9`. Stata imports and renames `q05` to `pred_r`, correlates with `qpred_q50`. The column name `q05` is what `import delimited` will produce from `q0.5` (dots become empty in Stata variable names). This is correct behavior for Stata's `import delimited`.
- Test 10 (probability): R produces `class0, class1, class2`. Stata renames `class0` to `prob_r`.
- Test 11 (instrumental): R produces `late, variance`. Stata renames `late` to `pred_r`.
- Test 12 (survival): R produces `t1..t20`. Stata renames `t1` to `pred_r`.
- Test 13 (expected survival): R produces `expected_time`. Stata renames `expected_time` to `et_r`.
- Test 14 (APE): R produces `estimate, std.err`. Stata reads `estimate[1]` and `std_err[1]`. Note: `std.err` becomes `std_err` in Stata import. Correct.
- Test 15 (test calibration): R produces `term, estimate, std.error, t.value, p.value`. Stata reads `estimate[j]` and `std_error[j]`. Note: `std.error` becomes `std_error` in Stata import. Correct.
- Test 18 (RATE AUTOC): R produces `estimate, std.err`. Stata reads `estimate[1]` and `std_err[1]`. Correct.
- Test 19 (LL regression): R produces `prediction`. Stata renames `prediction` to `pred_r`. Correct.
- Test 20 (LM forest): R produces `coef1, coef2, ...`. Stata renames `coef1` to `pred_r`. Correct.
- Test 21 (multi-arm): R produces `contrast1, contrast2, ...`. Stata renames `contrast1` to `pred_r`. Correct.
- Test 22 (multi-regression): R produces `pred_y1, pred_y2`. Stata renames `pred_y1` to `pred_r`. Correct.
- Test 23 (causal survival): R produces `cate`. Stata renames `cate` to `pred_r`. R also produces `causal_survival_horizon.csv` with `horizon` column, which Stata reads correctly at line 868.
- Test 24 (boosted regression): R produces `prediction`. Stata renames `prediction` to `pred_r`. Correct.

### Stata forest command parameters

All Stata calls use `ntrees(2000) seed(42)`. Honesty and min.node.size are presumably defaults in the Stata wrapper (matching grf defaults honesty=TRUE, min.node.size=5). The causal_survival test at line 872 does NOT explicitly set `min.node.size(15)` but does pass `horizon()`, so this depends on the Stata wrapper's default. This is noted as Issue #2 below.

---

## 4. Completeness Check: Unconsumed CSVs

Cross-referencing all 41 CSVs from generate_reference.R against consumers in test_fidelity.do:

| CSV | Consumed? | Consumer test |
|-----|-----------|---------------|
| regression_input.csv | YES | Test 1 |
| regression_output.csv | YES | Test 1 |
| regression_variable_importance.csv | **NO** | -- |
| causal_input.csv | YES | Tests 2-8, 15-18 |
| causal_output.csv | **NO** | -- |
| causal_ate.csv | YES | Test 2 |
| causal_variable_importance.csv | YES | Test 7 |
| causal_blp.csv | YES | Test 6 |
| causal_test_calibration.csv | YES | Test 15 |
| quantile_input.csv | YES | Test 9 |
| quantile_output.csv | YES | Test 9 |
| instrumental_input.csv | YES | Test 11 |
| instrumental_output.csv | YES | Test 11 |
| probability_input.csv | YES | Test 10 |
| probability_output.csv | YES | Test 10 |
| survival_input.csv | YES | Tests 12-13 |
| survival_failure_times.csv | **NO** | -- |
| survival_output.csv | YES | Test 12 |
| causal_ate_treated.csv | YES | Test 3 |
| causal_ate_control.csv | YES | Test 4 |
| causal_ate_overlap.csv | YES | Test 5 |
| causal_ape_input.csv | YES | Test 14 |
| causal_ape.csv | YES | Test 14 |
| causal_blp_hc0.csv | YES | Test 16 |
| causal_blp_hc3.csv | YES | Test 17 |
| causal_scores.csv | YES | Test 8 |
| survival_expected.csv | YES | Test 13 |
| causal_rate_autoc.csv | YES | Test 18 |
| ll_regression_input.csv | YES | Test 19 |
| ll_regression_output.csv | YES | Test 19 |
| lm_forest_input.csv | YES | Test 20 |
| lm_forest_output.csv | YES | Test 20 |
| multi_arm_input.csv | YES | Test 21 |
| multi_arm_output.csv | YES | Test 21 |
| multi_regression_input.csv | YES | Test 22 |
| multi_regression_output.csv | YES | Test 22 |
| causal_survival_input.csv | YES | Test 23 |
| causal_survival_output.csv | YES | Test 23 |
| causal_survival_horizon.csv | YES | Test 23 |
| boosted_regression_input.csv | YES | Test 24 |
| boosted_regression_output.csv | YES | Test 24 |

---

## Issues Found

### Issue #1 (Minor): 3 unconsumed CSVs

**Severity: Low** -- These are generated but never tested:

1. **regression_variable_importance.csv** (generate_reference.R:69) -- Test 7 tests causal variable importance but no test consumes regression variable importance. This was a pre-existing CSV from the original code, not introduced in Step 1, so arguably out of scope. However, for zero-unconsumed-CSVs the plan requires either a test or removing the generation.

2. **causal_output.csv** (generate_reference.R:107) -- Contains per-observation CATE predictions and variance. No test compares Stata CATE predictions directly against this file. Test 1 does a prediction correlation test for regression but there is no analogous CATE prediction correlation test for causal. Again, pre-existing.

3. **survival_failure_times.csv** (generate_reference.R:284) -- Contains the unique failure times from the survival forest. No test consumes this. Pre-existing.

**Recommendation:** Either add correlation tests for these 3 CSVs (matching the pattern used in tests 1, 12, etc.) or stop generating them. Strictly speaking, the Step 1C requirement says "zero unconsumed CSVs," so this is a gap.

### Issue #2 (Medium): causal_survival min.node.size mismatch risk

**File:** test_fidelity.do:871-872
```stata
grf_causal_survival_forest time w status x1-x5, gen(cspred_stata) ///
    ntrees(2000) seed(42) horizon(`cs_horizon')
```

The R script uses `min.node.size = 15` for causal_survival_forest (generate_reference.R:601), which differs from the default of 5 used everywhere else. The Stata test does NOT explicitly pass `min.node.size(15)` or equivalent. If the Stata wrapper defaults to `min.node.size = 5`, the comparison will be against a differently-configured forest, potentially causing the test to fail or pass spuriously.

**Recommendation:** Add `minnodesize(15)` (or whatever the Stata option name is) to the Stata causal_survival call to match the R reference exactly.

### Issue #3 (Informational): honesty.fraction and ci.group.size only on regression forest

**File:** generate_reference.R:50-53

The regression forest explicitly sets `honesty.fraction = 0.5, alpha = 0.05, ci.group.size = 2` but no other forest does. These are default values in grf, so omitting them elsewhere is harmless. However, the Stata regression test also does not explicitly set these. This is consistent (both use defaults), so no action needed. Noting for completeness only.

---

## Summary

| Criterion | Verdict |
|-----------|---------|
| Step 1A: Path consistency | PASS |
| Step 1B: All 7 new reference data types generated | PASS |
| Step 1C: 24 tests (16 new) | PASS |
| Step 1C: Graceful skip on missing files | PASS |
| Correctness: Column names match | PASS |
| Correctness: R function calls and params | PASS |
| Correctness: Thresholds appropriate | PASS |
| Consistency: seed=42, ntrees=2000, honesty=TRUE | PASS |
| Consistency: min.node.size=5 (15 for causal_survival) | **ISSUE #2** -- Stata test missing explicit min.node.size=15 |
| Completeness: Zero unconsumed CSVs | **ISSUE #1** -- 3 pre-existing CSVs lack consumers |

## Verdict

**NOT LGTM** -- two concrete issues should be addressed:

1. **Issue #1**: 3 CSVs generated with no consumer test: `regression_variable_importance.csv`, `causal_output.csv`, `survival_failure_times.csv`. Add tests or remove generation.
2. **Issue #2**: `test_fidelity.do` line 871-872: causal_survival_forest call is missing explicit `min.node.size(15)` to match the R reference's `min.node.size = 15`.
