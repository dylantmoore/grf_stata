# Consolidated Review: grf_stata Test Suite

**Date:** 2026-02-27
**Reviewers consolidated:** Claude Opus 4.6 ("Opus"), OpenAI Codex ("Codex"), Google Gemini ("Gemini")
**Reconciler:** Claude Opus 4.6 (with code verification)
**Scope:** 24 .ado commands, 39 .do test files, 2 .R support scripts at `/tmp/grf_stata/`

---

## 1. Executive Summary

The grf_stata test suite provides strong functional coverage: every one of the 24 commands has dedicated tests, all 12 forest types are exercised for fit/predict/tune, and per-option test files exist for all forest types plus post-estimation commands. The principal weakness is fidelity testing against R -- only 2 of 12 forest types (regression and causal) have R-vs-Stata numerical comparisons in `test_fidelity.do`, despite the R reference generator producing data for 4 additional forest types and several additional post-estimation quantities. The most significant functional gaps versus R's grf package are the absence of `clusters`, `sample.weights`, and MIA missing-data passthrough, all of which are important for applied research.

---

## 2. Command Coverage Matrix (Unified, Reconciled)

All three reviewers agreed that every command has at least one test. Disagreements arose around specific untested options. The table below reflects verified findings after checking the actual test code.

### 2.1 Forest Estimation Commands (12)

| Command | Basic Test | Options Test | Verified Untested Options |
|---------|-----------|-------------|--------------------------|
| `grf_regression_forest` | `test_regression.do` | `test_options_regression.do` | None |
| `grf_causal_forest` | `test_causal.do` | `test_options_causal.do` | None |
| `grf_quantile_forest` | `test_quantile.do` | `test_options_quantile.do` | None |
| `grf_probability_forest` | `test_probability.do` | `test_options_probability.do` | None |
| `grf_instrumental_forest` | `test_instrumental.do` | `test_options_instrumental.do` | None |
| `grf_survival_forest` | `test_survival.do` | `test_options_survival.do` | None |
| `grf_ll_regression_forest` | `test_ll_regression.do` | `test_options_ll_regression.do` | None |
| `grf_lm_forest` | `test_lm_forest.do` | `test_options_lm_forest.do` | None |
| `grf_multi_arm_causal_forest` | `test_multi_arm.do` | `test_options_multi_arm.do` | None |
| `grf_multi_regression_forest` | `test_multi_regression.do` | `test_options_multi_regression.do` | None |
| `grf_boosted_regression_forest` | `test_boosted_regression.do` | `test_options_boosted.do` | None |
| `grf_causal_survival_forest` | `test_causal_survival.do` | `test_options_causal_survival.do` | None [A] |

**[A] Reconciliation note:** Codex flagged `numer()` and `denom()` as untested for `grf_causal_survival_forest`. This is **incorrect**. `test_causal_survival.do` (the basic test, lines 59-63) explicitly generates `my_numer` and `my_denom` variables and passes them via `numer(my_numer) denom(my_denom)`. All three reviewers agreed that the per-option test files cover every documented option for each forest type.

### 2.2 Post-Estimation Commands (10)

| Command | Primary Test Files | Verified Untested Options |
|---------|-------------------|--------------------------|
| `grf_predict` | `test_predict.do`, `test_options_post_estimation.do` | None |
| `grf_tune` | `test_tune.do`, `test_tune_extended.do`, `test_options_post_estimation.do` | `nclasses()`, `numfailures()`, `predtype()`, `reducedformweight()`, `horizon()`, `target()` [B] |
| `grf_ate` | `test_post_estimation.do`, `test_options_post_estimation.do` | None |
| `grf_test_calibration` | `test_post_estimation.do`, `test_options_post_estimation.do` | None |
| `grf_variable_importance` | `test_post_estimation.do`, `test_options_post_estimation.do` | None |
| `grf_best_linear_projection` | `test_post_estimation.do`, `test_options_post_estimation.do` | None |
| `grf_get_scores` | `test_get_scores.do` | None |
| `grf_rate` | `test_rate.do`, `test_options_post_estimation.do` | `in` qualifier [C] |
| `grf_average_partial_effect` | `test_average_partial_effect.do` | `numtreesvariance()`, `in` qualifier |
| `grf_expected_survival` | `test_expected_survival.do` | None |

**[B] Reconciliation note on grf_tune:** Opus originally flagged only `numthreads()` as untested for tune. Codex flagged `nclasses()`, `numfailures()`, `predtype()`, `reducedformweight()`, `horizon()`, `target()`. Gemini flagged nothing. After verification: `numthreads(2)` IS tested in `test_options_post_estimation.do` Test 31 (Opus was wrong). However, the 6 forest-type-specific options (`nclasses`, `numfailures`, `predtype`, `reducedformweight`, `horizon`, `target`) are never passed with non-default values in any tune test (Codex was correct). The tune tests for probability, survival, instrumental, and causal_survival forests use only default values for these parameters.

**[C] Reconciliation note on grf_rate:** Opus flagged `quantiles()` and `compliancescore()` as untested in `test_rate.do`. While these options do NOT appear in `test_rate.do` (which has only 3 tests), they ARE thoroughly tested in `test_options_post_estimation.do`: `quantiles()` in Test 19 and `compliancescore()` in Tests 69-78 (10 tests). So these options are covered; the `in` qualifier remains the only untested option.

### 2.3 Data Generation Commands (2)

| Command | Test File | Verified Coverage |
|---------|----------|-------------------|
| `grf_generate_causal_data` | `test_generate_data.do` | All 11 DGPs tested; error cases tested (invalid DGP, minimum p) |
| `grf_generate_causal_survival_data` | `test_generate_data.do` | All 6 DGPs tested; no error-case tests (invalid DGP, minimum p) |

---

## 3. Forest Type Coverage Matrix (Unified)

All three reviewers agreed on fit, predict, and tune coverage. The key addition here is fidelity testing status, which all reviewers identified but with varying detail.

| Forest Type | Fit | Predict | Tune | Options | Fidelity vs R |
|------------|-----|---------|------|---------|---------------|
| regression | Yes | Yes | Yes | Yes | TESTED (correlation > 0.95) |
| causal | Yes | Yes | Yes | Yes | TESTED (ATE z-test, BLP, VI, DR scores) |
| quantile | Yes | Yes | Yes | Yes | NOT TESTED (ref data exists) |
| probability | Yes | Yes | Yes | Yes | NOT TESTED (ref data exists) |
| instrumental | Yes | Yes | Yes | Yes | NOT TESTED (ref data exists) |
| survival | Yes | Yes | Yes | Yes | NOT TESTED (ref data exists) |
| ll_regression | Yes | Yes | Yes | Yes | NOT TESTED (no ref data) |
| lm_forest | Yes | Yes | Yes | Yes | NOT TESTED (no ref data) |
| multi_arm_causal | Yes | Yes | Yes | Yes | NOT TESTED (no ref data) |
| multi_regression | Yes | Yes | Yes | Yes | NOT TESTED (no ref data) |
| boosted_regression | Yes | Negative-path only [D] | Yes | Yes | NOT TESTED (no ref data) |
| causal_survival | Yes | Yes | Yes | Yes | NOT TESTED (no ref data) |

**[D] Reconciliation note:** Codex correctly identified that the boosted regression predict test (Test 13 in `test_predict.do`) is negative-path only -- it asserts that `grf_predict` errors out after fitting a boosted regression forest, rather than verifying successful prediction. Gemini marked predict as "Yes" without this nuance. Opus listed it as a standard predict test. After verification: `test_predict.do` lines 540-551 train a boosted regression forest, then run `capture grf_predict` and `assert _rc != 0`, confirming this is an error-case test only.

**Reviewer agreement:** All three reviewers agreed that fidelity testing covers only regression and causal forests. All identified the gap where `generate_reference.R` produces reference data for quantile, probability, instrumental, and survival forests that `test_fidelity.do` does not consume.

---

## 4. Fidelity Coverage Assessment

### 4.1 Fidelity Tests in `test_fidelity.do` (8 comparisons)

| # | Description | R Reference File | Method | Threshold |
|---|------------|------------------|--------|-----------|
| 1 | Regression OOB predictions | `regression_output.csv` | Pearson correlation | > 0.95 |
| 2 | Causal ATE (all) | `causal_ate.csv` | z-test | < 3 |
| 3 | Causal ATE (treated) | `causal_ate_treated.csv` | z-test | < 3 |
| 4 | Causal ATE (control) | `causal_ate_control.csv` | z-test | < 3 |
| 5 | Causal ATE (overlap) | `causal_ate_overlap.csv` | z-test | < 3 |
| 6 | Causal BLP coefficients | `causal_blp.csv` | Per-coefficient z-test | < 3 |
| 7 | Variable importance rank corr. | `causal_variable_importance.csv` | Spearman rank correlation | > 0.70 |
| 8 | DR scores correlation | `causal_scores.csv` | Pearson correlation | > 0.90 |

All tests gracefully skip if reference files are missing.

### 4.2 Reference Data Generated by `generate_reference.R` but NOT Consumed

All three reviewers identified this gap. The following CSV files are produced but have no fidelity comparison:

| Reference File | Content | Reviewer Agreement |
|---------------|---------|-------------------|
| `regression_variable_importance.csv` | Regression forest VI | Opus only |
| `causal_output.csv` | CATE predictions | Opus only |
| `causal_test_calibration.csv` | Test calibration output | All three |
| `causal_blp_hc0.csv` | BLP with HC0 vcov | Opus, Gemini |
| `causal_blp_hc3.csv` | BLP with HC3 vcov | Opus, Gemini |
| `causal_ape.csv` / `causal_ape_input.csv` | Average partial effect | All three |
| `quantile_input.csv` / `quantile_output.csv` | Quantile predictions | All three |
| `instrumental_input.csv` / `instrumental_output.csv` | Instrumental LATE | All three |
| `probability_input.csv` / `probability_output.csv` | Probability class probs | All three |
| `survival_input.csv` / `survival_output.csv` / `survival_failure_times.csv` | Survival curves | All three |
| `survival_expected.csv` | Expected survival E[T\|X] | All three |

### 4.3 Additional Fidelity Context

Codex noted that some per-forest test files (e.g., `test_regression.do`) contain optional R comparison blocks that print PASS/MARGINAL/FAILED text but do not hard-assert failure. These use a `tests/ref/` path, while `test_fidelity.do` reads from `ref/`. This path inconsistency means the two mechanisms cannot share the same reference data directory without manual adjustment. Opus and Gemini did not flag this path issue.

### 4.4 Post-Estimation Fidelity Summary

| Stata Command | R Equivalent | Fidelity Status |
|---------------|-------------|-----------------|
| `grf_ate` (all/treated/control/overlap) | `average_treatment_effect()` | TESTED (4 target samples) |
| `grf_best_linear_projection` (default) | `best_linear_projection()` | TESTED (coefficients) |
| `grf_best_linear_projection` (HC0/HC3) | `best_linear_projection(vcov.type=...)` | NOT TESTED (ref data exists) |
| `grf_variable_importance` | `variable_importance()` | TESTED (rank correlation) |
| `grf_get_scores` | `get_scores()` | TESTED (DR score correlation) |
| `grf_test_calibration` | `test_calibration()` | NOT TESTED (ref data exists) |
| `grf_average_partial_effect` | `average_partial_effect()` [E] | NOT TESTED (ref data exists) |
| `grf_rate` | `rank_average_treatment_effect()` | NOT TESTED (no ref data generated) |
| `grf_expected_survival` | Manual trapezoidal integration | NOT TESTED (ref data exists) |

**[E] Note:** R's grf has **removed** `average_partial_effect()` in favor of `average_treatment_effect()`. The Stata port retains it as a separate command. This divergence from the current R API is worth documenting.

---

## 5. R grf Functions Not Ported to Stata

The table below merges findings from all three reviewers, verified against the [R grf reference index](https://grf-labs.github.io/grf/reference/index.html).

| R Function | Purpose | Impact | Flagged By |
|-----------|---------|--------|-----------|
| `get_tree()` | Extract individual tree structure | Medium | Opus, Codex |
| `get_leaf_node()` | Return leaf node IDs for observations | Low | Opus, Codex |
| `get_forest_weights()` | Alpha/kernel weights per observation | Medium | Opus, Codex |
| `merge_forests()` | Combine multiple forest objects | Low | All three |
| `split_frequencies()` | Split variable frequency by depth | Low | Opus, Codex |
| `get_scores.multi_arm_causal_forest()` | DR scores for multi-arm forests | Medium | Codex |
| `get_scores.causal_survival_forest()` | DR scores for causal survival forests | Medium | Codex |
| `print.grf()` / `print.grf_tree()` | Formatted forest summaries | Low | Opus, Gemini |
| `plot.grf_tree()` | Tree visualization | Low | Opus, Gemini |
| `plot.rank_average_treatment_effect()` | RATE visualization | Low | None (verified) |
| `grf_options()` | Package-level options | Low | None (verified) |

**Reconciliation notes:**
- Codex listed `ll_causal_forest` and `tune_ll_causal_forest` as unported. These are **not** standalone functions in R grf's current API. Local-linear causal prediction is accessed through `predict.causal_forest()` with `linear.correction.variables`. Codex was incorrect here.
- Gemini listed `custom_forest()` as unported. This function does **not** exist in R grf's reference index. Gemini was incorrect.
- Codex listed `test_regularity`, `compute_kernel_weights`, `leaf_stats`, and `mse.convergence` as unported. These are **not** exported functions in the current grf reference index. They may be internal or from older versions. These should not be counted as missing user-facing functions.

---

## 6. R grf Options Not Available in Stata

### 6.1 Global Options Missing Across All Forest Types

All three reviewers agreed on these:

| R Option | Purpose | Impact |
|---------|---------|--------|
| `clusters` | Cluster-robust splitting and SEs | **High** |
| `equalize.cluster.weights` | Equal cluster contribution | Medium |
| `sample.weights` | Observation-level sampling weights | **High** |
| `compute.oob.predictions` | Toggle OOB predictions on/off | Low |

Gemini additionally identified:

| R Option | Purpose | Impact |
|---------|---------|--------|
| `missing.action` (MIA) | Native handling of missing covariates in splitting | **High** [F] |

**[F] Verified:** The C++ plugin (`grf_plugin.cpp`, lines 274-287) performs casewise deletion of observations with any missing values before sending data to the grf backend. The underlying C++ grf library supports MIA natively, but the Stata wrapper does not pass missing values through. This was uniquely identified by Gemini and is a significant gap.

### 6.2 Tuning Integration

| R Feature | Stata Status | Impact |
|-----------|-------------|--------|
| `tune.parameters` (integrated tuning during forest fit) | Not available; `grf_tune` is a separate command | Medium |
| `tune.num.trees`, `tune.num.reps`, `tune.num.draws` | Covered by `grf_tune` options | Low |

### 6.3 User-Supplied Nuisance Estimates

| Forest Type | R Option | Stata Status | Impact |
|------------|---------|-------------|--------|
| `causal_forest` | `W.hat`, `Y.hat` | Not available; always internally estimated | Medium |
| `causal_forest` | `orthog.boosting` | Not ported | Low |
| `instrumental_forest` | `Z.hat` | Not available; always internally estimated | Medium |
| `causal_survival_forest` | `W.hat`, `S.hat`, `C.hat`, `Y.hat` | Partially available via `numer()`/`denom()` | Medium |
| `multi_arm_causal_forest` | `W.hat` | Not available | Medium |

### 6.4 Other Per-Forest Missing Options

| Forest Type | R Option | Stata Status |
|------------|---------|-------------|
| `quantile_forest` | `regression.splitting` | Not ported |
| `ll_regression_forest` | `ll.split.variables` (variable list) | `llsplit` is boolean only, not a variable list |
| `survival_forest` | `failure.times` (arbitrary grid) | `numfailures` provides count, not arbitrary times |
| All forests | `honesty.prune.method` (enum) | `nohonestyprune` (binary toggle) |

### 6.5 Post-Estimation Missing Options

| R Function | R Option | Stata Status |
|-----------|---------|-------------|
| `best_linear_projection()` | `target.sample = "treated"/"control"` | Only "all" and "overlap" supported |
| `rank_average_treatment_effect()` | `subset` | Not available |
| `get_scores()` | Multi-arm and causal-survival forest support | Only causal and instrumental |

---

## 7. Gaps and Recommendations (Prioritized, Deduplicated)

### Priority 1: Critical

**Gap 1 -- Incomplete fidelity test consumption (all three reviewers agreed)**

`generate_reference.R` produces 22+ CSV reference files. `test_fidelity.do` consumes only 8. The following should be added:
- Quantile forest predictions vs `quantile_output.csv`
- Probability forest class probabilities vs `probability_output.csv`
- Instrumental forest LATE vs `instrumental_output.csv`
- Survival forest curves vs `survival_output.csv`
- Expected survival vs `survival_expected.csv`
- APE vs `causal_ape.csv`
- Test calibration vs `causal_test_calibration.csv`
- BLP HC0/HC3 standard errors vs `causal_blp_hc0.csv` / `causal_blp_hc3.csv`

This requires only new Stata code in `test_fidelity.do`; no R changes needed.

**Gap 2 -- No RATE fidelity data (Opus; confirmed)**

Neither `generate_reference.R` nor `test_fidelity.do` include RATE (AUTOC/QINI) reference data. R's `rank_average_treatment_effect()` should be called in `generate_reference.R` and compared in `test_fidelity.do`.

**Gap 3 -- Missing `clusters` and `sample.weights` support (all three reviewers agreed)**

These are the most significant functional gaps versus R's grf for applied research with clustered or survey data.

**Gap 4 -- Missing MIA passthrough (Gemini only; verified)**

The C++ plugin drops observations with missing covariates instead of passing them as NaN to the grf backend, which supports MIA splitting natively. For datasets with incomplete covariates, this causes silent data loss.

### Priority 2: Moderate

**Gap 5 -- No fidelity data for 6 forest types (Opus)**

`generate_reference.R` does not generate reference data for ll_regression, lm_forest, multi_arm_causal, multi_regression, boosted_regression, or causal_survival forests. At minimum, prediction-level correlations would strengthen confidence in the C++ plugin for these types.

**Gap 6 -- `grf_tune` forest-specific options untested (Codex; verified)**

`nclasses()`, `numfailures()`, `predtype()`, `reducedformweight()`, `horizon()`, `target()` are accepted by `grf_tune` but never passed with non-default values in any tune test. These options are tested in their respective forest option files, but not in the context of tuning.

**Gap 7 -- `grf_average_partial_effect` `numtreesvariance()` untested (Opus; verified)**

The option is defined in `grf_average_partial_effect.ado` (line 11) but appears in zero test files.

**Gap 8 -- `best_linear_projection` `target.sample` limited (Opus; verified)**

R supports "all", "treated", "control", "overlap". Stata only supports "all" and "overlap" (enforced at line 45 of `grf_best_linear_projection.ado`).

**Gap 9 -- `grf_get_scores` scope limited (Codex; verified)**

R's `get_scores()` supports causal, instrumental, multi_arm_causal, and causal_survival forests. Stata's `grf_get_scores` only accepts causal and instrumental (enforced at lines 12-16 of `grf_get_scores.ado`).

**Gap 10 -- Fidelity path inconsistency (Codex only; verified)**

`test_fidelity.do` reads from `ref/...` while per-forest test files read from `tests/ref/...`. `generate_reference.R` writes to `ref/` relative to working directory. This means running `test_fidelity.do` from the project root works, but the optional per-forest R comparisons expect a different path.

**Gap 11 -- `average_partial_effect` is removed from R grf (verified via web)**

R grf has removed `average_partial_effect()` in favor of `average_treatment_effect()`. The Stata port retains it as `grf_average_partial_effect`. This divergence should be documented, and consideration given to whether the Stata command should be deprecated or updated to match current R API.

### Priority 3: Low

**Gap 12 -- `in` qualifier untested for `grf_rate` and `grf_average_partial_effect` (Codex; verified)**

The `in` qualifier is not tested for these commands (though `if` is tested for both).

**Gap 13 -- Survival DGP error handling (Opus)**

`test_generate_data.do` tests error handling for `grf_generate_causal_data` (invalid DGP, minimum p) but not for `grf_generate_causal_survival_data`.

**Gap 14 -- No explicit seed-reproducibility test for forest estimation (Opus)**

No test verifies that running the same forest command twice with the same seed produces identical predictions. The fidelity tests implicitly test against a fixed R reference, but a Stata-to-Stata reproducibility test is absent.

**Gap 15 -- `ll_regression_forest` `llsplit` is boolean, not a variable list (Opus; verified)**

R's `ll.split.variables` accepts a vector of variable indices. Stata's `llsplit` is a binary toggle enabling LL splitting on all variables.

**Gap 16 -- User-supplied nuisance estimates not exposed (Opus, Codex)**

`W.hat`, `Y.hat`, `Z.hat` for causal/instrumental/multi-arm forests cannot be user-supplied. These must always be internally estimated.

---

## 8. Reviewer Agreement Summary

| Topic | Agreement |
|-------|-----------|
| Every command has tests | All three agreed |
| All 12 forest types have fit/predict/tune tests | All three agreed |
| Fidelity covers only regression + causal | All three agreed |
| Unused reference CSVs exist | All three agreed |
| `clusters` and `sample.weights` missing | All three agreed |
| `merge_forests` not ported | All three agreed |
| `get_tree`, `get_forest_weights` not ported | Opus and Codex (Gemini did not mention) |
| MIA missing-data handling not available | Gemini only (verified correct) |
| `numer()`/`denom()` untested | Codex only (verified incorrect) |
| `quantiles()`/`compliancescore()` untested | Opus only (verified incorrect -- tested in `test_options_post_estimation.do`) |
| `grf_tune` forest-specific options untested | Codex only (verified correct) |
| Fidelity path inconsistency | Codex only (verified correct) |
| Boosted regression predict is negative-path only | Codex only (verified correct) |
| `ll_causal_forest` not ported | Codex only (verified incorrect -- function does not exist in R grf) |
| `custom_forest()` not ported | Gemini only (verified incorrect -- function does not exist in R grf) |

---

## Appendix: File Inventory

**Commands:** 24 .ado files in `/tmp/grf_stata/`
**Tests:** 39 .do files + 2 .R files in `/tmp/grf_stata/tests/`
- 12 basic forest tests
- 13 per-option tests (12 forest types + 1 post-estimation omnibus)
- 8 post-estimation tests (predict, tune, tune_extended, post_estimation, get_scores, rate, average_partial_effect, expected_survival)
- 1 data generation test
- 1 fidelity test
- 4 other (features, full_pipeline, speed, debug_survival)

**Estimated total individual test cases:** ~365 (across all assertion-based .do files, excluding speed benchmarks and debug utilities)
