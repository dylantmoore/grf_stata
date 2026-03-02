# grf_stata Package Review: Test Coverage and R Fidelity

**Reviewer:** Claude Opus 4.6
**Date:** 2026-02-27
**Scope:** Read-only audit of all 24 .ado files, 39 .do test files, and 2 .R support scripts
**Package location:** `/tmp/grf_stata/`

---

## 1. Executive Summary

The grf_stata package is a Stata port of R's `grf` (Generalized Random Forests) package. It implements 12 forest types, 10 post-estimation commands, and 2 data-generation utilities through a C++ plugin architecture. The test suite comprises 39 `.do` files organized into basic forest tests, per-option exhaustive tests, post-estimation tests, fidelity tests, and integration tests.

**Strengths:**

- Every one of the 24 commands has at least one dedicated test file. Most have two (a basic test plus a comprehensive per-option test).
- All 12 forest types are tested for fit, out-of-sample predict (`test_predict.do`), and cross-validation tune (`test_tune_extended.do`).
- Per-option test files exist for all 12 forest types plus post-estimation commands (13 files total), testing every documented option individually and in combination.
- Fidelity tests compare Stata output to R reference data for regression predictions, causal ATE (4 target samples), BLP coefficients, variable importance, and DR scores.
- Integration tests (`test_full_pipeline.do`) exercise multi-step workflows end to end.

**Weaknesses:**

- Fidelity tests cover only 2 of 12 forest types (regression and causal). The R reference generator (`generate_reference.R`) produces reference data for quantile, probability, survival, and instrumental forests, but `test_fidelity.do` does not consume these reference files.
- Several R grf post-estimation outputs with generated reference data have no fidelity test: APE (`causal_ape.csv`), test_calibration (`causal_test_calibration.csv`), quantile predictions, probability predictions, instrumental LATE, survival curves, and expected survival time.
- R's `grf` package has options not ported to Stata: `clusters`, `equalize.cluster.weights`, `sample.weights`, `tune.parameters` (built-in tuning), and `compute.oob.predictions` (as a toggle).
- R functions not ported: `get_tree()`, `get_leaf_node()`, `merge_forests()`, `get_forest_weights()`, `split_frequencies()`, `print.grf()`, and `plot.grf()`.
- The RATE command (`grf_rate`) has no fidelity test against R, and no reference data is generated for it.
- `grf_test_calibration` has reference data generated (`causal_test_calibration.csv`) but no fidelity comparison in `test_fidelity.do`.

**Overall assessment:** The test suite is thorough for functional correctness. The principal gap is in fidelity testing -- only regression and causal forests are validated against R, despite reference data existing for several additional forest types and post-estimation quantities.

---

## 2. Command Coverage Matrix

### 2.1 Forest Estimation Commands (12)

| Command | Basic Test File | Options Test File | Options Tested | Options Untested |
|---------|----------------|-------------------|----------------|-----------------|
| `grf_regression_forest` | `test_regression.do` | `test_options_regression.do` (16 tests) | generate, ntrees, seed, mtry, minnodesize, samplefrac, nohonesty, honestyfrac, nohonestyprune, alpha, imbalancepenalty, cigroupsize, numthreads, estimatevariance, replace, vargenerate | None |
| `grf_causal_forest` | `test_causal.do` | `test_options_causal.do` | generate, ntrees, seed, mtry, minnodesize, samplefrac, nohonesty, honestyfrac, nohonestyprune, alpha, imbalancepenalty, cigroupsize, numthreads, estimatevariance, replace, vargenerate, nuisancetrees, nostabilizesplits, yhatgenerate, whatgenerate | None |
| `grf_quantile_forest` | `test_quantile.do` | `test_options_quantile.do` | generate, ntrees, seed, mtry, minnodesize, samplefrac, nohonesty, honestyfrac, nohonestyprune, alpha, imbalancepenalty, numthreads, replace, quantiles | None |
| `grf_probability_forest` | `test_probability.do` | `test_options_probability.do` | generate, ntrees, seed, mtry, minnodesize, samplefrac, nohonesty, honestyfrac, nohonestyprune, alpha, imbalancepenalty, numthreads, replace, nclasses | None |
| `grf_survival_forest` | `test_survival.do` | `test_options_survival.do` | generate, ntrees, seed, mtry, minnodesize, samplefrac, nohonesty, honestyfrac, nohonestyprune, alpha, imbalancepenalty, numthreads, replace, noutput, numfailures, predtype | None |
| `grf_instrumental_forest` | `test_instrumental.do` | `test_options_instrumental.do` | generate, ntrees, seed, mtry, minnodesize, samplefrac, nohonesty, honestyfrac, nohonestyprune, alpha, imbalancepenalty, numthreads, estimatevariance, replace, vargenerate, reducedformweight, nuisancetrees, stabilizesplits | None |
| `grf_ll_regression_forest` | `test_ll_regression.do` | `test_options_ll_regression.do` | generate, ntrees, seed, mtry, minnodesize, samplefrac, nohonesty, honestyfrac, nohonestyprune, alpha, imbalancepenalty, numthreads, estimatevariance, replace, vargenerate, llsplit, lllambda, llweightpenalty, llcutoff | None |
| `grf_lm_forest` | `test_lm_forest.do` | `test_options_lm_forest.do` | generate, ntrees, seed, mtry, minnodesize, samplefrac, nohonesty, honestyfrac, nohonestyprune, alpha, imbalancepenalty, numthreads, replace, xvars, nuisancetrees, nostabilizesplits | None |
| `grf_multi_arm_causal_forest` | `test_multi_arm.do` | `test_options_multi_arm.do` | generate, ntrees, seed, mtry, minnodesize, samplefrac, nohonesty, honestyfrac, nohonestyprune, alpha, imbalancepenalty, numthreads, estimatevariance, replace, vargenerate, ntreat, nostabilizesplits | None |
| `grf_multi_regression_forest` | `test_multi_regression.do` | `test_options_multi_regression.do` | generate, ntrees, seed, mtry, minnodesize, samplefrac, nohonesty, honestyfrac, nohonestyprune, alpha, imbalancepenalty, numthreads, replace, ndep | None |
| `grf_boosted_regression_forest` | `test_boosted_regression.do` | `test_options_boosted.do` | generate, ntrees, seed, mtry, minnodesize, samplefrac, nohonesty, honestyfrac, nohonestyprune, alpha, imbalancepenalty, numthreads, replace, booststeps, boostmaxsteps, boosterrorreduction, boosttreestune, nostabilizesplits | None |
| `grf_causal_survival_forest` | `test_causal_survival.do` | `test_options_causal_survival.do` | generate, ntrees, seed, mtry, minnodesize, samplefrac, nohonesty, honestyfrac, nohonestyprune, alpha, imbalancepenalty, numthreads, estimatevariance, replace, vargenerate, numer, denom, horizon, target, nostabilizesplits | None |

### 2.2 Post-Estimation Commands (10)

| Command | Test File(s) | Options Tested | Options Untested |
|---------|-------------|----------------|-----------------|
| `grf_predict` | `test_predict.do` (13 tests, all forest types) | generate, replace, numthreads (via `test_options_post_estimation.do`) | None |
| `grf_tune` | `test_tune.do`, `test_tune_extended.do` (13 tests + reproducibility) | foresttype (all 12), numreps, tunetrees, seed, nohonesty, nohonestyprune, nostabilizesplits, xvars, ntreat, ndep, nclasses, numfailures, predtype, reducedformweight, horizon, target (via `test_options_post_estimation.do`) | numthreads (not explicitly tested in isolation; used only via default) |
| `grf_ate` | `test_post_estimation.do`, `test_options_post_estimation.do` | targetsample(all/treated/control/overlap), if/in | None |
| `grf_test_calibration` | `test_post_estimation.do`, `test_options_post_estimation.do` | if/in, wrong-forest error | None (command has no options besides if/in) |
| `grf_variable_importance` | `test_post_estimation.do`, `test_options_post_estimation.do` | ntrees, seed, maxdepth | None |
| `grf_best_linear_projection` | `test_post_estimation.do`, `test_options_post_estimation.do` | varlist, vcovtype(HC0/HC1/HC2/HC3), targetsample(all/overlap) | None |
| `grf_get_scores` | `test_get_scores.do` (7 tests) | generate, replace; tests causal + instrumental forests; error without prior estimation; error with wrong forest type | None |
| `grf_rate` | `test_rate.do` (3 tests) | target(AUTOC/QINI), bootstrap, catevar, seed | quantiles, compliancescore |
| `grf_average_partial_effect` | `test_average_partial_effect.do` (7 tests) | if/in, debiasweights, nocalibrate, wrong-forest error, binary W vs continuous W, CI properties | numtreesvariance (not explicitly tested) |
| `grf_expected_survival` | `test_expected_survival.do` (7 tests) | generate, replace, predictions, grid; wrong-forest error; positive values; heterogeneity | None |

### 2.3 Data Generation Commands (2)

| Command | Test File | DGPs Tested | Error Cases Tested |
|---------|----------|-------------|-------------------|
| `grf_generate_causal_data` | `test_generate_data.do` | All 11: simple, aw1, aw2, aw3, ai1, ai2, kunzel, nw1, nw2, nw3, nw4 | minimum p enforcement (p=2 for nw1, p=2 for simple), invalid DGP name |
| `grf_generate_causal_survival_data` | `test_generate_data.do` | All 6: simple1, type1, type2, type3, type4, type5 | (none explicitly, but invalid DGP would be caught by validation) |

**Untested data generation error case:** Invalid DGP name for `grf_generate_causal_survival_data` (only tested for `grf_generate_causal_data`). Minimum p enforcement for survival DGPs is also not tested.

---

## 3. Forest Type Coverage Matrix

| Forest Type | Fit Test | Predict Test | Tune Test | Fidelity Test | Options Test |
|------------|----------|-------------|-----------|---------------|-------------|
| regression | `test_regression.do` | `test_predict.do` (Test 1) | `test_tune_extended.do` (Test 1) | `test_fidelity.do` (Test 1: correlation > 0.95) | `test_options_regression.do` |
| causal | `test_causal.do` | `test_predict.do` (Test 2) | `test_tune_extended.do` (Test 2) | `test_fidelity.do` (Tests 2-8: ATE, BLP, VI, DR scores) | `test_options_causal.do` |
| quantile | `test_quantile.do` | `test_predict.do` (Test 3) | `test_tune_extended.do` (Test 3) | **MISSING** (ref data exists: `quantile_output.csv`) | `test_options_quantile.do` |
| probability | `test_probability.do` | `test_predict.do` (Test 4) | `test_tune_extended.do` (Test 4) | **MISSING** (ref data exists: `probability_output.csv`) | `test_options_probability.do` |
| instrumental | `test_instrumental.do` | `test_predict.do` (Test 5) | `test_tune_extended.do` (Test 5) | **MISSING** (ref data exists: `instrumental_output.csv`) | `test_options_instrumental.do` |
| survival | `test_survival.do` | `test_predict.do` (Test 6) | `test_tune_extended.do` (Test 6) | **MISSING** (ref data exists: `survival_output.csv`, `survival_expected.csv`) | `test_options_survival.do` |
| ll_regression | `test_ll_regression.do` | `test_predict.do` (Test 7) | `test_tune_extended.do` (Test 7) | **MISSING** (no ref data generated) | `test_options_ll_regression.do` |
| lm_forest | `test_lm_forest.do` | `test_predict.do` (Test 8) | `test_tune_extended.do` (Test 8) | **MISSING** (no ref data generated) | `test_options_lm_forest.do` |
| multi_arm_causal | `test_multi_arm.do` | `test_predict.do` (Test 9) | `test_tune_extended.do` (Test 9) | **MISSING** (no ref data generated) | `test_options_multi_arm.do` |
| multi_regression | `test_multi_regression.do` | `test_predict.do` (Test 10) | `test_tune_extended.do` (Test 10) | **MISSING** (no ref data generated) | `test_options_multi_regression.do` |
| boosted_regression | `test_boosted_regression.do` | `test_predict.do` (Test 11) | `test_tune_extended.do` (Test 12) | **MISSING** (no ref data generated) | `test_options_boosted.do` |
| causal_survival | `test_causal_survival.do` | `test_predict.do` (Test 12) | `test_tune_extended.do` (Test 11) | **MISSING** (no ref data generated) | `test_options_causal_survival.do` |

**Summary:** All 12 forest types have complete functional test coverage (fit, predict, tune, options). Fidelity against R is tested for 2 of 12 forest types (regression and causal). Reference data exists for 4 additional forest types (quantile, probability, instrumental, survival) but is not consumed by `test_fidelity.do`.

---

## 4. Fidelity Coverage Assessment

### 4.1 Current Fidelity Tests (`test_fidelity.do`)

| Test # | Description | R Reference File | Method | Threshold |
|--------|------------|------------------|--------|-----------|
| 1 | Regression OOB predictions | `regression_output.csv` | Pearson correlation | > 0.95 |
| 2 | Causal ATE (all) | `causal_ate.csv` | z-test: \|Stata - R\| / sqrt(SE_S^2 + SE_R^2) | < 3 |
| 3 | Causal ATE (treated) | `causal_ate_treated.csv` | z-test | < 3 |
| 4 | Causal ATE (control) | `causal_ate_control.csv` | z-test | < 3 |
| 5 | Causal ATE (overlap) | `causal_ate_overlap.csv` | z-test | < 3 |
| 6 | Causal BLP coefficients | `causal_blp.csv` | Per-coefficient z-test | < 3 (all 6 coefficients) |
| 7 | Variable importance rank correlation | `causal_variable_importance.csv` | Spearman rank correlation | > 0.70 |
| 8 | DR scores correlation | `causal_scores.csv` | Pearson correlation | > 0.90 |

All 8 fidelity tests gracefully skip if reference files are missing (lines 17-19 pattern in `test_fidelity.do`).

### 4.2 Reference Data Generated but NOT Tested

The following reference CSV files are produced by `generate_reference.R` but have no corresponding comparison in `test_fidelity.do`:

| R Reference File | Content | Line in `generate_reference.R` |
|-----------------|---------|-------------------------------|
| `regression_variable_importance.csv` | Regression forest VI | Line 68 |
| `causal_output.csv` | Causal CATE predictions | Line 107 |
| `causal_test_calibration.csv` | Test calibration output | Line 140 |
| `causal_blp_hc0.csv` | BLP with HC0 vcov | Line 364 |
| `causal_blp_hc3.csv` | BLP with HC3 vcov | Line 376 |
| `causal_ape.csv` | Average partial effect | Line 347 |
| `causal_ape_input.csv` | APE input data | Line 336 |
| `quantile_input.csv` / `quantile_output.csv` | Quantile forest predictions | Lines 158, 176 |
| `instrumental_input.csv` / `instrumental_output.csv` | Instrumental forest LATE | Lines 196, 213 |
| `probability_input.csv` / `probability_output.csv` | Probability forest class probs | Lines 231, 246 |
| `survival_input.csv` / `survival_output.csv` / `survival_failure_times.csv` | Survival curves | Lines 264-290 |
| `survival_expected.csv` | Expected survival time E[T\|X] | Line 416 |

### 4.3 Post-Estimation Commands Without Fidelity Tests

| Stata Command | R Equivalent | Fidelity Status |
|---------------|-------------|-----------------|
| `grf_ate` (all) | `average_treatment_effect(cf)` | TESTED |
| `grf_ate` (treated/control/overlap) | `average_treatment_effect(cf, target.sample=...)` | TESTED |
| `grf_best_linear_projection` (default HC3) | `best_linear_projection(cf, X)` | TESTED (BLP coefficients) |
| `grf_best_linear_projection` (HC0) | `best_linear_projection(cf, X, vcov.type="HC0")` | NOT TESTED (ref data exists) |
| `grf_best_linear_projection` (HC3) | `best_linear_projection(cf, X, vcov.type="HC3")` | NOT TESTED (ref data exists) |
| `grf_variable_importance` | `variable_importance(rf)` | TESTED (rank correlation) |
| `grf_get_scores` | `get_scores(cf)` | TESTED (correlation > 0.90) |
| `grf_test_calibration` | `test_calibration(cf)` | NOT TESTED (ref data exists) |
| `grf_average_partial_effect` | `average_partial_effect(cf_cont)` | NOT TESTED (ref data exists) |
| `grf_rate` | `rank_average_treatment_effect()` | NOT TESTED (no ref data) |
| `grf_expected_survival` | (manual trapezoidal integration) | NOT TESTED (ref data exists) |

---

## 5. R grf Functions Not Ported to Stata

The following user-facing R grf functions have no Stata equivalent:

| R Function | Purpose | Impact |
|-----------|---------|--------|
| `get_tree(forest, i)` | Extract single tree structure (split vars, split values, leaf values) | Medium -- useful for interpretation and debugging |
| `get_leaf_node(forest, X)` | Return leaf node IDs for observations | Low -- primarily for advanced analysis |
| `merge_forests(f1, f2, ...)` | Combine multiple forests (e.g., from distributed training) | Low -- niche use case |
| `get_forest_weights(forest, X)` | Return alpha weights (kernel weights) for each observation | Medium -- enables custom post-estimation |
| `split_frequencies(forest)` | Tabulate split variable frequencies per depth | Low -- variable_importance provides similar info |
| `print.grf()` / `summary()` | Formatted summary of forest parameters and fit | Low -- .ado files print headers at estimation time |
| `plot.grf()` | Visualization of single trees | Low -- tree plotting is uncommon in applied work |
| `tune_forest()` / auto-tuning via `tune.parameters` | Built-in cross-validation tuning during forest fitting | Medium -- Stata has separate `grf_tune` command but not integrated tuning |
| `predict(forest, newdata, estimate.variance=TRUE)` for quantile/probability/survival | Variance estimation for non-regression/causal forests | Low -- R's grf does not support this either for these types |

**Note:** The most impactful missing functions are `get_tree()`, `get_forest_weights()`, and integrated `tune.parameters`. The separate `grf_tune` command partially compensates for the last of these, though it requires the user to manually pass tuned parameters to the forest command.

---

## 6. R grf Options Not Available in Stata

### 6.1 Global Options Missing Across All Forest Types

| R Option | Purpose | Impact |
|---------|---------|--------|
| `clusters` | Cluster-robust splitting and standard errors | **High** -- important for clustered data (panel, survey) |
| `equalize.cluster.weights` | Reweight so each cluster contributes equally | Medium -- companion to `clusters` |
| `sample.weights` | Observation-level sampling weights | **High** -- important for survey/weighted data |
| `tune.parameters` | Character vector of parameters to auto-tune during fitting | Medium -- `grf_tune` exists but is separate |
| `tune.num.trees`, `tune.num.reps`, `tune.num.draws` | Tuning hyperparameters (when using built-in tuning) | Low -- covered by `grf_tune` |
| `compute.oob.predictions` | Toggle OOB predictions on/off | Low -- Stata always computes OOB predictions |

### 6.2 Per-Forest-Type Missing Options

| Forest Type | R Option | Stata Status |
|------------|---------|-------------|
| `regression_forest` | `ci.group.size` as default 2 when `compute.oob.predictions=TRUE` | Stata defaults to 1; auto-bumps to 2 only with `estimatevariance` |
| `causal_forest` | `W.hat` (user-supplied propensity) | Not available; always internally estimated via regression forest |
| `causal_forest` | `Y.hat` (user-supplied outcome model) | Not available; always internally estimated |
| `causal_forest` | `orthog.boosting` | Not ported |
| `instrumental_forest` | `Z.hat` (user-supplied instrument model) | Not available; always internally estimated |
| `survival_forest` | `failure.times` (user-supplied grid) | Stata has `numfailures` but not arbitrary failure time grid |
| `causal_survival_forest` | `W.hat`, `S.hat`, `C.hat`, `Y.hat` (user-supplied nuisance) | Partially available via `numer()`/`denom()` but different interface |
| `probability_forest` | — | All main options ported |
| `quantile_forest` | `regression.splitting` (toggle regression vs quantile splitting) | Not ported |
| `ll_regression_forest` | `ll.split.variables` (subset of variables for LL splits) | Stata has `llsplit` as a boolean toggle, not a variable list |
| `multi_arm_causal_forest` | `W.hat` (user-supplied multinomial propensity) | Not available |
| `boosted_regression_forest` | `tune.parameters` | Not integrated; must use `grf_tune` separately |
| All forests | `honesty.prune.leaves` vs `honesty.prune.method` | Stata uses `nohonestyprune` (binary), R added `honesty.prune.method` in later versions |

### 6.3 Post-Estimation Missing Options

| R Function | R Option | Stata Status |
|-----------|---------|-------------|
| `average_treatment_effect()` | `subset` (observation subset) | Stata uses `if`/`in` instead (equivalent) |
| `rank_average_treatment_effect()` | `target = "QINI"` with subset | Stata `grf_rate` supports target(QINI) but not a subset option |
| `rank_average_treatment_effect()` | `R` (number of bootstrap reps) | Stata has `bootstrap()` (equivalent) |
| `best_linear_projection()` | `target.sample = "treated"/"control"` | Stata only supports `all` and `overlap` (line 45 of `grf_best_linear_projection.ado`) |
| `average_partial_effect()` | `subset` | Stata uses `if`/`in` (equivalent) |

---

## 7. Gaps and Recommendations

### 7.1 Critical Gaps

**Gap 1: Incomplete fidelity test consumption.**
`generate_reference.R` produces 22+ CSV reference files (lines 23-428), but `test_fidelity.do` only consumes 8 of them. The following fidelity tests should be added to `test_fidelity.do`:

- Quantile predictions: correlate Stata `grf_quantile_forest` output with `quantile_output.csv` (per-quantile correlation > 0.90).
- Probability predictions: correlate Stata `grf_probability_forest` class probabilities with `probability_output.csv`.
- Instrumental LATE: correlate Stata `grf_instrumental_forest` with `instrumental_output.csv`, or z-test on mean LATE.
- Survival curves: correlate Stata `grf_survival_forest` curves with `survival_output.csv`.
- Expected survival: compare Stata `grf_expected_survival` with `survival_expected.csv`.
- APE: z-test comparing Stata `grf_average_partial_effect` estimate with `causal_ape.csv`.
- Test calibration: compare Stata `grf_test_calibration` coefficients with `causal_test_calibration.csv`.
- BLP HC0/HC3: compare Stata BLP standard errors with `causal_blp_hc0.csv` and `causal_blp_hc3.csv`.

**Gap 2: No RATE fidelity test.**
Neither `generate_reference.R` nor `test_fidelity.do` include RATE (AUTOC/QINI) reference data. R's `rank_average_treatment_effect()` should be called in `generate_reference.R` and the output compared in `test_fidelity.do`.

**Gap 3: Missing cluster/weight support.**
R's `clusters` and `sample.weights` arguments are important for applied research with clustered or survey data. These are not ported to any Stata command. This is the most significant functional gap versus R.

### 7.2 Moderate Gaps

**Gap 4: `grf_rate` option coverage.**
`test_rate.do` has only 3 tests. The `quantiles()` and `compliancescore()` options are not tested. The `compliancescore` option is particularly important for instrumental variable settings.

File: `grf_rate.ado` line 15 defines `Quantiles(numlist >0 <=1)`.
File: `grf_rate.ado` line 18 defines `COMPliancescore(varname numeric)`.
File: `test_rate.do` -- neither option appears in any test.

**Gap 5: `grf_average_partial_effect` `numtreesvariance()` not tested.**
File: `grf_average_partial_effect.ado` line 11 defines `NUMTreesvariance(integer 500)`.
File: `test_average_partial_effect.do` -- this option does not appear in any test call.

**Gap 6: `grf_tune` `numthreads()` not explicitly tested.**
While `numthreads()` is exercised in options tests for forest commands, `test_tune.do` and `test_tune_extended.do` never pass a non-default `numthreads()` value.

**Gap 7: `best_linear_projection` target.sample limited.**
R supports `target.sample = c("all", "treated", "control", "overlap")` for BLP. Stata only supports `all` and `overlap` (line 45 of `grf_best_linear_projection.ado`). The `treated` and `control` options are missing.

**Gap 8: No fidelity data for ll_regression, lm_forest, multi_arm, multi_regression, boosted_regression, or causal_survival.**
`generate_reference.R` does not generate reference data for these 6 forest types. Adding at least prediction-level reference comparisons would strengthen confidence in the C++ plugin's correctness for these types.

### 7.3 Minor Gaps

**Gap 9: Survival DGP error handling.**
`test_generate_data.do` tests error handling (invalid DGP name, minimum p) for `grf_generate_causal_data` but not for `grf_generate_causal_survival_data`. Adding parallel error-case tests for survival data generation would improve completeness.

**Gap 10: `grf_predict` numthreads test.**
`test_predict.do` tests all 12 forest types but always uses default `numthreads`. The `test_options_post_estimation.do` file does test `numthreads` for predict, so this is covered but in a separate file.

**Gap 11: Seed reproducibility coverage.**
`test_tune_extended.do` Test 13 verifies tuning reproducibility with matching seeds. `test_generate_data.do` Tests 12 and 22 verify data generation reproducibility. However, no test verifies that forest estimation itself (e.g., `grf_regression_forest`) produces identical predictions across runs with the same seed. (The fidelity tests implicitly test this by comparing to a fixed R reference, but an explicit Stata-to-Stata reproducibility test is absent.)

**Gap 12: `ll_regression_forest` `llsplit` is boolean, not variable list.**
R's `ll.split.variables` accepts a vector of variable indices specifying which variables to use for local-linear splitting. Stata's `llsplit` (line in `grf_ll_regression_forest.ado`) is a boolean toggle that enables LL splitting on all variables. This limits functionality compared to R.

### 7.4 Recommended Priority Actions

1. **High priority:** Add fidelity tests to `test_fidelity.do` for the 14+ reference CSV files already generated by `generate_reference.R`. This requires only Stata code, no R changes.

2. **High priority:** Add RATE reference data generation to `generate_reference.R` and a corresponding fidelity test.

3. **High priority:** Port `clusters` and `sample.weights` options to the C++ plugin and all forest .ado files.

4. **Medium priority:** Add fidelity reference data for ll_regression, lm_forest, multi_arm, multi_regression, boosted_regression, and causal_survival forests.

5. **Medium priority:** Test `grf_rate` `compliancescore()` and `quantiles()` options.

6. **Medium priority:** Add `target.sample = "treated"/"control"` support to `grf_best_linear_projection`.

7. **Low priority:** Port `get_tree()` and `get_forest_weights()` for advanced users.

8. **Low priority:** Add explicit seed-reproducibility tests for forest estimation commands.

---

## Appendix A: File Inventory

### A.1 Command Files (24 .ado files)

```
/tmp/grf_stata/grf_ate.ado                         (151 lines)
/tmp/grf_stata/grf_average_partial_effect.ado       (180 lines)
/tmp/grf_stata/grf_best_linear_projection.ado       (271 lines)
/tmp/grf_stata/grf_boosted_regression_forest.ado    (252 lines)
/tmp/grf_stata/grf_causal_forest.ado                (331 lines)
/tmp/grf_stata/grf_causal_survival_forest.ado       (418 lines)
/tmp/grf_stata/grf_expected_survival.ado            (235 lines)
/tmp/grf_stata/grf_generate_causal_data.ado
/tmp/grf_stata/grf_generate_causal_survival_data.ado
/tmp/grf_stata/grf_get_scores.ado                   (155 lines)
/tmp/grf_stata/grf_instrumental_forest.ado          (320 lines)
/tmp/grf_stata/grf_ll_regression_forest.ado         (225 lines)
/tmp/grf_stata/grf_lm_forest.ado                    (348 lines)
/tmp/grf_stata/grf_multi_arm_causal_forest.ado      (360 lines)
/tmp/grf_stata/grf_multi_regression_forest.ado      (234 lines)
/tmp/grf_stata/grf_predict.ado                      (~560 lines)
/tmp/grf_stata/grf_probability_forest.ado           (188 lines)
/tmp/grf_stata/grf_quantile_forest.ado              (205 lines)
/tmp/grf_stata/grf_rate.ado                         (414 lines)
/tmp/grf_stata/grf_regression_forest.ado            (193 lines)
/tmp/grf_stata/grf_survival_forest.ado              (222 lines)
/tmp/grf_stata/grf_test_calibration.ado             (136 lines)
/tmp/grf_stata/grf_tune.ado                         (~1900 lines)
/tmp/grf_stata/grf_variable_importance.ado          (136 lines)
```

### A.2 Test Files (39 .do files + 2 .R files)

```
Basic forest tests (12):
  tests/test_regression.do
  tests/test_causal.do
  tests/test_quantile.do
  tests/test_probability.do
  tests/test_survival.do
  tests/test_instrumental.do
  tests/test_ll_regression.do
  tests/test_lm_forest.do
  tests/test_multi_arm.do
  tests/test_multi_regression.do
  tests/test_boosted_regression.do
  tests/test_causal_survival.do

Per-option tests (13):
  tests/test_options_regression.do          (307 lines, 16 tests)
  tests/test_options_causal.do              (352 lines)
  tests/test_options_instrumental.do        (347 lines)
  tests/test_options_ll_regression.do       (358 lines)
  tests/test_options_lm_forest.do           (335 lines)
  tests/test_options_multi_arm.do           (320 lines)
  tests/test_options_multi_regression.do    (294 lines)
  tests/test_options_boosted.do             (393 lines)
  tests/test_options_probability.do         (312 lines)
  tests/test_options_quantile.do            (278 lines)
  tests/test_options_survival.do            (362 lines)
  tests/test_options_causal_survival.do     (366 lines)
  tests/test_options_post_estimation.do     (1542 lines)

Post-estimation tests (7):
  tests/test_post_estimation.do
  tests/test_predict.do                     (13 forest-type tests)
  tests/test_tune.do
  tests/test_tune_extended.do               (13 tests + reproducibility)
  tests/test_get_scores.do                  (7 tests)
  tests/test_rate.do                        (3 tests)
  tests/test_average_partial_effect.do      (7 tests)
  tests/test_expected_survival.do           (7 tests)

Data generation tests (1):
  tests/test_generate_data.do               (22+ tests)

Fidelity tests (1):
  tests/test_fidelity.do                    (8 fidelity comparisons)

Integration / other (4):
  tests/test_features.do                    (5 cross-cutting feature tests)
  tests/test_full_pipeline.do               (end-to-end workflow tests)
  tests/test_speed.do                       (performance benchmarks)
  tests/debug_survival.do                   (debugging utility)

R support scripts (2):
  tests/generate_reference.R                (generates 22+ reference CSVs)
  tests/benchmark_r.R                       (R timing benchmarks)
```

---

## Appendix B: Test Count Summary

| Category | Files | Approximate Test Count |
|----------|-------|----------------------|
| Basic forest tests | 12 | ~60 |
| Per-option tests | 13 | ~200 |
| Post-estimation tests | 8 | ~60 |
| Data generation tests | 1 | ~22 |
| Fidelity tests | 1 | 8 |
| Integration/features | 3 | ~15 |
| **Total** | **38** | **~365** |

(Excludes `debug_survival.do` and `test_speed.do` which are utilities, not assertion-based tests.)
