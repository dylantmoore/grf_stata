# Consolidated Fidelity Discrepancy Report: `grf_stata` vs R `grf` 2.5.0

**Date:** 2026-02-28
**Scope:** 15 independent test suites, 295 total test cases
**R:** grf 2.5.0 (R 4.5.2)
**Stata:** StataNow 19.5 MP, grf_stata v0.1.0–v0.3.0

---

## Aggregate Results

| Suite | Command(s) | Tests | Pass | Fail | Pass Rate |
|-------|-----------|:-----:|:----:|:----:|:---------:|
| 01 | regression_forest | 21 | 21 | 0 | 100% |
| 02 | causal_forest (fitting) | 19 | 18 | 1 | 95% |
| 03 | average_treatment_effect | 19 | 16 | 3 | 84% |
| 04 | instrumental_forest | 18 | 17 | 1 | 94% |
| 05 | quantile_forest | 17 | 15 | 2 | 88% |
| 06 | probability_forest | 13 | 13 | 0 | 100% |
| 07 | survival_forest + expected_survival | 17 | 7 | 10 | 41%* |
| 08 | causal_survival_forest | 16 | 3 | 13 | 19% |
| 09 | multi_arm_causal_forest | 36 | 33 | 3 | 92% |
| 10 | boosted_regression_forest | 17 | 17 | 0 | 100% |
| 11 | ll_regression_forest | 18 | 18 | 0 | 100% |
| 12 | lm_forest | 17 | 17 | 0 | 100% |
| 13 | multi_regression_forest | 14 | 14 | 0 | 100% |
| 14 | BLP + calibration + get_scores | 29 | 27 | 2 | 93% |
| 15 | variable_importance + RATE + tune | 24 | 14 | 10 | 58% |
| **Total** | | **295** | **250** | **45** | **85%** |

*Suite 07: 9 of 10 "failures" are a narrow-variance column-2 artifact (see D-14); 1 is an API design difference (see D-13).

---

## Discrepancy Index

| ID | Category | Severity | Component | Summary |
|----|----------|----------|-----------|---------|
| D-1 | Bug | HIGH | grf_ate.ado | TMLE + clusters not blocked |
| D-2 | Bug | HIGH | grf_best_linear_projection.ado | debiasing.weights formula wrong |
| D-3 | Bug | MEDIUM | grf_ate.ado | debiasing.weights SE formula wrong |
| D-4 | Bug | HIGH | grf_causal_survival_forest.ado | cluster() not propagated to nuisance |
| D-5 | Architectural | HIGH | grf_causal_survival_forest.ado | Simplified IPCW pipeline |
| D-6 | Architectural | MEDIUM | grf_variable_importance.ado | Re-trains forest vs reading existing |
| D-7 | Architectural | LOW | All forests | Variance estimate CI group RNG divergence |
| D-8 | API Diff | LOW | causal/instrumental/lm/multi-arm | Partial nuisance supply not supported |
| D-9 | API Diff | INFO | Multiple forests | nuisancetrees() is Stata-only |
| D-10 | API Diff | INFO | quantile_forest, ll_regression_forest | weights() is Stata extension |
| D-11 | API Diff | INFO | grf_best_linear_projection | target.sample treated/control is Stata extension |
| D-12 | API Diff | LOW | grf_survival_forest | Default fast.logrank differs |
| D-13 | API Diff | LOW | grf_survival_forest | numfailures(K) vs failure.times design |
| D-14 | Statistical | INFO | grf_survival_forest | Column-2 near-constant artifact |
| D-15 | API Diff | LOW | grf_multi_arm_causal_forest | W.hat dimension convention differs |
| D-16 | API Diff | INFO | grf_regression_forest | nomia has no R equivalent in grf 2.5.0 |
| D-17 | API Diff | LOW | All forests | Default ci.group.size differs (R=2, Stata=1) |
| D-18 | API Diff | INFO | grf_ate.ado | TMLE + debiasing: Stata correctly errors, R silently ignores |
| D-19 | Naming | LOW | grf_tune, grf_regression_forest | e(min_node) vs option name minnodesize() |
| D-20 | Behavioral | LOW | grf_rate.ado | AUTOC ~8-10% underestimate vs R |
| D-21 | Statistical | INFO | grf_causal_forest | Null effect (tau=0) low correlation expected |
| D-22 | Statistical | INFO | grf_quantile_forest | Extreme quantiles (0.01/0.99) low correlation |
| D-23 | Behavioral | LOW | grf_boosted_regression_forest | Auto-tune step count may differ |
| D-24 | Behavioral | INFO | grf_tune | Tuned parameters differ (stochastic search) |
| D-25 | API Diff | INFO | grf_test_calibration | One-sided vs two-sided p-values |
| D-26 | API Diff | INFO | grf_test_calibration | Coefficient scaling convention differs |

---

## Category 1: Bugs (Require Code Fixes)

### D-1: TMLE + Clusters Not Blocked in Stata [HIGH]

**Suite:** 03 (ATE), Test 16
**File:** `grf_ate.ado`

**Problem:** R's `grf::average_treatment_effect()` raises an error when `method = "TMLE"` is called on a forest fitted with `clusters`:

```
Error: TMLE has not yet been implemented with clustered observations.
```

Stata's `grf_ate` proceeds silently and returns results. The TMLE clever covariate computation does not account for cluster structure, so the returned estimates may have incorrect standard errors.

**R behavior:** Error (correct)
**Stata behavior:** Runs silently (incorrect)

**Recommended fix:**
```stata
if "`method'" == "TMLE" & "`e(cluster_var)'" != "" {
    display as error "TMLE has not been implemented for clustered observations"
    exit 198
}
```

---

### D-2: BLP debiasing.weights Formula Mismatch [HIGH]

**Suite:** 14 (BLP/Calibration/Scores), Tests 15-16
**File:** `grf_best_linear_projection.ado`

**Problem:** R and Stata apply debiasing weights differently in the DR score construction for BLP:

| Step | R formula | Stata formula |
|------|-----------|---------------|
| DR score | `tau_hat + debiasing_weight * (Y_resid) / var_term` | `(standard_DR_score) * debiasingweights` |

R replaces the propensity-based weight `(W-W_hat)/Var(W-W_hat)` with the user-supplied `debiasing_weight` inside the DR score numerator. Stata multiplies the already-computed standard DR score by the weight as a post-hoc scaling.

**Impact:** Large coefficient differences:

| Coefficient | R | Stata | Difference |
|-------------|---|-------|-----------|
| Intercept | 0.122 | 0.191 | +56% |
| X1 | 0.964 | 1.251 | +30% |
| X2 | 0.874 | 1.223 | +40% |

**Recommended fix:** Modify Stata's DR score computation to match R's formula: replace the `(W-W_hat)/Var(W-W_hat)` weight with the user-supplied debiasing weight before computing the DR score, rather than multiplying the final score.

---

### D-3: ATE debiasing.weights SE Formula Differs [MEDIUM]

**Suite:** 03 (ATE), Test 5
**File:** `grf_ate.ado`

**Problem:** When `debiasingweights()` is specified, Stata uses a different standard error formula than R:

| | Formula | SE value |
|---|---------|---------|
| R | `sd(weighted_scores) / sqrt(n)` | 0.024 |
| Stata | `sqrt(sum(w^2 * (score - ATE)^2)) / sum(w)` | 0.056 |

The point estimates agree (z-test passes), but the SE ratio is 2.3x, exceeding the 50% threshold.

**Impact:** Confidence intervals with debiasing weights are too wide in Stata (conservative).

---

### D-4: Causal Survival Forest cluster() Not Propagated [HIGH]

**Suite:** 08 (Causal Survival Forest), Test 6
**File:** `grf_causal_survival_forest.ado`

**Problem:** When `cluster()` is specified, the cluster variable is passed to the main causal forest call but is NOT propagated to the internal nuisance regression forests (Y.hat, W.hat estimation). This causes the nuisance forests to ignore cluster structure, producing incorrect orthogonalization.

**Evidence:** Test 6 (cluster) achieves r = 0.12, far below the 0.85 threshold — the worst single-test result across all 295 tests.

**Recommended fix:** Pass the cluster variable to all internal `grf_regression_forest` calls in the nuisance estimation stage.

---

## Category 2: Architectural Gaps (Design Differences)

### D-5: Causal Survival Forest Simplified IPCW Pipeline [HIGH]

**Suite:** 08, 13/16 tests fail
**File:** `grf_causal_survival_forest.ado`

**Problem:** R's `causal_survival_forest()` (Cui et al. 2023) uses a full survival-forest-based IPCW nuisance pipeline:

1. Fit `survival_forest` for S.hat(t|X) — conditional survival
2. Fit `survival_forest` for C.hat(t|X) — conditional censoring
3. Compute IPCW weights using S.hat and C.hat evaluated at observation-specific times
4. Pass IPCW-transformed outcomes to causal forest

Stata uses a simplified regression-forest-based pipeline:
1. Fit `regression_forest` for Y.hat (outcome regression)
2. Fit `regression_forest` for W.hat (propensity)
3. Apply fixed Kaplan-Meier censoring weights
4. Pass to causal forest

**Impact:**
- Stata predictions are compressed (SD consistently 1.3-3.2x lower than R)
- Pearson correlations range 0.51-0.93 (only 3/16 tests achieve r > 0.85)
- Correlation improves with longer horizons (more complete data)
- Tests with `target = "survival.probability"` or long horizons pass

**This is the single largest fidelity gap in the package.**

---

### D-6: Variable Importance Re-Trains Forest [MEDIUM]

**Suite:** 15, Tests 1-9 (8/9 fail Spearman threshold)
**File:** `grf_variable_importance.ado`

**Problem:** R's `variable_importance()` reads split frequencies from an already-fitted forest object in memory. Stata's `grf_variable_importance` trains a brand-new forest from scratch via the C++ plugin, since the Stata plugin architecture doesn't persist fitted forest objects.

**Impact:**
- Spearman rank correlations: 0.42-0.69 (threshold 0.70) — 8/9 fail
- But: top-2 variable identification is correct in every test (x1 and x2 always ranked #1 and #2)
- Directional sensitivity to decay.exponent and max.depth is correct
- Irrelevant variables (x3-x10) have noisy rankings because their VI scores are O(0.01)

**Assessment:** Functionally correct for the primary use case (identifying important variables). Magnitude comparison with R is inherently limited by the architectural constraint.

---

### D-7: Variance Estimate CI Group RNG Divergence [LOW]

**Suites:** 01, 02, 04, 09 (all forests with variance estimation)

**Problem:** Per-observation variance estimates show low Pearson correlation between R and Stata (r = 0.04-0.28) despite having matching distributions (mean ratios within 10%). This affects regression, causal, instrumental, multi-arm, and lm forests.

**Root cause:** Variance estimation uses leave-one-CI-group-out jackknife. The CI group assignments are determined by an internal RNG path that diverges between R's Rcpp wrapper and Stata's plugin wrapper, even with the same seed. This means different observations are grouped together, producing different per-observation variance estimates.

**Evidence:** Within-R cross-seed variance correlation (seed 42 vs seed 99) is similarly low (r = 0.17-0.24), confirming this is inherent jackknife noise, not an implementation bug.

**Impact:** Individual variance estimates are not comparable between R and Stata. Aggregate variance statistics (means, distributions) match well.

---

## Category 3: API Differences (Documented, Not Bugs)

### D-8: Partial Nuisance Supply Not Supported [LOW]

**Suites:** 02, 04, 09, 12
**Files:** `grf_causal_forest.ado`, `grf_instrumental_forest.ado`, `grf_lm_forest.ado`, `grf_multi_arm_causal_forest.ado`

R allows supplying only `Y.hat` or only `W.hat` (or only `Z.hat` for IV); the other nuisance estimates are computed internally. Stata requires either all nuisance inputs or none — `yhatinput()` and `whatinput()` must both be specified.

For instrumental forest, all three (`yhatinput`, `whatinput`, `zhatinput`) must be provided together.

---

### D-9: nuisancetrees() Is Stata-Only [INFO]

**Suites:** 02, 12
**Files:** All orthogonalized forest .ado files

R's grf 2.5.0 does not expose a `nuisance.trees` argument. Stata's `nuisancetrees()` option controls the tree count for internal nuisance regression forests. No R equivalent exists for direct comparison.

---

### D-10: weights() Is Stata Extension for Some Forests [INFO]

**Suites:** 05, 11

- `quantile_forest`: R's grf 2.5.0 does not support `sample.weights` for quantile forests. Stata accepts `weights()` and passes it to the C++ plugin.
- `ll_regression_forest`: R raises an error ("sample.weights are currently not supported for local linear forests"). Stata supports `weights()` via the C++ plugin's `weight_col_idx`.

---

### D-11: BLP target.sample=treated/control Is Stata Extension [INFO]

**Suite:** 14, Tests 17-18
**File:** `grf_best_linear_projection.ado`

R's grf 2.5.0 supports only `target.sample = c("all", "overlap")`. Stata extends this with `targetsample(treated)` (W.hat-weighted WLS) and `targetsample(control)` ((1-W.hat)-weighted WLS).

---

### D-12: Default fast.logrank Differs [LOW]

**Suite:** 07
**File:** `grf_survival_forest.ado`

R grf 2.5.0 defaults to `fast.logrank = FALSE` (standard log-rank splitting). Stata defaults to fast log-rank ON. This causes minor fidelity differences at the earliest survival time columns.

**Workaround:** Pass `nofastlogrank` in Stata to match R defaults.

---

### D-13: numfailures(K) vs failure.times Design Difference [LOW]

**Suite:** 07, Test 16
**File:** `grf_survival_forest.ado`

R's `failure.times` accepts a vector of specific time values. Stata's `numfailures(K)` accepts an integer and selects K evenly-spaced failure times via quantile indices:

```c
for (int k = 0; k < K; k++) {
    int idx = (int)((double)k / K * num_failures);
    subsampled.push_back(failure_times_vec[idx]);
}
```

When R uses `failure.times = sort(unique(Y[D==1]))[1:K]` (first K consecutive times) and Stata uses `numfailures(K)` (K evenly-spaced times), column j refers to different time points. Correlation drops to near-zero on a per-column basis. When R is given the Stata-equivalent times, correlations jump to 0.94-0.98.

---

### D-14: Survival Curve Column-2 Near-Constant Artifact [INFO]

**Suite:** 07, Tests 01-05, 09, 11-13

Column 2 of the survival curve (second unique failure time) consistently shows r = 0.84, just below the 0.85 threshold. This is because S(t) at very early times is ~1.0 for nearly all observations (range [0.977, 1.000]), leaving almost no variance for Pearson correlation to measure. All other columns (3-20) achieve r > 0.90.

---

### D-15: Multi-Arm W.hat Dimension Convention [LOW]

**Suite:** 09
**File:** `grf_multi_arm_causal_forest.ado`

R requires W.hat as an n x (K+1) matrix with column names matching treatment levels. Stata takes K binary propensity variables (one per non-reference treatment arm), and the reference arm probability is computed as 1 minus the sum.

---

### D-16: nomia Has No R Equivalent [INFO]

**Suite:** 01
**File:** `grf_regression_forest.ado`

R's grf 2.5.0 removed the `enable.missing.indicator` argument — MIA is always active when NAs are present. Stata's `nomia` flag can disable MIA, but this has no R equivalent to compare against.

---

### D-17: Default ci.group.size Differs [LOW]

**Suites:** 01, 02, 04, 09, 12

R defaults to `ci.group.size = 2`. Stata defaults to `cigroupsize(1)`. Users must explicitly set `cigroupsize(2)` in Stata to match R's variance estimation behavior.

---

### D-18: TMLE + Debiasing Weights Handling [INFO]

**Suite:** 03, Test 17

R silently ignores `debiasing.weights` when `method = "TMLE"` (falls through to standard TMLE). Stata correctly raises an error: "TMLE does not support debiasing weights." Stata's behavior is arguably better here.

---

### D-19: e(min_node) Naming Inconsistency [LOW]

**Suite:** 15
**Files:** `grf_regression_forest.ado`, `grf_tune.ado`

The stored result is named `e(min_node)` but the option is `minnodesize()` and R's argument is `min.node.size`. Users accessing `e(min_node_size)` will get missing. Should add an alias or rename.

---

### D-20: RATE AUTOC ~8-10% Underestimate [LOW]

**Suite:** 15, Tests 12, 14, 16

Stata's AUTOC estimates are consistently ~8-10% lower than R's (e.g., 1.034 vs 1.134). Both are within 3 standard errors (|diff|/SE = 0.46-0.76). The gap likely arises from different DR score computation paths: R uses `get_scores()` from the C++ engine with regularized AIPW; Stata recomputes DR scores from stored regression outputs.

QINI estimates show near-perfect agreement (|diff|/SE = 0.02).

---

### D-23: Boosted Forest Auto-Tune Step Count May Differ [LOW]

**Suite:** 10, Tests 01, 05, 07

In auto-tune mode (`booststeps(0)`), R may select 3 boosting steps while Stata selects 2, due to Monte Carlo variance in the pilot CV forests (`boosttreestune` default = 10 trees). Manual step specification eliminates this discrepancy. Predictions remain highly correlated (r > 0.995) regardless.

**Workaround:** Use `boosttreestune(50)` for more stable auto-tune agreement.

---

### D-24: Tuned Parameters Differ Between R and Stata [INFO]

**Suite:** 15, Tests 18, 20-21

R and Stata find different optimal hyperparameters via `tune.parameters` / `grf_tune` because they use different random search strategies and RNGs. For example: R finds mtry=10, min.node.size=5; Stata finds mtry=8, min.node.size=13. This is expected behavior for stochastic search with limited iterations.

---

### D-25: Calibration P-Values One-Sided vs Two-Sided [INFO]

**Suite:** 14, Tests 21-26
**File:** `grf_test_calibration.ado`

R's `test_calibration()` returns one-sided p-values (H0: coef <= 0). Stata returns two-sided p-values. When t > 0: Stata p ~ 2 * R p. T-statistics match within 1%.

---

### D-26: Calibration Coefficient Scaling Convention [INFO]

**Suite:** 14, Tests 21-23
**File:** `grf_test_calibration.ado`

R uses `(W-W_hat) * mean_tau_hat` as the first regressor (coefficient = ATE/mean_tau). Stata uses `(W-W_hat)` directly (coefficient = ATE). T-statistics are identical (scale-invariant). This is an intentional parametrization difference, not a bug.

---

## Category 4: Statistical Limitations (Expected Behavior)

### D-21: Null Treatment Effect Low Correlation [INFO]

**Suite:** 02, Test 19

When true tau(X) = 0 (no treatment effect), both R and Stata predict near-constant values around 0. With near-zero variance, Pearson correlation is dominated by noise (r = 0.81). Both implementations correctly identify the null effect (mean predictions < 0.01 in absolute value).

---

### D-22: Extreme Quantiles Low Correlation [INFO]

**Suite:** 05, Tests 10-11

At quantiles 0.01 and 0.99, predictions cluster near the tail of the distribution with limited variation. Pearson correlations are 0.78-0.82 (below 0.90 threshold). This is a statistical limitation at n=500, not an implementation divergence.

---

## Commands with 100% Pass Rate

Six commands achieved perfect fidelity across all tested configurations:

| Command | Tests | Min Correlation | Mean Correlation |
|---------|:-----:|:---------------:|:----------------:|
| regression_forest | 21 | 0.9529 | 0.9917 |
| probability_forest | 13 | 0.9913 | 0.9971 |
| boosted_regression_forest | 17 | 0.9008 | 0.9950 |
| ll_regression_forest | 18 | 0.9876 | 0.9962 |
| lm_forest | 17 | 0.8784* | 0.9892 |
| multi_regression_forest | 14 | 0.9768 | 0.9966 |

*lm_forest Test 3 (K=3, beta_3) has r=0.8784, above the 0.85 threshold for that suite.

---

## Priority Action Items

### Must Fix (3 items)

1. **D-1:** Add TMLE + clusters guard in `grf_ate.ado` — simple conditional error
2. **D-2:** Rewrite BLP debiasing.weights to use R's replacement-weight formula
3. **D-4:** Propagate cluster() to all nuisance regression forest calls in causal survival forest

### Should Fix (3 items)

4. **D-3:** Align ATE debiasing.weights SE formula with R's `sd(weighted)/sqrt(n)`
5. **D-5:** Implement full survival-forest-based IPCW pipeline for causal survival forest (major work)
6. **D-19:** Add `e(min_node_size)` alias for consistency

### Nice to Have (2 items)

7. **D-12:** Document that `nofastlogrank` matches R defaults in help file
8. **D-20:** Investigate AUTOC DR score computation alignment with R

---

## Files

| File | Description |
|------|-------------|
| `01_regression_forest.md` | 21 tests, 718 lines |
| `02_causal_forest_fit.md` | 19 tests, 665 lines |
| `03_ate.md` | 19 tests, 584 lines |
| `04_instrumental_forest.md` | 18 tests, 726 lines |
| `05_quantile_forest.md` | 17 tests, 540 lines |
| `06_probability_forest.md` | 13 tests, 400 lines |
| `07_survival_expected.md` | 17 tests, 214 lines |
| `08_causal_survival_forest.md` | 16 tests, 641 lines |
| `09_multi_arm_causal_forest.md` | 36 checks, 633 lines |
| `10_boosted_regression_forest.md` | 17 tests, 634 lines |
| `11_ll_regression_forest.md` | 18 tests, 509 lines |
| `12_lm_forest.md` | 17 tests, 480 lines |
| `13_multi_regression_forest.md` | 14 tests, 404 lines |
| `14_blp_calibration_scores.md` | 29 checks, 473 lines |
| `15_vi_rate_tune.md` | 24 tests, 361 lines |
| **Total** | **295 tests, ~7,982 lines** |
