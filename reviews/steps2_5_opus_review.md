# Review of grf_stata Steps 2-5 Implementation

**Reviewer:** Claude Opus 4.6
**Date:** 2026-02-27
**Repository:** https://github.com/dylantmoore/grf_stata (commit f41d6b4)
**Verdict:** NOT LGTM -- Steps 2, 3, and 4 are entirely unimplemented. Step 5 is partially implemented with correctness issues.

---

## Executive Summary

The gap-closing plan described Steps 2-5 as implemented. Upon thorough review of the
codebase at commit f41d6b4, **Steps 2, 3, and 4 have not been implemented at all**, and
Step 5 is only partially implemented with several deviations from R's grf. Additionally,
the existing DR score computation across all post-estimation commands uses an incorrect
formula that diverges from R's grf.

---

## Step 2: MIA Missing Data Passthrough -- NOT IMPLEMENTED

**Status: No changes found.**

### What was claimed:
- `grf_plugin.cpp` modified to pass NaN for missing covariate values
- `allow_missing_x` flag at argv[20]
- All 12 forest .ado files: `noMIA` syntax option
- `tests/test_mia.do`

### What exists:
- **grf_plugin.cpp**: The data reading loop (lines 274-296) still performs full casewise
  deletion across ALL columns (including X). No NaN passthrough. No `allow_missing_x`
  parameter. argv[20+] remain forest-specific parameters (e.g., `stabilize_splits`
  for causal, `quantiles` for quantile).
- **Forest .ado files**: No `noMIA` syntax option in any of the 12 forest .ado files.
  Confirmed by grep for `noMIA`, `allow_missing`, `MIA` -- zero matches.
- **tests/test_mia.do**: File does not exist.

### Impact:
The entire argv layout described in the spec (argv[20]=allow_missing_x,
argv[21]=cluster_col_idx, argv[22]=weight_col_idx, forest-specific at argv[23+])
does not exist. The current layout has forest-specific params starting at argv[20].

---

## Step 3: Cluster Support -- NOT IMPLEMENTED

**Status: No changes found.**

### What was claimed:
- `grf_plugin.cpp`: `cluster_col_idx` at argv[21], cluster vector passed to ForestOptions
- All 12 forest .ado files: `CLuster(varname numeric)` syntax
- Post-estimation commands updated for cluster-robust SEs
- `tests/test_clusters.do`

### What exists:
- **grf_plugin.cpp** (line 366-367): `clusters` vector is declared as empty and
  `samples_per_cluster = 0`. These are never populated from any input parameter.
  ForestOptions always receives empty clusters.
- **Forest .ado files**: No `CLuster` syntax in any .ado file. No `e(cluster_var)` stored.
- **Post-estimation commands**: No cluster-robust SE paths in grf_ate.ado,
  grf_best_linear_projection.ado, or grf_test_calibration.ado.
- **tests/test_clusters.do**: File does not exist.

---

## Step 4: Sample Weights Support -- NOT IMPLEMENTED

**Status: No changes found.**

### What was claimed:
- `grf_plugin.cpp`: `weight_col_idx` at argv[22], weight vector passed as sample_weights
- All 12 forest .ado files: `WEIghts(varname numeric)` syntax
- `tests/test_weights.do`

### What exists:
- **grf_plugin.cpp**: No `weight_col_idx` parameter. No sample weights passed to any
  trainer. The grf C++ library accepts sample weights via `ForestTrainer::train(data, options)`
  with an optional `sample_weights` vector, but this is never populated.
- **Forest .ado files**: No `WEIghts` syntax in any .ado file. No `e(weight_var)` stored.
- **tests/test_weights.do**: File does not exist.

---

## Step 5: Extended Post-Estimation -- PARTIALLY IMPLEMENTED with issues

### 5a. grf_get_scores.ado

**Status: Implemented for causal and instrumental forests only. Multi-arm causal and
causal survival are NOT implemented.**

#### What was claimed:
- Extended to support `multi_causal` (per-arm DR scores)
- Extended to support `causal_survival` (simplified DR scores)

#### What exists (grf_get_scores.ado, lines 11-16):
The forest type check only allows "causal" and "instrumental":
```stata
if "`forest_type'" != "causal" & "`forest_type'" != "instrumental" {
    display as error "grf_get_scores requires prior estimation by" ///
        " grf_causal_forest or grf_instrumental_forest"
    exit 301
}
```
No multi_arm_causal or causal_survival code paths exist.

#### DR Score Formula Issue (affects ALL post-estimation commands):

**grf_get_scores.ado** (lines 93-110), **grf_ate.ado** (lines 61-74), and
**grf_best_linear_projection.ado** (lines 90-103) all use the SAME formula:

```stata
gen double `w_resid' = `treatvar' - `whatvar' if `touse'
summarize `w_resid' if `touse'
local w_resid_var = r(Var)
...
gen double `dr_score' = `tauvar' + (`w_resid' / `w_resid_var') * `y_resid'
```

This computes: `DR_i = tau_hat_i + (W_i - W_hat_i) / Var_sample(W - W_hat) * Y_resid_i`

**R's grf** (get_scores.R, line 63) for binary treatment computes:
```r
debiasing.weights <- (W.orig - W.hat) / (W.hat * (1 - W.hat))
```

This computes: `DR_i = tau_hat_i + (W_i - W_hat_i) / (W_hat_i * (1 - W_hat_i)) * Y_resid_i`

**The difference**: Stata uses a **scalar** sample variance of W-W.hat as the denominator
for all observations, while R uses the **observation-specific** conditional variance
`W_hat_i * (1 - W_hat_i)`. These are mathematically different:
- R: heterogeneous debiasing weights that vary by propensity score
- Stata: homogeneous debiasing weights (same denominator for every observation)

For binary treatment with constant propensity (randomized experiment), these converge.
But with heterogeneous propensities, the Stata formula is incorrect. R's approach produces
the standard AIPW/Horvitz-Thompson-style estimator; the Stata approach loses the
inverse-propensity weighting structure entirely.

**Severity: HIGH** -- This affects grf_get_scores, grf_ate, and grf_best_linear_projection.

For continuous treatment, R trains an auxiliary variance forest to estimate E[(W-W_hat)^2|X],
giving heterogeneous V_hat. The Stata code's use of sample variance is a cruder approximation
but more defensible in the continuous case as a simplification. However, for binary treatment
it is clearly wrong.

### 5b. grf_best_linear_projection.ado

**Status: Partially implemented. "treated" and "control" target.sample NOT added.**

#### What was claimed:
- Extended target.sample to accept "treated" (weight=W_hat) and "control" (weight=1-W_hat)

#### What exists (grf_best_linear_projection.ado, lines 41-48):
```stata
if "`targetsample'" == "" {
    local targetsample "all"
}
if "`targetsample'" != "all" & "`targetsample'" != "overlap" {
    display as error "target.sample must be one of: all, overlap"
    exit 198
}
```
Only "all" and "overlap" are supported.

#### R grf reference check:
R's grf `best_linear_projection()` (forest_summary.R, line 171) also ONLY supports
`target.sample = c("all", "overlap")`. It does NOT support "treated" or "control" for BLP.
The spec's claim that BLP should support "treated" and "control" is itself inconsistent
with R's grf -- this was a spec error, not an implementation gap.

### 5c. grf_ate.ado target.sample

**Status: "treated", "control", "overlap", and "all" ARE implemented.** This appears to have
been completed in commit f41d6b4.

#### Correctness issue with target.sample weighting:

grf_ate.ado (lines 79-90) uses:
```stata
if "`targetsample'" == "treated" {
    gen double `tsweight' = `treatvar' if `touse'
}
else if "`targetsample'" == "control" {
    gen double `tsweight' = (1 - `treatvar') if `touse'
}
```

R's grf ATE for "treated" and "control" uses a more complex TMLE/AIPW procedure
(average_treatment_effect.R, lines 298-338) that involves:
1. Computing a raw weighted average of tau_hat on the treated/control subset
2. Computing correction terms using `gamma` weights:
   - For treated: `gamma_control = W_hat / (1 - W_hat)`, `gamma_treated = 1`
   - For control: `gamma_control = 1`, `gamma_treated = (1 - W_hat) / W_hat`
3. A doubly-robust correction using these gamma weights

The Stata implementation weights DR scores by `W_i` (treated) or `1 - W_i` (control),
which is a simpler "Horvitz-Thompson-like" approach. This is a valid simplification
that gives a consistent estimator under correct specification of either the outcome or
propensity model, but it differs from R's more efficient AIPW/TMLE approach for ATT/ATC.

**Severity: MEDIUM** -- The estimator is consistent but less efficient than R's approach.

### 5d. grf_test_calibration.ado

**Status: Implemented, no cluster VCE path.** The current implementation uses `vce(hc3)`
(line 91). No cluster-robust option exists since clusters are not implemented (Step 3).

### 5e. Tests for Step 5

**test_options_post_estimation.do**: Tests go up to Test 78 (lines 1513-1527). Tests 79-86
for BLP treated/control and get_scores multi-arm/causal-survival do NOT exist.

**test_get_scores.do**: Tests 1-7 exist and cover basic causal and instrumental forest
DR score extraction, replace option, error handling, and mean-approximates-ATE property.
No tests for multi-arm or causal-survival.

---

## Summary of Issues

### Critical (must fix):

| # | Issue | File | Lines |
|---|-------|------|-------|
| 1 | **Steps 2-4 entirely unimplemented** | All files | N/A |
| 2 | **DR score formula uses scalar sample variance instead of observation-specific W_hat*(1-W_hat)** | grf_get_scores.ado | 95-110 |
| 3 | Same DR formula bug | grf_ate.ado | 61-74 |
| 4 | Same DR formula bug | grf_best_linear_projection.ado | 90-103 |

### High (should fix):

| # | Issue | File | Lines |
|---|-------|------|-------|
| 5 | grf_get_scores does not support multi_arm_causal or causal_survival | grf_get_scores.ado | 11-16 |
| 6 | grf_ate ATT/ATC uses simple weighting instead of R's AIPW/TMLE approach | grf_ate.ado | 79-90 |

### Spec issues (not code bugs):

| # | Issue | Notes |
|---|-------|-------|
| 7 | Spec claims BLP should support "treated"/"control" target.sample | R's grf only supports "all" and "overlap" for BLP |
| 8 | Spec claims tests 79-86 exist | test_options_post_estimation.do ends at test 78 |

### Missing files:

- `tests/test_mia.do` -- does not exist
- `tests/test_clusters.do` -- does not exist
- `tests/test_weights.do` -- does not exist

---

## Recommendations

1. **Implement Steps 2-4 before reviewing them.** The current codebase has no MIA, cluster,
   or sample weight support. This requires:
   - Shifting the argv layout: common params at [0-19], new params at [20-22],
     forest-specific at [23+]
   - Updating all 12 forest .ado files and their nuisance regression calls
   - Updating grf_plugin.cpp to read cluster/weight columns and pass them through

2. **Fix the DR score formula.** For binary treatment, replace:
   ```stata
   summarize `w_resid' if `touse'
   local w_resid_var = r(Var)
   gen double `dr_score' = `tauvar' + (`w_resid' / `w_resid_var') * `y_resid'
   ```
   with:
   ```stata
   gen double `v_hat' = `whatvar' * (1 - `whatvar') if `touse'
   gen double `dr_score' = `tauvar' + (`w_resid' / `v_hat') * `y_resid' if `touse'
   ```
   This should apply to grf_get_scores.ado, grf_ate.ado, and grf_best_linear_projection.ado.

   For continuous treatment, the ideal approach is to train an auxiliary variance forest
   as R does, but the sample variance approximation is a reasonable simplification if
   documented. A runtime check for binary vs continuous W should branch between the two
   formulas.

3. **Drop BLP "treated"/"control" from the spec.** R's grf does not support these. Keep
   BLP at "all" and "overlap" only, matching R.

4. **Implement multi-arm and causal-survival DR scores** in grf_get_scores.ado to match
   R's get_scores.multi_arm_causal_forest and get_scores.causal_survival_forest methods.

5. **Write the missing test files** (test_mia.do, test_clusters.do, test_weights.do) and
   add tests 79+ for the extended post-estimation features.
