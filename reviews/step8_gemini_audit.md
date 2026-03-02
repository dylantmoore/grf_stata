# Comprehensive Fidelity Audit of `grf_stata` against R's `grf` package

This report details a comprehensive review of the `grf_stata` implementation against the canonical R `grf` package. The audit evaluates feature coverage, parameter parity, default values, edge case handling, and test fidelity.

## 1. Summary Table

| R Function | Stata Equivalent | Test Files | Status |
|---|---|---|---|
| `regression_forest` | `grf_regression_forest` | `test_regression.do`, `test_options_regression.do`, `test_fidelity.do` | Excellent |
| `causal_forest` | `grf_causal_forest` | `test_causal.do`, `test_options_causal.do`, `test_fidelity.do` | Good (Missing nuisance inputs) |
| `quantile_forest` | `grf_quantile_forest` | `test_quantile.do`, `test_options_quantile.do`, `test_fidelity.do` | Excellent |
| `instrumental_forest` | `grf_instrumental_forest` | `test_instrumental.do`, `test_options_instrumental.do`, `test_fidelity.do` | Good (Missing nuisance outputs/inputs) |
| `probability_forest` | `grf_probability_forest` | `test_probability.do`, `test_options_probability.do`, `test_fidelity.do` | Excellent |
| `survival_forest` | `grf_survival_forest` | `test_survival.do`, `test_options_survival.do`, `test_fidelity.do` | Excellent |
| `causal_survival_forest` | `grf_causal_survival_forest` | `test_causal_survival.do`, `test_options_causal_survival.do`, `test_fidelity.do` | Good (Missing nuisance outputs/inputs) |
| `multi_arm_causal_forest` | `grf_multi_arm_causal_forest` | `test_multi_arm.do`, `test_options_multi_arm.do`, `test_fidelity.do` | Good (Missing nuisance outputs/inputs) |
| `multi_regression_forest` | `grf_multi_regression_forest` | `test_multi_regression.do`, `test_options_multi_regression.do`, `test_fidelity.do` | Excellent |
| `ll_regression_forest` | `grf_ll_regression_forest` | `test_ll_regression.do`, `test_options_ll_regression.do`, `test_fidelity.do` | Good (Limited llsplit support) |
| `boosted_regression_forest` | `grf_boosted_regression_forest` | `test_boosted_regression.do`, `test_options_boosted.do`, `test_fidelity.do` | Excellent |
| `lm_forest` | `grf_lm_forest` | `test_lm_forest.do`, `test_options_lm_forest.do`, `test_fidelity.do` | Excellent |
| `predict` | `grf_predict` | `test_predict.do`, `test_fidelity.do` | Excellent |
| `average_treatment_effect` | `grf_ate` | `test_post_estimation.do`, `test_options_post_estimation.do` | Excellent |
| `best_linear_projection` | `grf_best_linear_projection` | `test_post_estimation.do`, `test_options_post_estimation.do` | Excellent |
| `test_calibration` | `grf_test_calibration` | `test_post_estimation.do`, `test_options_post_estimation.do` | Excellent |
| `variable_importance` | `grf_variable_importance` | `test_post_estimation.do`, `test_options_post_estimation.do` | Excellent |
| `rank_average_treatment_effect` | `grf_rate` | `test_rate.do`, `test_options_post_estimation.do` | Excellent |
| `get_scores` | `grf_get_scores` | `test_get_scores.do`, `test_options_post_estimation.do` | Excellent |
| `average_partial_effect` | `grf_average_partial_effect` | `test_average_partial_effect.do`, `test_options_post_estimation.do` | Excellent |

## 2. Parameter-by-Parameter Comparison & Divergence

The Stata implementation successfully wraps the C++ layer with high parity, but deviates from R in a few deliberate and unintentional ways:

**Common Parameter Deviations:**
* `seed`: R uses a randomly drawn integer by default (`runif(1, 0, .Machine$integer.max)`). Stata uses a fixed default of `42`. This ensures out-of-the-box reproducibility in Stata (which matches standard Stata expectations), but differs from R's default stochasticity.
* `ci.group.size`: R defaults to `2` to enable variance estimates automatically. Stata defaults to `1` (which disables variance estimation), requiring users to explicitly pass `estimatevariance` or `cigroupsize(2)` to compute standard errors.
* `equalize.cluster.weights`: Available in R (default `FALSE`), but **completely missing** from the Stata implementation syntax across all models.

**Forest-Specific Parameter Deviations:**
* **`causal_forest`**: 
  * R allows users to input precomputed `Y.hat` and `W.hat` variables. Stata only provides `yhatgenerate(varname)` and `whatgenerate(varname)` to output them, but does not allow passing custom nuisance parameters to override the internal regression step.
  * `orthog.boosting` is missing from Stata.
* **`instrumental_forest`**: 
  * Similar to `causal_forest`, `Y.hat`, `W.hat`, and `Z.hat` are generated internally using tempvars, but Stata exposes no options to output them (no `*generate()` options) or input them.
* **`causal_survival_forest`**:
  * Missing options to input or explicitely generate `Y.hat`, `W.hat`, `S.hat`, `C.hat`.
* **`multi_arm_causal_forest`**:
  * Missing options to input/output nuisance parameters. It silently drops variables like `_grf_mac_yhat` globally into the dataset without a prefix or option to name them cleanly.
  * `orthog.boosting` is missing from Stata.
* **`ll_regression_forest`**:
  * R's `ll.split.variables` allows selecting a specific subset of variables to use for local linear splits. Stata's `llsplit` option is only a boolean toggle that hardcodes splitting on *all* variables.
* **`quantile_forest`**:
  * Missing the `regression.splitting` option present in R.

## 3. Issues Found & Severity Ratings

| Issue | Severity | Description |
|---|---|---|
| Missing Nuisance Inputs | **Major** | Users cannot supply precomputed `Y.hat`, `W.hat`, etc., for `causal_forest`, `instrumental_forest`, or `multi_arm_causal_forest`. This breaks a common advanced workflow where users fit nuisance parameters using highly tuned models (e.g., neural networks or lasso outside of `grf`) and pass them in. |
| Missing `equalize.cluster.weights` | **Major** | This parameter is completely absent from the `.ado` definitions. |
| Sloppy Global Output in `multi_arm_causal_forest` | **Minor** | It drops variables named `_grf_mac_yhat` and `_grf_mac_what*` into the dataset directly rather than using standard `generate()` conventions as `causal_forest` does. |
| Missing `ll.split.variables` selection | **Minor** | The boolean `llsplit` limits the flexibility of the local linear forest by forcing it on all variables. |
| Missing `regression.splitting` in Quantile Forest | **Minor** | Option not exposed in the Stata wrapper. |
| Differing defaults (`seed` and `ci.group.size`) | **Info** | Stata forces `seed(42)` and `cigroupsize(1)`. This is more of a design choice appropriate for Stata, but should be distinctly highlighted in the documentation for users migrating from R. |

## 4. Overall Fidelity Score: 8.5 / 10

**Justification:** 
The `grf_stata` library is an extremely robust and reliable wrapper. The test suite is meticulously constructed with over thousands of lines in `tests/test_fidelity.do` directly correlating outputs against `.csv` dumps from the R implementation, effectively ensuring the C++ binding achieves the same underlying statistical computations. Edge cases around missing data (`MIA`), weights, out-of-sample prediction (`test_predict.do`), and restricted estimations (`if`/`in`) are validated exceptionally well.

However, it loses 1.5 points strictly on parameter flexibility. R's `grf` strongly encourages modular workflows (passing in custom `Y.hat` and `W.hat`). By hardcoding these to be internally generated and failing to expose them as inputs, the Stata implementation limits the canonical "double machine learning" approach that researchers often use `grf` for. 

## 5. Recommendations

1. **Implement Nuisance Parameter Inputs (`yhat()`, `what()`, `zhat()`):** Update the syntax definitions in `causal_forest.ado`, `instrumental_forest.ado`, `causal_survival_forest.ado`, and `multi_arm_causal_forest.ado` to accept existing variables for the nuisance parameters, bypassing the internal nuisance forest estimation when provided.
2. **Standardize Nuisance Outputs:** Align `multi_arm_causal_forest` and `instrumental_forest` to behave like `causal_forest` by introducing `*generate()` options instead of polluting the workspace with hardcoded `_grf_*` variable names.
3. **Add Missing Parameters:** Plumb `equalizeclusterweights`, `regression_splitting`, and specific `llsplitvars(varlist)` arguments through the Stata `.ado` parsers into the C++ plugin layer.
4. **Documentation Highlight:** Explicitly note in `grf.sthlp` that `grf_stata` defaults to a fixed seed (42) and disables variance estimation by default (cigroupsize=1), unlike R.
