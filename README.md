# grf_stata: Generalized Random Forests for Stata

Stata implementation of the [grf](https://github.com/grf-labs/grf) R package (v2.5.0) by Athey, Tibshirani, and Wager. Wraps the same C++ backend used by the R package, compiled as a native Stata plugin for full performance.

## Installation

```stata
* macOS:
net install grf_stata_mac, from("https://raw.githubusercontent.com/dylantmoore/grf_stata/main") replace

* Linux:
net install grf_stata_linux, from("https://raw.githubusercontent.com/dylantmoore/grf_stata/main") replace

* Windows:
net install grf_stata_win, from("https://raw.githubusercontent.com/dylantmoore/grf_stata/main") replace

* All platforms (larger download):
net install grf_stata, from("https://raw.githubusercontent.com/dylantmoore/grf_stata/main") replace
```

## Forest Types

| Command | Description |
|---------|-------------|
| `grf_regression_forest` | Regression forest for conditional mean estimation |
| `grf_causal_forest` | Causal forest for heterogeneous treatment effects (CATE) |
| `grf_quantile_forest` | Quantile forest for conditional quantile estimation |
| `grf_instrumental_forest` | Instrumental forest for LATE with instruments |
| `grf_probability_forest` | Probability forest for class probability estimation |
| `grf_survival_forest` | Survival forest for conditional survival functions |
| `grf_causal_survival_forest` | Causal survival forest for treatment effects on survival |
| `grf_multi_arm_causal_forest` | Multi-arm causal forest for multi-valued treatments |
| `grf_multi_regression_forest` | Multi-outcome regression forest |
| `grf_ll_regression_forest` | Local linear regression forest |
| `grf_boosted_regression_forest` | Boosted regression forest (iterative residual fitting) |
| `grf_lm_forest` | Linear model forest (conditional linear coefficients) |

## Post-Estimation Commands

| Command | Description |
|---------|-------------|
| `grf_ate` | Average treatment effect (AIPW doubly-robust estimator) |
| `grf_average_partial_effect` | Average partial effect for continuous treatments |
| `grf_best_linear_projection` | Project CATEs onto covariates (BLP) |
| `grf_test_calibration` | Calibration test for forest predictions |
| `grf_variable_importance` | Variable importance via weighted split frequencies |
| `grf_get_scores` | Extract doubly-robust/IPCW score variables |
| `grf_rate` | Rank-weighted average treatment effect (AUTOC/QINI) |
| `grf_predict` | Predict on new (out-of-sample) data |
| `grf_tune` | Cross-validation tuning of forest parameters |

Utility/introspection commands:

- `grf_forest_summary`
- `grf_tree_summary`
- `grf_get_tree`
- `grf_get_leaf_node`
- `grf_get_forest_weights`
- `grf_merge_forests`
- `grf_split_frequencies`
- `grf_plot_tree`

## Quick Start

```stata
* Simulate data with heterogeneous treatment effects
clear
set obs 2000
set seed 12345
gen x1 = rnormal()
gen x2 = rnormal()
gen w = rbinomial(1, normal(x1))
gen y = 2*x1 + w*(1 + x2) + rnormal()

* Estimate CATEs
grf_causal_forest y w x1 x2, gen(tau_hat)

* Average treatment effect
grf_ate

* Best linear projection
grf_best_linear_projection x1 x2

* Calibration test
grf_test_calibration

* Variable importance
grf_variable_importance x1 x2

* Rank-weighted ATE
grf_rate tau_hat
```

## Features

- **Full grf C++ backend**: Same algorithms as the R package, not a reimplementation
- **12 forest types**: Regression, causal, quantile, instrumental, probability, survival, causal survival, multi-arm causal, multi-regression, local linear, boosted, LM forest
- **Post-estimation suite**: ATE, APE, BLP, calibration test, variable importance, DR/IPCW scores, RATE, predict, tune
- **Utility/introspection suite**: forest/tree summaries, tree/leaf extraction helpers, proxy weights, split-frequency proxies, merge helpers, and plot wrapper
- **Honest estimation**: Split-selection and leaf-estimation on disjoint subsamples (default)
- **Variance estimation**: Out-of-bag variance estimates for CATEs
- **Prediction on new data**: Append test observations and predict with `grf_predict`
- **Cross-validation tuning**: Automatic parameter tuning via `grf_tune`
- **Cross-platform**: macOS (ARM64, x86_64), Linux (x86_64), Windows (x86_64)
- **Standard Stata interface**: `if`/`in` restrictions, `replace` option, `e()` stored results

## Common Options

All forest commands share these options:

| Option | Default | Description |
|--------|---------|-------------|
| `generate(newvar)` | required | Name for prediction variable |
| `ntrees(#)` | 2000 | Number of trees |
| `seed(#)` | 42 | Random seed |
| `mtry(#)` | 0 (sqrt(p)) | Variables per split |
| `minnodesize(#)` | 5 | Minimum leaf size |
| `samplefrac(#)` | 0.5 | Subsample fraction |
| `honesty` / `nohonesty` | honesty | Honest estimation |
| `honestyfrac(#)` | 0.5 | Honest hold-out fraction |
| `alpha(#)` | 0.05 | Split imbalance bound |
| `numthreads(#)` | 0 (all) | Number of threads |
| `estimatevariance` | off | Compute variance estimates |
| `replace` | off | Overwrite existing variables |
| `cluster(varname)` | none | Cluster variable for cluster-robust forests |
| `weights(varname)` | none | Sample weights variable (maps to R `sample.weights`) |
| `nomia` | MIA on | Disable Missing Indicator Action (casewise deletion) |

## Stored Results

All forest commands store results in `e()`:

- `e(N)` — number of observations
- `e(n_trees)` — number of trees
- `e(model_id)` — session-local model identifier (increments each fit)
- `e(ate)` / `e(ate_se)` — average treatment effect and SE (causal forests)
- `e(cmd)` — command name
- `e(forest_type)` — forest type identifier
- `e(depvar)`, `e(treatvar)`, `e(indepvars)` — variable names

See `help` for each command for complete stored results.

## Known Limitations vs R's grf

- **User-supplied nuisance estimates**: supported for several forests (`yhatinput()`, `whatinput()`, `zhatinput()`, and related generators), but not every R nuisance hook is exposed for every forest family.
- **Some forest-inspection utilities**: full tree/forest object persistence is constrained by the plugin architecture.
- **APE deprecation**: R has deprecated `average_partial_effect()`. Stata retains `grf_average_partial_effect` for backward compatibility.
- **Causal survival scores**: `grf_get_scores` uses `tau + psi / V.hat` with stored nuisance moments (`psi = numer - denom * tau`), matching current upstream score construction at the wrapper level.
- **Package-level options API**: R's `grf_options()` is not mirrored; Stata uses per-command options by design.
- **OOB prediction toggle**: R's `compute.oob.predictions` option is not exposed; Stata fit commands always materialize prediction outputs by command design.
- **Upstream-constrained parity**: options such as `orthog.boosting` and enum-style `honesty.prune.method` are not exposed under the current vendored core API.

### Introspection Coverage Notes

- Partial introspection utilities are implemented and usable in Stata (`grf_forest_summary`, `grf_tree_summary`, `grf_get_tree`, `grf_get_leaf_node`, `grf_get_forest_weights`, `grf_split_frequencies`, `grf_plot_tree`).
- `grf_forest_summary, all` now exposes all stored `e()` scalar/macro names for model-level inspection.
- `grf_get_forest_weights` now supports combined prediction-space and feature-space proxy weights via `xvars()` and `predweight()`.
- Remaining non-equivalent internals (exact node structures, exact forest weights, exact terminal-node IDs) are documented in:
  `reviews/introspection_discrepancies.md`.

## Parity Baseline Artifacts (Gaps 4/5/6)

- Baseline pin: `reviews/parity_baseline.md`
- Machine-extracted upstream API manifest: `reviews/r_api_manifest.json`
- Generator script: `tools/extract_r_api_manifest.R`
- Scope check script: `tools/check_parity_scope.R`

Run:

```bash
Rscript tools/extract_r_api_manifest.R --tag v2.5.0 --out reviews/r_api_manifest.json
Rscript tools/check_parity_scope.R reviews/r_api_manifest.json
```

## Migration Note (Causal Survival)

- `grf_causal_survival_forest` now enforces explicit nuisance modes.
- Partial nuisance override is no longer accepted. If using nuisance inputs, provide all of:
  `whatinput()`, `yhatinput()`, `shatinput()`, and `chatinput()`.
- `numer()/denom()` is now an explicit `moment_input` override and is mutually exclusive with nuisance-input mode.

## Testing

Run the standard regression suite from the project root:

```stata
do tests/run_all.do
```

## Requirements

- Stata 14.0 or later
- No external dependencies (grf C++ library is statically linked)

## Platform Support

| Platform | Binary |
|----------|--------|
| macOS (Apple Silicon) | `grf_plugin_macosx.plugin` |
| Linux x86_64 | `grf_plugin_unix.plugin` |
| Windows x86_64 | `grf_plugin_windows.plugin` |

## References

Athey, S., J. Tibshirani, and S. Wager. 2019. "Generalized Random Forests." *Annals of Statistics* 47(2): 1148-1178.

Wager, S. and S. Athey. 2018. "Estimation and Inference of Heterogeneous Treatment Effects using Random Forests." *Journal of the American Statistical Association* 113(523): 1228-1242.

## License

GPL-3.0. Wraps the grf C++ library ([grf-labs/grf](https://github.com/grf-labs/grf)) by Susan Athey, Julie Tibshirani, Stefan Wager, and contributors.
