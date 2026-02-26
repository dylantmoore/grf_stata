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
| `grf_best_linear_projection` | Project CATEs onto covariates (BLP) |
| `grf_test_calibration` | Calibration test for forest predictions |
| `grf_variable_importance` | Variable importance via weighted split frequencies |
| `grf_rate` | Rank-weighted average treatment effect (AUTOC/QINI) |
| `grf_predict` | Predict on new (out-of-sample) data |
| `grf_tune` | Cross-validation tuning of forest parameters |

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
- **7 post-estimation commands**: ATE, BLP, calibration test, variable importance, RATE, predict, tune
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

## Stored Results

All forest commands store results in `e()`:

- `e(N)` — number of observations
- `e(n_trees)` — number of trees
- `e(ate)` / `e(ate_se)` — average treatment effect and SE (causal forests)
- `e(cmd)` — command name
- `e(forest_type)` — forest type identifier
- `e(depvar)`, `e(treatvar)`, `e(indepvars)` — variable names

See `help` for each command for complete stored results.

## Requirements

- Stata 14.0 or later
- No external dependencies (grf C++ library is statically linked)

## Platform Support

| Platform | Binary |
|----------|--------|
| macOS ARM64 (Apple Silicon) | `grf_plugin.darwin-arm64.plugin` |
| macOS x86_64 (Intel) | `grf_plugin.darwin-x86_64.plugin` |
| Linux x86_64 | `grf_plugin.linux-x86_64.plugin` |
| Windows x86_64 | `grf_plugin.windows-x86_64.plugin` |

## References

Athey, S., J. Tibshirani, and S. Wager. 2019. "Generalized Random Forests." *Annals of Statistics* 47(2): 1148-1178.

Wager, S. and S. Athey. 2018. "Estimation and Inference of Heterogeneous Treatment Effects using Random Forests." *Journal of the American Statistical Association* 113(523): 1228-1242.

## License

GPL-3.0. Wraps the grf C++ library ([grf-labs/grf](https://github.com/grf-labs/grf)) by Susan Athey, Julie Tibshirani, Stefan Wager, and contributors.
