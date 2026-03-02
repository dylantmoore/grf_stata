# Comprehensive Fidelity Audit: grf_stata

**Auditor:** Claude Opus 4.6
**Date:** 2026-02-27
**Scope:** Full parameter-level comparison of R's grf package vs. grf_stata

---

## 1. Function Mapping: R grf --> Stata Commands --> Test Files

### 1.1 Forest Training Functions

| R Function | Stata Command | Test File(s) | Status |
|---|---|---|---|
| `causal_forest()` | `grf_causal_forest` | test_causal.do, test_options_causal.do | **Implemented** |
| `regression_forest()` | `grf_regression_forest` | test_regression.do, test_options_regression.do | **Implemented** |
| `quantile_forest()` | `grf_quantile_forest` | test_quantile.do, test_options_quantile.do | **Implemented** |
| `instrumental_forest()` | `grf_instrumental_forest` | test_instrumental.do, test_options_instrumental.do | **Implemented** |
| `probability_forest()` | `grf_probability_forest` | test_probability.do, test_options_probability.do | **Implemented** |
| `survival_forest()` | `grf_survival_forest` | test_survival.do, test_options_survival.do | **Implemented** |
| `causal_survival_forest()` | `grf_causal_survival_forest` | test_causal_survival.do, test_options_causal_survival.do | **Implemented** (partial nuisance) |
| `multi_arm_causal_forest()` | `grf_multi_arm_causal_forest` | test_multi_arm.do, test_options_multi_arm.do | **Implemented** |
| `multi_regression_forest()` | `grf_multi_regression_forest` | test_multi_regression.do, test_options_multi_regression.do | **Implemented** |
| `ll_regression_forest()` | `grf_ll_regression_forest` | test_ll_regression.do, test_options_ll_regression.do | **Implemented** |
| `boosted_regression_forest()` | `grf_boosted_regression_forest` | test_boosted_regression.do, test_options_boosted.do | **Implemented** |
| `lm_forest()` | `grf_lm_forest` | test_lm_forest.do, test_options_lm_forest.do | **Implemented** |

### 1.2 Post-Estimation / Analysis Functions

| R Function | Stata Command | Test File(s) | Status |
|---|---|---|---|
| `average_treatment_effect()` | `grf_ate` | test_post_estimation.do, test_options_post_estimation.do | **Implemented** |
| `best_linear_projection()` | `grf_best_linear_projection` | test_post_estimation.do, test_options_post_estimation.do | **Implemented** |
| `test_calibration()` | `grf_test_calibration` | test_post_estimation.do, test_options_post_estimation.do | **Implemented** |
| `variable_importance()` | `grf_variable_importance` | test_post_estimation.do, test_options_post_estimation.do | **Implemented** |
| `rank_average_treatment_effect()` | `grf_rate` | test_rate.do, test_options_post_estimation.do | **Implemented** |
| `get_scores()` | `grf_get_scores` | test_get_scores.do | **Implemented** |
| `predict.*()` (all forest types) | `grf_predict` | test_predict.do | **Implemented** |

### 1.3 Utility Functions

| R Function | Stata Command | Status |
|---|---|---|
| `get_forest_weights()` | -- | **Missing** (documented R-only: returns kernel weight matrix; memory-intensive for large N) |
| `get_leaf_node()` | -- | **Missing** (tree introspection not exposed) |
| `get_tree()` | -- | **Missing** (tree introspection not exposed) |
| `merge_forests()` | -- | **Missing** (forest object persistence not supported in plugin architecture) |
| `split_frequencies()` | Partially via `grf_variable_importance` | **Partial** (variable_importance uses split_frequencies internally) |
| `plot.grf_tree()` | -- | **Missing** (not applicable to Stata workflow) |
| `plot.rank_average_treatment_effect()` | -- | **Missing** (Stata graphing could be added but is not) |
| `print.*()` (all types) | Display output in each .ado | **Implemented** (via display statements) |
| `grf_options()` | -- | **Missing** (no global option registry) |

### 1.4 Additional Stata-only Commands

| Stata Command | R Equivalent | Status |
|---|---|---|
| `grf_tune` | `tune.parameters` arg in forest functions | **Stata-specific** (standalone tuning command) |
| `grf_average_partial_effect` | Deprecated `average_partial_effect()` | **Implemented** (backward compat) |
| `grf_expected_survival` | Manual computation of E[T|X] | **Stata-specific** (convenience helper) |
| `grf_generate_causal_data` | R example DGPs | **Stata-specific** (data generation) |
| `grf_generate_causal_survival_data` | R example DGPs | **Stata-specific** (data generation) |

---

## 2. Parameter-by-Parameter Comparison

### 2.1 regression_forest()

| R Parameter | R Default | Stata Option | Stata Default | Match? |
|---|---|---|---|---|
| `X` | required | varlist (indep) | required | Yes |
| `Y` | required | varlist (dep) | required | Yes |
| `num.trees` | 2000 | `ntrees()` | 2000 | Yes |
| `sample.weights` | NULL | `weights()` | none | Yes |
| `clusters` | NULL | `cluster()` | none | Yes |
| `equalize.cluster.weights` | FALSE | -- | -- | **Missing** |
| `sample.fraction` | 0.5 | `samplefrac()` | 0.5 | Yes |
| `mtry` | min(ceil(sqrt(p)+20),p) | `mtry()` | 0 (auto in C++) | Yes (0=auto) |
| `min.node.size` | 5 | `minnodesize()` | 5 | Yes |
| `honesty` | TRUE | `nohonesty` toggle | TRUE | Yes |
| `honesty.fraction` | 0.5 | `honestyfrac()` | 0.5 | Yes |
| `honesty.prune.leaves` | TRUE | `nohonestyprune` toggle | TRUE | Yes |
| `alpha` | 0.05 | `alpha()` | 0.05 | Yes |
| `imbalance.penalty` | 0 | `imbalancepenalty()` | 0.0 | Yes |
| `ci.group.size` | 2 | `cigroupsize()` | 1 | **Mismatch** (R=2, Stata=1; see note) |
| `tune.parameters` | "none" | `grf_tune` (separate) | N/A | **Different approach** |
| `tune.num.trees` | 50 | `grf_tune tunetrees()` | 200 | Separate command |
| `tune.num.reps` | 100 | `grf_tune numreps()` | 50 | Separate command |
| `tune.num.draws` | 1000 | -- | -- | **Missing** (not exposed) |
| `compute.oob.predictions` | TRUE | Always TRUE | TRUE | Yes |
| `num.threads` | NULL (all) | `numthreads()` | 0 (all) | Yes |
| `seed` | random | `seed()` | 42 | **Mismatch** (R=random, Stata=42) |

**Note on ci.group.size:** Stata defaults to 1 (no variance estimation) and bumps to 2 when `estimatevariance` is specified. R defaults to 2 always. This means R computes variance by default while Stata does not. This is a deliberate design choice documented in Stata -- users opt in to variance estimation.

**Note on seed:** Stata uses a deterministic default seed (42) for reproducibility. R uses a random seed. Both approaches are valid; the Stata approach is more reproducible-by-default.

### 2.2 causal_forest()

| R Parameter | R Default | Stata Option | Stata Default | Match? |
|---|---|---|---|---|
| `X` | required | varlist (indep) | required | Yes |
| `Y` | required | varlist (dep) | required | Yes |
| `W` | required | varlist (treat) | required | Yes |
| `Y.hat` | NULL | Internal nuisance | Auto-computed | Yes |
| `W.hat` | NULL | Internal nuisance | Auto-computed | Yes |
| `num.trees` | 2000 | `ntrees()` | 2000 | Yes |
| `sample.weights` | NULL | `weights()` | none | Yes |
| `clusters` | NULL | `cluster()` | none | Yes |
| `equalize.cluster.weights` | FALSE | -- | -- | **Missing** |
| `sample.fraction` | 0.5 | `samplefrac()` | 0.5 | Yes |
| `mtry` | auto | `mtry()` | 0 (auto) | Yes |
| `min.node.size` | 5 | `minnodesize()` | 5 | Yes |
| `honesty` | TRUE | `nohonesty` | TRUE | Yes |
| `honesty.fraction` | 0.5 | `honestyfrac()` | 0.5 | Yes |
| `honesty.prune.leaves` | TRUE | `nohonestyprune` | TRUE | Yes |
| `alpha` | 0.05 | `alpha()` | 0.05 | Yes |
| `imbalance.penalty` | 0 | `imbalancepenalty()` | 0.0 | Yes |
| `stabilize.splits` | TRUE | `nostabilizesplits` | TRUE | Yes |
| `ci.group.size` | 2 | `cigroupsize()` | 1 | **Mismatch** (see note above) |
| `tune.parameters` | "none" | `grf_tune` (separate) | N/A | Different approach |
| `tune.num.trees` | 200 | `grf_tune tunetrees()` | 200 | Separate command |
| `tune.num.reps` | 50 | `grf_tune numreps()` | 50 | Separate command |
| `tune.num.draws` | 1000 | -- | -- | **Missing** |
| `compute.oob.predictions` | TRUE | Always TRUE | TRUE | Yes |
| `num.threads` | NULL | `numthreads()` | 0 | Yes |
| `seed` | random | `seed()` | 42 | **Mismatch** (deterministic default) |

**Additional Stata options:** `nuisancetrees()` (default 500), `yhatgenerate()`, `whatgenerate()`, `estimatevariance`, `vargenerate()`, `replace`, `nomia`. These are Stata-specific conveniences.

### 2.3 quantile_forest()

| R Parameter | R Default | Stata Option | Stata Default | Match? |
|---|---|---|---|---|
| `quantiles` | c(0.1, 0.5, 0.9) | `quantiles()` | "0.1 0.5 0.9" | Yes |
| `regression.splitting` | FALSE | -- | -- | **Missing** |
| All standard params | (same as regression) | (same) | (same) | Yes |

### 2.4 instrumental_forest()

| R Parameter | R Default | Stata Option | Stata Default | Match? |
|---|---|---|---|---|
| `Y.hat, W.hat, Z.hat` | NULL | Internal nuisance | Auto-computed | Yes |
| `reduced.form.weight` | 0 | `reducedformweight()` | 0.0 | Yes |
| `stabilize.splits` | TRUE | `stabilizesplits` (opt-in) | FALSE | **Mismatch** (R=TRUE, Stata=FALSE) |
| All standard params | (same) | (same) | (same) | Yes |

**Note:** Stata's instrumental_forest has `stabilizesplits` as an opt-in flag (default OFF), while R defaults to TRUE. This is a default value mismatch.

### 2.5 probability_forest()

| R Parameter | R Default | Stata Option | Stata Default | Match? |
|---|---|---|---|---|
| `ci.group.size` | 2 | (not exposed) | 1 | **Missing** (always 1) |
| All standard params | (same) | (same) | (same) | Yes |

### 2.6 survival_forest()

| R Parameter | R Default | Stata Option | Stata Default | Match? |
|---|---|---|---|---|
| `failure.times` | NULL | `numfailures()` (count) | 0 (auto) | **Partial** (count only, not explicit times) |
| `num.trees` | 1000 | `ntrees()` | 2000 | **Mismatch** (R=1000, Stata=2000) |
| `min.node.size` | 15 | `minnodesize()` | 15 | Yes |
| `prediction.type` | "Kaplan-Meier" | `predtype()` | 1 (KM) | Yes |
| `fast.logrank` | FALSE | -- | -- | **Missing** |
| `sample.weights` | NULL | `weights()` | none | Yes |
| All other standard params | (same) | (same) | (same) | Yes |

### 2.7 causal_survival_forest()

| R Parameter | R Default | Stata Option | Stata Default | Match? |
|---|---|---|---|---|
| `target` | "RMST" | `target()` | 1 (RMST) | Yes |
| `horizon` | NULL | `horizon()` | 0 (auto-median) | Yes |
| `failure.times` | NULL | -- | -- | **Missing** |
| `W.hat` | NULL | Internal nuisance | Auto-computed | Yes |
| `fast.logrank` | FALSE | -- | -- | **Missing** |
| `stabilize.splits` | TRUE | `nostabilizesplits` | TRUE | Yes |
| `tune.parameters` | "none" | `grf_tune` (separate) | N/A | Different approach |
| Nuisance pipeline | Full IPCW | Simplified IPCW | Approximate | **Partial** (see note) |

**Note:** The causal survival forest nuisance pipeline uses a simplified IPCW approximation rather than the full censoring survival function estimation that R implements. This is explicitly documented in the source code with a recommendation to pre-compute nuisance columns in R for exact results. The `numer()` and `denom()` options allow users to supply pre-computed nuisance estimates.

### 2.8 multi_arm_causal_forest()

| R Parameter | R Default | Stata Option | Stata Default | Match? |
|---|---|---|---|---|
| `Y.hat, W.hat` | NULL | Internal nuisance | Auto-computed | Yes |
| `stabilize.splits` | TRUE | `nostabilizesplits` | TRUE | Yes |
| All standard params | (same) | (same) | (same) | Yes |

**Note:** R accepts W as a factor with K levels. Stata accepts K-1 binary treatment indicator columns via `ntreat()`. The interface is adapted for Stata conventions but functionally equivalent.

### 2.9 multi_regression_forest()

| R Parameter | R Default | Stata Option | Stata Default | Match? |
|---|---|---|---|---|
| All standard params | (same) | (same) | (same) | Yes |

**Note:** Variance estimation is correctly documented as unsupported, matching R's behavior.

### 2.10 ll_regression_forest()

| R Parameter | R Default | Stata Option | Stata Default | Match? |
|---|---|---|---|---|
| `enable.ll.split` | FALSE | `llsplit` | FALSE | Yes |
| `ll.split.weight.penalty` | FALSE | `llweightpenalty` | FALSE | Yes |
| `ll.split.lambda` | 0.1 | `lllambda()` | 0.1 | Yes |
| `ll.split.variables` | NULL (all) | -- | -- | **Partial** (boolean toggle, not per-variable) |
| `ll.split.cutoff` | NULL | `llcutoff()` | 0 | Yes |
| All standard params | (same) | (same) | (same) | Yes |

**Known limitation (documented in source):** R's `ll.split.variables` accepts a vector of variable indices. Stata's `llsplit` is a boolean toggle only.

### 2.11 boosted_regression_forest()

| R Parameter | R Default | Stata Option | Stata Default | Match? |
|---|---|---|---|---|
| `boost.steps` | NULL (auto) | `booststeps()` | 0 (auto) | Yes |
| `boost.error.reduction` | 0.97 | `boosterrorreduction()` | 0.97 | Yes |
| `boost.max.steps` | 5 | `boostmaxsteps()` | 5 | Yes |
| `boost.trees.tune` | 10 | `boosttreestune()` | 10 | Yes |
| All standard params | (same) | (same) | (same) | Yes |

### 2.12 lm_forest()

| R Parameter | R Default | Stata Option | Stata Default | Match? |
|---|---|---|---|---|
| `Y.hat, W.hat` | NULL | Internal nuisance | Auto-computed | Yes |
| `gradient.weights` | NULL | -- | -- | **Missing** |
| `stabilize.splits` | FALSE | `nostabilizesplits` | TRUE | **Mismatch** (R=FALSE, Stata=TRUE) |
| All standard params | (same) | (same) | (same) | Yes |

**Note:** lm_forest `stabilize.splits` default mismatch. R defaults to FALSE, Stata defaults to TRUE. The source code comments acknowledge this: "R's lm_forest defaults stabilize.splits = FALSE, but we keep the consistent noSTABilizesplits interface (default ON)."

### 2.13 average_treatment_effect()

| R Parameter | R Default | Stata Option | Stata Default | Match? |
|---|---|---|---|---|
| `forest` | required | Uses e() from causal forest | required | Yes |
| `target.sample` | "all" | `targetsample()` | "all" | Yes |
| `method` | "AIPW" | AIPW only | AIPW | **Partial** (no TMLE option) |
| `subset` | NULL | `if`/`in` | none | Yes |
| `debiasing.weights` | NULL | -- | -- | **Missing** |
| `compliance.score` | NULL | -- | -- | **Missing** |
| `num.trees.for.weights` | 500 | -- | -- | **Missing** |

### 2.14 best_linear_projection()

| R Parameter | R Default | Stata Option | Stata Default | Match? |
|---|---|---|---|---|
| `A` | NULL | varlist | e(indepvars) | Yes |
| `subset` | NULL | `if`/`in` | none | Yes |
| `debiasing.weights` | NULL | -- | -- | **Missing** |
| `compliance.score` | NULL | -- | -- | **Missing** |
| `num.trees.for.weights` | 500 | -- | -- | **Missing** |
| `vcov.type` | "HC3" | `vcovtype()` | "HC3" | Yes |
| `target.sample` | "all"/"overlap" | `targetsample()` | "all" | **Extended** (Stata also supports "treated"/"control") |

### 2.15 test_calibration()

| R Parameter | R Default | Stata Option | Stata Default | Match? |
|---|---|---|---|---|
| `forest` | required | Uses e() | required | Yes |
| `vcov.type` | "HC3" | HC3 hardcoded | "HC3" | Yes (but not configurable) |

### 2.16 variable_importance()

| R Parameter | R Default | Stata Option | Stata Default | Match? |
|---|---|---|---|---|
| `forest` | required | varlist (refits internally) | required | **Different** (refits) |
| `decay.exponent` | 2 | -- | -- | **Missing** (hardcoded) |
| `max.depth` | 4 | `maxdepth()` | 4 | Yes |

**Note:** R's `variable_importance()` operates on an existing forest object. Stata's `grf_variable_importance` fits a new regression forest internally, since the plugin architecture doesn't persist forest objects. The `decay.exponent` parameter is not exposed.

### 2.17 rank_average_treatment_effect()

| R Parameter | R Default | Stata Option | Stata Default | Match? |
|---|---|---|---|---|
| `forest` | required | Uses e() | required | Yes |
| `priorities` | required | varname | required | Yes |
| `target` | c("AUTOC","QINI") | `target()` | "AUTOC" | Yes |
| `q` | seq(0.1,1,by=0.1) | `quantiles()` | same grid | Yes |
| `R` | 200 | `bootstrap()` | 200 | Yes |
| `subset` | NULL | `if`/`in` | none | Yes |
| `debiasing.weights` | NULL | -- | -- | **Missing** |
| `compliance.score` | NULL | `compliancescore()` | none | Yes |
| `num.trees.for.weights` | 500 | -- | -- | **Missing** |

### 2.18 get_scores()

| R Method | Stata Support | Status |
|---|---|---|
| `get_scores.causal_forest` | `grf_get_scores` | **Implemented** |
| `get_scores.multi_arm_causal_forest` | `grf_get_scores` | **Implemented** |
| `get_scores.causal_survival_forest` | `grf_get_scores` (plug-in scores) | **Partial** (not true AIPW) |
| `get_scores.instrumental_forest` | `grf_get_scores` | **Implemented** |

---

## 3. Issues Found

### 3.1 Critical Issues (None)

No critical issues found. The package correctly implements the core GRF algorithms via the same C++ backend.

### 3.2 High Severity Issues

| # | Issue | Severity | Details |
|---|---|---|---|
| H1 | `equalize.cluster.weights` not exposed | High | R's parameter to equalize weighting across clusters of different sizes is missing from all forest types. Users with unbalanced clusters may get different results. |
| H2 | Causal survival nuisance pipeline is approximate | High | The full IPCW censoring hazard pipeline is not implemented; uses a simplified approximation. Explicitly documented with workaround (pre-compute in R). |
| H3 | `debiasing.weights` not exposed in ATE/BLP/RATE | High | Advanced users cannot supply custom debiasing weights for doubly-robust inference. Affects `grf_ate`, `grf_best_linear_projection`, and `grf_rate`. |

### 3.3 Medium Severity Issues

| # | Issue | Severity | Details |
|---|---|---|---|
| M1 | `instrumental_forest` `stabilize.splits` default mismatch | Medium | R defaults to TRUE, Stata defaults to FALSE. Opt-in vs opt-out. |
| M2 | `lm_forest` `stabilize.splits` default mismatch | Medium | R defaults to FALSE, Stata defaults to TRUE. Opposite of R. |
| M3 | `survival_forest` `num.trees` default mismatch | Medium | R defaults to 1000, Stata defaults to 2000. Stata uses more trees. |
| M4 | `method="TMLE"` not available in `grf_ate` | Medium | R supports TMLE as alternative to AIPW. Stata only supports AIPW. |
| M5 | `tune.num.draws` not exposed | Medium | The tuning random search draw count is not configurable. |
| M6 | `regression.splitting` for quantile forest not exposed | Medium | R allows switching to regression-based splitting. |
| M7 | `gradient.weights` for lm_forest not exposed | Medium | Advanced parameter for gradient weighting. |
| M8 | `fast.logrank` not exposed for survival/causal_survival | Medium | Performance optimization parameter missing. |
| M9 | `variable_importance` `decay.exponent` not exposed | Medium | Hardcoded rather than configurable. |
| M10 | `variable_importance` refits forest instead of using existing | Medium | Different from R API (which operates on an existing forest object). |

### 3.4 Low Severity Issues

| # | Issue | Severity | Details |
|---|---|---|---|
| L1 | `ci.group.size` defaults to 1 (not 2) | Low | Deliberate design choice for Stata (opt-in variance estimation). |
| L2 | `seed` defaults to 42 (not random) | Low | Deliberate choice for reproducibility. Acceptable. |
| L3 | `ll.split.variables` is boolean, not per-variable | Low | Documented limitation in source code. |
| L4 | `test_calibration` `vcov.type` not configurable | Low | Always uses HC3 (matches R default). |
| L5 | `num.trees.for.weights` not exposed | Low | Internal parameter for debiasing weight computation. |
| L6 | `failure.times` only as count, not explicit times | Low | `numfailures()` specifies count; R allows explicit time vector. |
| L7 | Tree introspection not available | Low | `get_tree()`, `get_leaf_node()`, `get_forest_weights()` -- not feasible in plugin architecture. |
| L8 | `merge_forests()` not available | Low | Forest objects are not persisted between calls. |
| L9 | `grf_options()` not implemented | Low | No global option registry needed for Stata's per-command design. |
| L10 | Plot functions not implemented | Low | Stata uses different graphing conventions; users can graph manually. |

---

## 4. Test Coverage Assessment

### 4.1 Forest Type Test Coverage

| Forest Type | Basic Test | Options Test | Fidelity Test | Cluster Test | Weight Test | MIA Test | Seed Test |
|---|---|---|---|---|---|---|---|
| regression | Yes | Yes (16 tests) | Yes (corr>0.95) | Yes | Yes | Yes | Yes |
| causal | Yes | Yes (19 tests) | Yes (corr>0.90, ATE z<3) | Yes | Yes | Yes | Yes |
| quantile | Yes | Yes | Yes (corr>0.90) | Yes | Yes | Yes | Yes |
| instrumental | Yes | Yes | Yes (corr>0.90) | Yes | Yes | Yes | Yes |
| probability | Yes | Yes | Yes (corr>0.90) | Yes | Yes | Yes | -- |
| survival | Yes | Yes | Yes (corr>0.90) | Yes | Yes | Yes | Yes |
| causal_survival | Yes | Yes | Yes (corr>0.90) | Yes | -- | Yes | -- |
| multi_arm_causal | Yes | Yes | Yes (corr>0.90) | Yes | -- | Yes | -- |
| multi_regression | Yes | Yes | Yes (corr>0.90) | Yes | -- | Yes | -- |
| ll_regression | Yes | Yes | Yes (corr>0.90) | Yes | -- | Yes | -- |
| boosted_regression | Yes | Yes | Yes (corr>0.90) | Yes | -- | Yes | -- |
| lm_forest | Yes | Yes | Yes (corr>0.90) | -- | -- | Yes | -- |

### 4.2 Post-Estimation Test Coverage

| Command | Basic Test | Options Test | Fidelity Test |
|---|---|---|---|
| `grf_ate` | Yes | Yes (4 target.sample variants) | Yes (z<3 vs R, all 4 targets) |
| `grf_best_linear_projection` | Yes | Yes (vcov types, target samples) | Yes (per-coef z<3, HC0/HC3) |
| `grf_test_calibration` | Yes | Yes | Yes (per-coef z<3 vs R) |
| `grf_variable_importance` | Yes | Yes | Yes (Spearman > 0.70) |
| `grf_rate` | Yes | Yes (AUTOC/QINI, bootstrap) | Yes (z<3 vs R, AUTOC) |
| `grf_get_scores` | Yes (causal, multi-arm, instrumental, csf) | Yes | Yes (corr>0.90 vs R) |
| `grf_predict` | Yes | Yes | -- |
| `grf_tune` | Yes (all forest types) | Yes (extended) | -- |

### 4.3 Cross-Cutting Feature Tests

| Feature | Test File | Coverage |
|---|---|---|
| Cluster-robust estimation | test_clusters.do | 7+ tests across forest types, ATE, BLP, calibration |
| Sample weights | test_weights.do | 5+ tests: regression, causal, produces different results |
| Missing data (MIA) | test_mia.do | 6+ tests: MIA enabled vs nomia, all obs used with MIA |
| Seed reproducibility | test_seed_reproducibility.do | 6+ tests: same seed => identical, different seed => different |
| Tuning | test_tune.do, test_tune_extended.do | 8+ tests: all forest types, apply tuned params |
| Fidelity vs R | test_fidelity.do | 27 tests: predictions, ATE, BLP, calibration, scores, VI |
| Full pipeline | test_full_pipeline.do | End-to-end workflow test |
| Data generation | test_generate_data.do | DGP validation |

### 4.4 Test Count Summary

- **Total test files:** 43
- **Forest option test files:** 12 (one per forest type)
- **Post-estimation test files:** 2 (basic + options)
- **Fidelity test comparisons:** 27 (vs R reference data)
- **Cross-cutting test files:** 7 (clusters, weights, MIA, seed, tune, features, pipeline)

---

## 5. Default Value Alignment Summary

Across all 12 forest types, the parameter defaults are audited below:

| Parameter | R Default | Stata Default | Aligned? |
|---|---|---|---|
| `num.trees` | 2000 (1000 for survival) | 2000 (all types) | Partial (survival mismatch) |
| `sample.fraction` | 0.5 | 0.5 | Yes |
| `mtry` | auto formula | 0 (auto in C++) | Yes |
| `min.node.size` | 5 (15 for survival/csf) | 5 (15 for survival, 15 for csf) | Yes |
| `honesty` | TRUE | TRUE | Yes |
| `honesty.fraction` | 0.5 | 0.5 | Yes |
| `honesty.prune.leaves` | TRUE | TRUE | Yes |
| `alpha` | 0.05 | 0.05 | Yes |
| `imbalance.penalty` | 0 | 0.0 | Yes |
| `ci.group.size` | 2 | 1 | No (deliberate) |
| `stabilize.splits` | TRUE (causal/csf/mac) | TRUE (causal/csf/mac) | Mostly (instrumental/lm mismatch) |
| `seed` | random | 42 | No (deliberate) |
| `num.threads` | NULL (all) | 0 (all) | Yes |

---

## 6. Documentation Quality

### 6.1 Help Files

Every forest type and post-estimation command has a dedicated `.sthlp` file:
- 12 forest help files
- 7 post-estimation help files
- 1 overview help file (`grf.sthlp`)

### 6.2 README.md

Comprehensive README with:
- Installation instructions (multi-platform)
- Forest type table
- Post-estimation command table
- Quick start example
- Full parameter reference

### 6.3 Source Code Documentation

Each `.ado` file includes:
- Version header with citation
- Comments explaining the nuisance estimation pipeline
- Clear variable ordering documentation for plugin calls
- Known limitations documented in-line

---

## 7. Architecture Assessment

### 7.1 Plugin Design

The architecture uses a single C++ plugin (`grf_plugin.cpp`, 1564 lines) that:
- Statically links the full grf C++ library
- Dispatches on `forest_type` string argument
- Handles all 12 forest types + variable importance
- Manages train/test splitting for `grf_predict`
- Supports cluster and weight column passthrough

### 7.2 Nuisance Estimation

Forest types requiring orthogonalization (causal, instrumental, multi-arm, lm, causal_survival) implement the nuisance pipeline in Stata:
- Fit regression forests for Y.hat, W.hat (and Z.hat for instrumental)
- Center outcomes and treatments
- Pass centered data to the C++ plugin

This is correct and matches R's internal pipeline.

### 7.3 Post-Estimation

Post-estimation commands read from `e()` results stored by forest commands:
- `grf_ate`: Computes AIPW DR scores in Stata, supports cluster-robust SE
- `grf_best_linear_projection`: OLS on DR scores with HC0-HC3 via Mata
- `grf_test_calibration`: Calibration regression with HC3
- `grf_rate`: Bootstrap RATE computation in pure Stata
- `grf_get_scores`: DR score computation for causal/instrumental/multi-arm/csf

---

## 8. Overall Fidelity Score

### Score: 8.5 / 10

### Justification

**Strengths (pushing toward 10):**
- All 12 R forest types are implemented with the same C++ backend
- Core parameters (num.trees, mtry, min.node.size, sample.fraction, honesty, alpha, imbalance.penalty) are faithfully mapped across all forest types
- Nuisance estimation pipelines correctly implement orthogonalization
- All 7 post-estimation functions are implemented (ATE, BLP, calibration, VI, RATE, scores, predict)
- Comprehensive test suite: 43 test files, 27 R-vs-Stata fidelity comparisons
- Test fidelity thresholds are appropriate (correlation > 0.90, z-test < 3)
- MIA (Missing Indicator Action) support is a notable feature matching R's native handling
- Cluster and sample weight support throughout
- Seed reproducibility verified
- Separate tuning command (`grf_tune`) covers all forest types
- Well-documented with help files, README, and inline comments

**Weaknesses (pulling below 10):**
- `equalize.cluster.weights` missing across all forest types (-0.3)
- Causal survival nuisance pipeline is approximate, not exact (-0.3)
- `debiasing.weights` not exposed in post-estimation commands (-0.2)
- Several default value mismatches (instrumental stabilize.splits, lm stabilize.splits, survival num.trees) (-0.2)
- TMLE method not available for ATE (-0.1)
- A handful of advanced parameters not exposed (gradient.weights, regression.splitting, fast.logrank, decay.exponent, tune.num.draws) (-0.2)
- Tree introspection functions (get_tree, get_leaf_node, get_forest_weights, merge_forests) not available (-0.1) -- architectural limitation
- Variable importance refits rather than using existing forest (-0.1)

### Comparison with Perfect Score

A perfect 10 would require:
1. Every single R parameter exposed with matching defaults
2. Tree introspection capabilities
3. Full causal survival nuisance pipeline
4. TMLE support
5. Forest object persistence for `merge_forests()` and `variable_importance()` on existing forests

Given the inherent constraints of Stata's plugin architecture (no persistent C++ objects between calls), achieving items 2, 5, and parts of 3 would require fundamental architectural changes. The package makes intelligent trade-offs within these constraints.

---

## 9. Recommendations

### Priority 1 (High Impact, Moderate Effort)
1. **Add `equalize.cluster.weights` option** to all forest commands. This is a straightforward boolean passthrough to the C++ backend.
2. **Fix `stabilize.splits` defaults** for `instrumental_forest` (should default ON) and document the `lm_forest` deviation clearly in the help file.
3. **Fix `survival_forest` num.trees default** to 1000 to match R.

### Priority 2 (Medium Impact)
4. **Expose `debiasing.weights`** in `grf_ate`, `grf_best_linear_projection`, and `grf_rate`.
5. **Add `regression.splitting` option** to quantile forest.
6. **Expose `fast.logrank`** for survival and causal survival forests.
7. **Expose `decay.exponent`** in `grf_variable_importance`.

### Priority 3 (Low Impact / Nice-to-Have)
8. **Add TMLE method** to `grf_ate` (significant implementation effort).
9. **Expose `gradient.weights`** in `grf_lm_forest`.
10. **Add `ll.split.variables` vector support** to local linear regression forest.
11. **Make `test_calibration` `vcov.type` configurable**.

---

*End of audit.*
