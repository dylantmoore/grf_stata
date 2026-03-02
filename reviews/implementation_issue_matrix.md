# Implementation Issue Matrix (Consolidated Report)

## Current baseline scope (requested gaps 4/5/6)

| Gap | Resolution Type | Status | Primary Files |
|---|---|---|---|
| 4. `orthog.boosting` | Manifest-backed closure | Complete (non-applicable in upstream v2.5.0) | reviews/r_api_manifest.json, tools/extract_r_api_manifest.R, tools/check_parity_scope.R |
| 5. `honesty.prune.method` enum | Manifest-backed closure | Complete (non-applicable in upstream v2.5.0) | reviews/r_api_manifest.json, tools/extract_r_api_manifest.R, tools/check_parity_scope.R |
| 6. causal-survival nuisance parity + score construction | Implement | Complete | grf_causal_survival_forest.ado, grf_get_scores.ado, tests/test_options_causal_survival.do, tests/test_get_scores.do, tests/test_fidelity.do |
| Advanced introspection parity | Partial implement + document | Complete | grf_forest_summary.ado, grf_get_forest_weights.ado, tests/test_forest_utilities.do, reviews/introspection_discrepancies.md |

## Consolidated report scope (1-16)

| Gap | Resolution Type | Status | Primary Files |
|---|---|---|---|
| 1. Fidelity consumption | Verify-and-lock | Complete | tests/test_fidelity.do |
| 2. RATE fidelity data | Verify-and-lock | Complete | tests/generate_reference.R, tests/test_fidelity.do |
| 3. Clusters/sample weights | Verify-and-lock | Complete | grf_*_forest.ado, README.md, grf.sthlp |
| 4. MIA passthrough | Verify-and-lock | Complete | grf_plugin.cpp, tests/test_mia.do |
| 5. Fidelity for remaining forests | Verify-and-lock | Complete | tests/test_fidelity.do |
| 6. Tune option coverage | Verify-and-lock | Complete | tests/test_tune_extended.do |
| 7. APE numtreesvariance coverage | Verify-and-lock | Complete | grf_average_partial_effect.ado, tests/test_average_partial_effect.do |
| 8. BLP target.sample + docs | Verify-and-lock + docs | Complete | grf_best_linear_projection.ado, grf_best_linear_projection.sthlp |
| 9. get_scores causal-survival parity | Implement | Complete | grf_get_scores.ado, tests/test_get_scores.do |
| 10. ref path consistency | Implement | Complete | tests/*.do |
| 11. APE deprecation docs | Verify-and-lock | Complete | grf_average_partial_effect.ado, grf_average_partial_effect.sthlp, README.md, grf.sthlp |
| 12. `in` qualifier coverage (RATE/APE) | Verify-and-lock | Complete | tests/test_options_post_estimation.do, tests/test_average_partial_effect.do |
| 13. Survival DGP error handling | Implement | Complete | tests/test_generate_data.do |
| 14. Seed reproducibility runner | Verify-and-lock + runner | Complete | tests/test_seed_reproducibility.do, tests/run_all.do |
| 15. llvars parity docs | Verify-and-lock | Complete | grf_ll_regression_forest.ado, grf_ll_regression_forest.sthlp |
| 16. User nuisance estimates (causal-survival extension) | Implement (partial) | Complete | grf_causal_survival_forest.ado, grf_causal_survival_forest.sthlp, tests/test_options_causal_survival.do |
