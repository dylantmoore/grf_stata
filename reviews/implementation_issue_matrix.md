# Implementation Issue Matrix (Consolidated Report)

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
