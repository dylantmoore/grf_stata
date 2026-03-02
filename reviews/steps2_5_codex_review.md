Not LGTM.

1. **Critical: clustered training is effectively broken (`samples_per_cluster` is hardcoded to 0).**  
   In the plugin, `samples_per_cluster` is initialized to `0` and passed into `ForestOptions`, and GRF’s sampler then draws `0` observations per sampled cluster.  
   [grf_plugin.cpp:396](/private/tmp/grf_stata/grf_plugin.cpp:396) [grf_plugin.cpp:442](/private/tmp/grf_stata/grf_plugin.cpp:442) [RandomSampler.cpp:101](/private/tmp/grf_stata/vendor/grf/core/src/sampling/RandomSampler.cpp:101) [RandomSampler.cpp:107](/private/tmp/grf_stata/vendor/grf/core/src/sampling/RandomSampler.cpp:107)

2. **Critical: sample weights are parsed but never propagated to GRF weighting.**  
   `sample_weights` is read/validated, but never attached to `grf::Data` (no `set_weight_index`) and never passed to trainer APIs; therefore Step 4 is functionally incomplete.  
   [grf_plugin.cpp:409](/private/tmp/grf_stata/grf_plugin.cpp:409) [grf_plugin.cpp:421](/private/tmp/grf_stata/grf_plugin.cpp:421) [grf_plugin.cpp:354](/private/tmp/grf_stata/grf_plugin.cpp:354) [grf_plugin.cpp:384](/private/tmp/grf_stata/grf_plugin.cpp:384)

3. **High: cluster/weight columns are added to the training matrix and can be used as split features.**  
   The `.ado` files append cluster/weight vars into plugin varlists; plugin includes all non-output columns in `Data`; only Y/W/Z are marked disallowed. Tree splitting samples from all non-disallowed columns, so cluster/weight can leak into splits (inconsistent with R `grf`).  
   [grf_regression_forest.ado:143](/private/tmp/grf_stata/grf_regression_forest.ado:143) [grf_plugin.cpp:227](/private/tmp/grf_stata/grf_plugin.cpp:227) [grf_plugin.cpp:361](/private/tmp/grf_stata/grf_plugin.cpp:361) [TreeTrainer.cpp:134](/private/tmp/grf_stata/vendor/grf/core/src/tree/TreeTrainer.cpp:134)

4. **High: MIA passthrough is likely defeated in wrappers by `marksample touse` usage.**  
   With standard Stata behavior, `marksample` drops varlist missings before plugin call, so missing-X rows may never reach the new NaN/MIA logic. This pattern appears across all 12 forest commands.  
   [grf_regression_forest.ado:89](/private/tmp/grf_stata/grf_regression_forest.ado:89) [grf_causal_forest.ado:106](/private/tmp/grf_stata/grf_causal_forest.ado:106) [grf_survival_forest.ado:92](/private/tmp/grf_stata/grf_survival_forest.ado:92)

5. **High: multi-arm DR score formula does not match GRF’s multi-treatment score construction.**  
   Current code uses per-arm scalar `Var(W_j - Ŵ_j)` correction, but GRF’s multi-causal relabeling uses the full inverse covariance matrix of centered treatment residuals (`A^{-1}`), i.e., cross-arm covariance terms matter.  
   [grf_get_scores.ado:55](/private/tmp/grf_stata/grf_get_scores.ado:55) [grf_get_scores.ado:84](/private/tmp/grf_stata/grf_get_scores.ado:84) [MultiCausalRelabelingStrategy.cpp:77](/private/tmp/grf_stata/vendor/grf/core/src/relabeling/MultiCausalRelabelingStrategy.cpp:77) [MultiCausalRelabelingStrategy.cpp:87](/private/tmp/grf_stata/vendor/grf/core/src/relabeling/MultiCausalRelabelingStrategy.cpp:87)

6. **High: causal-survival “DR scores” are not actually DR scores.**  
   The implementation computes `w_resid` and variance, then sets score equal to `tau_hat` only. This is not the documented correction and not consistent with GRF causal-survival nuisance structure.  
   [grf_get_scores.ado:158](/private/tmp/grf_stata/grf_get_scores.ado:158) [grf_get_scores.ado:165](/private/tmp/grf_stata/grf_get_scores.ado:165) [grf_get_scores.ado:173](/private/tmp/grf_stata/grf_get_scores.ado:173)

7. **Medium: new tests 79–86 have correctness issues (string assertions + weak checks).**  
   Several assertions compare string returns using expression form (`assert e(target_sample) == "treated"` / `assert r(forest_type) == ...`), which is not the usual safe Stata macro comparison form. Also tests named “coefficients differ” only assert non-missingness, not difference.  
   [test_options_post_estimation.do:1540](/private/tmp/grf_stata/tests/test_options_post_estimation.do:1540) [test_options_post_estimation.do:1556](/private/tmp/grf_stata/tests/test_options_post_estimation.do:1556) [test_options_post_estimation.do:1657](/private/tmp/grf_stata/tests/test_options_post_estimation.do:1657) [test_options_post_estimation.do:1735](/private/tmp/grf_stata/tests/test_options_post_estimation.do:1735) [test_options_post_estimation.do:1583](/private/tmp/grf_stata/tests/test_options_post_estimation.do:1583) [test_options_post_estimation.do:1612](/private/tmp/grf_stata/tests/test_options_post_estimation.do:1612)

8. **Medium: cluster test data generation appears incorrect (cluster effect mostly missing).**  
   `cluster_effect` is generated only at `_n==1`, then propagated within clusters from first value; for most clusters that first value is missing. This undermines the intended cluster-structure test.  
   [test_clusters.do:33](/private/tmp/grf_stata/tests/test_clusters.do:33) [test_clusters.do:34](/private/tmp/grf_stata/tests/test_clusters.do:34)