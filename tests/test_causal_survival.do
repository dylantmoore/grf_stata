* test_causal_survival.do -- Test grf_causal_survival_forest
* Run from the test-c-plugin-skill-grf/ directory

clear all
set more off

display as text ""
display as text "=============================================="
display as text " GRF Causal Survival Forest Tests"
display as text "=============================================="

* ---- Test 1: Basic functionality (synthetic data) ----
display as text ""
display as text "--- Test 1: Basic functionality ---"

clear
set obs 500
set seed 42

* Generate test data with treatment effect on survival
gen x1 = rnormal()
gen x2 = rnormal()
gen x3 = rnormal()
gen w = (runiform() > 0.5)
gen double time = exp(0.5 * x1 + 0.3 * w + rnormal())
gen status = (runiform() > 0.3)

* Run causal survival forest
grf_causal_survival_forest time status w x1 x2 x3, gen(cate_cs) ntrees(500) seed(42) replace

* Check predictions exist and are non-missing
quietly count if !missing(cate_cs)
local n_pred = r(N)
display as text "  CATE predictions: " as result `n_pred' " / 500 non-missing"
assert `n_pred' > 400

* Check reasonable range for causal survival effects
summarize cate_cs
display as text "  Mean CATE: " as result %9.4f r(mean)
display as text "  SD CATE:   " as result %9.4f r(sd)

* Treatment effect should have some variation
assert r(sd) > 0

* e() results should be stored
assert "`e(cmd)'" == "grf_causal_survival_forest"
assert "`e(forest_type)'" == "causal_survival"
assert e(N) == 500
assert e(n_trees) == 500
assert e(horizon) > 0

display as text "  PASSED"

* ---- Test 2: Pre-computed nuisance columns ----
display as text ""
display as text "--- Test 2: Pre-computed nuisance columns ---"

* Create fake nuisance columns (simplified for testing)
gen double my_numer = (w - 0.5) * status * min(time, 2.0)
gen double my_denom = 1.0

grf_causal_survival_forest time status w x1 x2 x3, gen(cate_pre) ///
    ntrees(500) seed(42) numer(my_numer) denom(my_denom) replace

quietly count if !missing(cate_pre)
local n_pred = r(N)
display as text "  CATE predictions (precomputed): " as result `n_pred' " / 500 non-missing"
assert `n_pred' > 400

* Should still produce valid predictions
summarize cate_pre
assert r(sd) > 0

display as text "  PASSED"

* ---- Test 3: Options ----
display as text ""
display as text "--- Test 3: Options (horizon, nohonesty) ---"

grf_causal_survival_forest time status w x1 x2 x3, gen(cate_opts) ///
    ntrees(200) seed(42) horizon(3.0) nohonesty replace

quietly count if !missing(cate_opts)
assert r(N) > 400

* Horizon should be stored
assert abs(e(horizon) - 3.0) < 0.001

display as text "  PASSED"

* ---- Test 4: if/in restrictions ----
display as text ""
display as text "--- Test 4: if/in restrictions ---"

grf_causal_survival_forest time status w x1 x2 x3 if x1 > 0, ///
    gen(cate_sub) ntrees(200) seed(42) replace

quietly count if !missing(cate_sub) & x1 > 0
local n_pred = r(N)
quietly count if x1 > 0
local n_eligible = r(N)
display as text "  Predictions: " as result `n_pred' " / " as result `n_eligible' " (if x1 > 0)"
assert `n_pred' > 0

* Predictions should only exist for x1 > 0
quietly count if !missing(cate_sub) & x1 <= 0
assert r(N) == 0
display as text "  No predictions for x1 <= 0: OK"

display as text "  PASSED"

* ---- Summary ----
display as text ""
display as text "=============================================="
display as text " All causal survival forest tests completed"
display as text "=============================================="
display as text ""
