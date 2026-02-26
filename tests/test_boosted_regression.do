* test_boosted_regression.do -- Test grf_boosted_regression_forest
* Run from the test-c-plugin-skill-grf/ directory

clear all
set more off

display as text ""
display as text "=============================================="
display as text " GRF Boosted Regression Forest Tests"
display as text "=============================================="

* ---- Test 1: Basic functionality (auto-tune) ----
display as text ""
display as text "--- Test 1: Basic functionality (auto-tune) ---"

clear
set obs 500
set seed 12345

gen x1 = rnormal()
gen x2 = rnormal()
gen x3 = rnormal()
gen y = 3*x1 + 2*x2^2 + x1*x3 + rnormal()

grf_boosted_regression_forest y x1 x2 x3, gen(pred_boost) ntrees(500) seed(42)

quietly count if !missing(pred_boost)
local n_pred = r(N)
display as text "  Predictions: " as result `n_pred' " / 500 non-missing"
assert `n_pred' > 400

summarize pred_boost
display as text "  Mean pred: " as result %9.4f r(mean)
assert r(sd) > 0

assert "`e(cmd)'" == "grf_boosted_regression_forest"
assert "`e(forest_type)'" == "boosted_regression"
assert e(N) == 500
assert e(boost_steps) >= 1

display as text "  Boosting steps used: " as result e(boost_steps)
display as text "  PASSED"

* ---- Test 2: Fixed number of steps ----
display as text ""
display as text "--- Test 2: Fixed 3 boosting steps ---"

grf_boosted_regression_forest y x1 x2 x3, gen(pred_fixed) ///
    ntrees(500) seed(42) booststeps(3) replace

quietly count if !missing(pred_fixed)
assert r(N) > 400

assert e(boost_steps) == 3

display as text "  PASSED"

* ---- Test 3: Boosted should improve over single forest ----
display as text ""
display as text "--- Test 3: Boosted vs single forest ---"

grf_regression_forest y x1 x2 x3, gen(pred_single) ntrees(500) seed(42) replace

* Compute RMSE for boosted vs single
quietly gen double resid_boost = (y - pred_boost)^2 if !missing(pred_boost)
quietly gen double resid_single = (y - pred_single)^2 if !missing(pred_single)
quietly summarize resid_boost
local rmse_boost = sqrt(r(mean))
quietly summarize resid_single
local rmse_single = sqrt(r(mean))

display as text "  RMSE single forest: " as result %9.4f `rmse_single'
display as text "  RMSE boosted:       " as result %9.4f `rmse_boost'

* Boosted should be at least as good (possibly better)
* Allow some tolerance since OOB predictions can vary
assert `rmse_boost' < `rmse_single' * 1.1

display as text "  PASSED"

* ---- Test 4: Options (nohonesty, max steps) ----
display as text ""
display as text "--- Test 4: Options ---"

grf_boosted_regression_forest y x1 x2 x3, gen(pred_opts) ///
    ntrees(200) seed(42) nohonesty boostmaxsteps(3) ///
    boosterrorreduction(0.5) replace

quietly count if !missing(pred_opts)
assert r(N) > 400
assert e(honesty) == 0
assert e(boost_steps) <= 3

display as text "  PASSED"

* ---- Test 5: if/in restrictions ----
display as text ""
display as text "--- Test 5: if/in restrictions ---"

grf_boosted_regression_forest y x1 x2 x3 if x1 > 0, gen(pred_sub) ///
    ntrees(200) seed(42) booststeps(2) replace

quietly count if !missing(pred_sub) & x1 > 0
local n_pred = r(N)
assert `n_pred' > 0

quietly count if !missing(pred_sub) & x1 <= 0
assert r(N) == 0

display as text "  PASSED"

* ---- Summary ----
display as text ""
display as text "=============================================="
display as text " All boosted regression forest tests completed"
display as text "=============================================="
display as text ""
