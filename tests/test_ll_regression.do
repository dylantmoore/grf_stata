* test_ll_regression.do -- Test grf_ll_regression_forest
* Run from the test-c-plugin-skill-grf/ directory

clear all
set more off

display as text ""
display as text "=============================================="
display as text " GRF Local Linear Regression Forest Tests"
display as text "=============================================="

* ---- Test 1: Basic functionality (no LL splits) ----
display as text ""
display as text "--- Test 1: Basic functionality (default, no LL splits) ---"

clear
set obs 500
set seed 12345

gen x1 = rnormal()
gen x2 = rnormal()
gen y = 2*x1 + x2^2 + rnormal()

grf_ll_regression_forest y x1 x2, gen(pred_ll) ntrees(500) seed(42)

quietly count if !missing(pred_ll)
local n_pred = r(N)
display as text "  Predictions: " as result `n_pred' " / 500 non-missing"
assert `n_pred' > 400

summarize pred_ll
display as text "  Mean pred: " as result %9.4f r(mean)
assert r(sd) > 0

assert "`e(cmd)'" == "grf_ll_regression_forest"
assert "`e(forest_type)'" == "ll_regression"
assert e(N) == 500
assert e(enable_ll_split) == 0

display as text "  PASSED"

* ---- Test 2: With LL splits enabled ----
display as text ""
display as text "--- Test 2: With LL splits enabled ---"

grf_ll_regression_forest y x1 x2, gen(pred_lls) ntrees(500) seed(42) ///
    llenable lllambda(0.5) replace

quietly count if !missing(pred_lls)
local n_pred = r(N)
display as text "  Predictions (LL split): " as result `n_pred' " / 500 non-missing"
assert `n_pred' > 400

assert e(enable_ll_split) == 1
assert abs(e(ll_lambda) - 0.5) < 0.001

display as text "  PASSED"

* ---- Test 3: Options (nohonesty, lambda, weight penalty) ----
display as text ""
display as text "--- Test 3: Options ---"

grf_ll_regression_forest y x1 x2, gen(pred_opt) ntrees(200) seed(42) ///
    nohonesty lllambda(1.0) llweightpenalty replace

quietly count if !missing(pred_opt)
assert r(N) > 400
assert e(honesty) == 0
assert abs(e(ll_lambda) - 1.0) < 0.001
assert e(ll_weight_penalty) == 1

display as text "  PASSED"

* ---- Test 4: if/in restrictions ----
display as text ""
display as text "--- Test 4: if/in restrictions ---"

grf_ll_regression_forest y x1 x2 if x1 > 0, gen(pred_sub) ///
    ntrees(200) seed(42) replace

quietly count if !missing(pred_sub) & x1 > 0
local n_pred = r(N)
quietly count if x1 > 0
local n_eligible = r(N)
display as text "  Predictions: " as result `n_pred' " / " as result `n_eligible' " (if x1 > 0)"
assert `n_pred' > 0

quietly count if !missing(pred_sub) & x1 <= 0
assert r(N) == 0
display as text "  No predictions for x1 <= 0: OK"

display as text "  PASSED"

* ---- Test 5: Comparison with regression forest ----
display as text ""
display as text "--- Test 5: LL vs standard regression (should differ) ---"

grf_regression_forest y x1 x2, gen(pred_reg) ntrees(500) seed(42) replace

correlate pred_ll pred_reg
display as text "  Correlation LL vs standard: " as result %6.4f r(rho)
* They should be correlated but not identical
assert r(rho) > 0.5
assert r(rho) < 1.0

display as text "  PASSED"

* ---- Summary ----
display as text ""
display as text "=============================================="
display as text " All local linear regression forest tests completed"
display as text "=============================================="
display as text ""
