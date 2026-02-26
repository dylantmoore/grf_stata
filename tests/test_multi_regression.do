* test_multi_regression.do -- Test grf_multi_regression_forest
* Run from the test-c-plugin-skill-grf/ directory

clear all
set more off

display as text ""
display as text "=============================================="
display as text " GRF Multi-Regression Forest Tests"
display as text "=============================================="

* ---- Test 1: Basic functionality (2 outcomes) ----
display as text ""
display as text "--- Test 1: Basic functionality (2 outcomes) ---"

clear
set obs 500
set seed 42

* Generate test data with 2 correlated outcomes
gen x1 = rnormal()
gen x2 = rnormal()
gen x3 = rnormal()
gen double y1 = 2 * x1 + x2 + rnormal()
gen double y2 = x1 + 3 * x2 + 0.5 * rnormal()

* Run multi-regression forest
grf_multi_regression_forest y1 y2 x1 x2 x3, gen(mreg) ndep(2) ntrees(500) seed(42) replace

* Check that mreg_y1 and mreg_y2 exist and have non-missing values
forvalues j = 1/2 {
    quietly count if !missing(mreg_y`j')
    local n_pred = r(N)
    display as text "  mreg_y`j': " as result `n_pred' " / 500 non-missing"
    assert `n_pred' > 400
}

* Predictions should correlate with true values
corr mreg_y1 y1
local r1 = r(rho)
corr mreg_y2 y2
local r2 = r(rho)
display as text "  Correlation mreg_y1 vs y1: " as result %6.4f `r1'
display as text "  Correlation mreg_y2 vs y2: " as result %6.4f `r2'
assert `r1' > 0.5
assert `r2' > 0.5

* e() results should be stored
assert "`e(cmd)'" == "grf_multi_regression_forest"
assert "`e(forest_type)'" == "multi_regression"
assert e(N) == 500
assert e(n_outcomes) == 2

display as text "  PASSED"

* ---- Test 2: 3 outcomes ----
display as text ""
display as text "--- Test 2: 3 outcomes ---"

gen double y3 = x1 * x2 + rnormal()

grf_multi_regression_forest y1 y2 y3 x1 x2 x3, gen(mreg3) ndep(3) ntrees(500) seed(42) replace

* Check all 3 output variables
forvalues j = 1/3 {
    quietly count if !missing(mreg3_y`j')
    local n_pred = r(N)
    display as text "  mreg3_y`j': " as result `n_pred' " / 500 non-missing"
    assert `n_pred' > 400
}

assert e(n_outcomes) == 3

display as text "  PASSED"

* ---- Test 3: Options ----
display as text ""
display as text "--- Test 3: Options (nohonesty, minnodesize) ---"

grf_multi_regression_forest y1 y2 x1 x2 x3, gen(mreg_opts) ndep(2) ///
    ntrees(200) seed(42) nohonesty minnodesize(10) replace

forvalues j = 1/2 {
    quietly count if !missing(mreg_opts_y`j')
    assert r(N) > 400
}

assert e(honesty) == 0
assert e(min_node) == 10

display as text "  PASSED"

* ---- Test 4: if/in restrictions ----
display as text ""
display as text "--- Test 4: if/in restrictions ---"

grf_multi_regression_forest y1 y2 x1 x2 x3 if x1 > 0, ///
    gen(mreg_sub) ndep(2) ntrees(200) seed(42) replace

quietly count if !missing(mreg_sub_y1) & x1 > 0
local n_pred = r(N)
quietly count if x1 > 0
local n_eligible = r(N)
display as text "  Predictions: " as result `n_pred' " / " as result `n_eligible' " (if x1 > 0)"
assert `n_pred' > 0

* Predictions should only exist for x1 > 0
quietly count if !missing(mreg_sub_y1) & x1 <= 0
assert r(N) == 0
display as text "  No predictions for x1 <= 0: OK"

display as text "  PASSED"

* ---- Test 5: ndep < 2 should fail ----
display as text ""
display as text "--- Test 5: Error handling (ndep < 2) ---"

capture noisily grf_multi_regression_forest y1 x1 x2 x3, gen(mreg_fail) ndep(1)
assert _rc != 0
display as text "  ndep(1) correctly rejected"
display as text "  PASSED"

* ---- Summary ----
display as text ""
display as text "=============================================="
display as text " All multi-regression forest tests completed"
display as text "=============================================="
display as text ""
