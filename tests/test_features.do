* test_features.do -- Test common features across all forest types
* Run from the test-c-plugin-skill-grf/ directory

clear all
set more off

display as text ""
display as text "=============================================="
display as text " GRF Common Feature Tests"
display as text "=============================================="

* ---- Setup: Generate shared dataset ----
clear
set obs 500
set seed 42

gen x1 = rnormal()
gen x2 = rnormal()
gen x3 = rnormal()
gen y = 2 * x1 + x2^2 + 0.5 * rnormal()

* ---- Test 1: replace option ----
display as text ""
display as text "--- Test 1: replace option ---"

* Run regression forest once
grf_regression_forest y x1 x2 x3, gen(rpred) ntrees(200) seed(42) replace

quietly count if !missing(rpred)
local n_first = r(N)
display as text "  First run predictions: " as result `n_first'
assert `n_first' > 400

* Run again with replace -- should overwrite without error
grf_regression_forest y x1 x2 x3, gen(rpred) ntrees(200) seed(99) replace

quietly count if !missing(rpred)
local n_second = r(N)
display as text "  Second run predictions: " as result `n_second'
assert `n_second' > 400

display as text "  PASSED"

* ---- Test 2: Missing data handling ----
display as text ""
display as text "--- Test 2: Missing data handling ---"

* Introduce missing values in y
gen y_miss = y
replace y_miss = . in 1/50

* Count non-missing observations
quietly count if !missing(y_miss)
local n_nonmiss = r(N)
display as text "  Non-missing obs: " as result `n_nonmiss'

* Run forest on data with missing y
grf_regression_forest y_miss x1 x2 x3, gen(pred_miss) ntrees(200) seed(42) replace

* Predictions should skip missing y observations
quietly count if !missing(pred_miss)
local n_pred = r(N)
display as text "  Predictions generated: " as result `n_pred'

* Should have predictions only for non-missing y
assert `n_pred' <= `n_nonmiss'
assert `n_pred' > 0
display as text "  Predictions <= non-missing obs: OK"

* Check that obs with missing y have missing predictions
quietly count if missing(y_miss) & !missing(pred_miss)
local n_leaked = r(N)
display as text "  Predictions on missing-y obs: " as result `n_leaked' " (expect 0)"
assert `n_leaked' == 0

display as text "  PASSED"

* ---- Test 3: Small sample ----
display as text ""
display as text "--- Test 3: Small sample (N=50) ---"

clear
set obs 50
set seed 42

gen x1 = rnormal()
gen x2 = rnormal()
gen y = x1 + 0.5 * rnormal()

* Should work even with only 50 observations
grf_regression_forest y x1 x2, gen(pred_small) ntrees(100) seed(42) replace

quietly count if !missing(pred_small)
local n_small = r(N)
display as text "  Predictions: " as result `n_small' " / 50"
assert `n_small' > 0

display as text "  PASSED"

* ---- Test 4: Single predictor ----
display as text ""
display as text "--- Test 4: Single predictor ---"

clear
set obs 300
set seed 42

gen x1 = rnormal()
gen y = 3 * x1 + rnormal()

* Run with just one predictor variable
grf_regression_forest y x1, gen(pred_single) ntrees(200) seed(42) replace

quietly count if !missing(pred_single)
local n_single = r(N)
display as text "  Predictions: " as result `n_single' " / 300"
assert `n_single' > 200

* Should still capture the relationship
corr y pred_single
local r_single = r(rho)
display as text "  Correlation(Y, pred): " as result %6.4f `r_single'
assert `r_single' > 0.3

display as text "  PASSED"

* ---- Test 5: In range ----
display as text ""
display as text "--- Test 5: In range (in 1/200) ---"

clear
set obs 500
set seed 42

gen x1 = rnormal()
gen x2 = rnormal()
gen x3 = rnormal()
gen y = 2 * x1 + x2^2 + 0.5 * rnormal()

* Run with in restriction
grf_regression_forest y x1 x2 x3 in 1/200, gen(pred_in) ntrees(200) seed(42) replace

* Predictions should exist only for obs 1-200
quietly count if !missing(pred_in) & _n <= 200
local n_in = r(N)
display as text "  Predictions in 1/200: " as result `n_in'
assert `n_in' > 0

quietly count if !missing(pred_in) & _n > 200
local n_out = r(N)
display as text "  Predictions outside range: " as result `n_out' " (expect 0)"
assert `n_out' == 0

display as text "  PASSED"

* ---- Summary ----
display as text ""
display as text "=============================================="
display as text " All common feature tests completed"
display as text "=============================================="
display as text ""
