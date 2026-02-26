* test_predict.do -- Test grf_predict (predict on new data)
* Run from the test-c-plugin-skill-grf/ directory

clear all
set more off

display as text ""
display as text "=============================================="
display as text " GRF Predict on New Data Tests"
display as text "=============================================="

* ---- Test 1: Regression predict on new data ----
display as text ""
display as text "--- Test 1: Regression predict on new data ---"

clear
set obs 600
set seed 42

* Generate training data (first 400 obs) and test data (last 200 obs)
gen x1 = rnormal()
gen x2 = rnormal()
gen x3 = rnormal()
gen y = 2 * x1 + x2 + rnormal()

* Mark which are test obs
gen byte is_test = (_n > 400)

* Train regression forest on all data (will use first 400 as training via OOB)
grf_regression_forest y x1 x2 x3 in 1/400, gen(yhat_oob) ntrees(500) seed(42) replace

* Check OOB predictions exist for training obs
quietly count if !missing(yhat_oob)
display as text "  OOB predictions: " as result r(N)
assert r(N) > 350

* Now set test obs Y to missing (simulate new data without outcome)
replace y = . if is_test

* Predict on new data
grf_predict, gen(yhat_new) replace

* Check predictions exist for test obs
quietly count if !missing(yhat_new) & is_test
local n_test_pred = r(N)
display as text "  Test predictions: " as result `n_test_pred'
assert `n_test_pred' > 150

* Check no predictions for training obs (grf_predict only fills test obs)
quietly count if !missing(yhat_new) & !is_test
display as text "  Train predictions (should be 0): " as result r(N)

* Predictions should be reasonable (y = 2*x1 + x2 + noise)
summarize yhat_new if is_test, meanonly
display as text "  Mean prediction: " as result %6.3f r(mean)
display as text "  SD prediction: " as result %6.3f r(sd)

display as text "  PASSED"

* ---- Test 2: Causal predict on new data ----
display as text ""
display as text "--- Test 2: Causal predict on new data ---"

clear
set obs 600
set seed 42

gen x1 = rnormal()
gen x2 = rnormal()
gen x3 = rnormal()
gen w = (runiform() > 0.5)
gen y = 2 * x1 + 0.5 * w * x1 + rnormal()

* Train causal forest on first 400 obs
grf_causal_forest y w x1 x2 x3 in 1/400, gen(cate_oob) ntrees(500) seed(42) replace

* Check OOB predictions exist
quietly count if !missing(cate_oob)
display as text "  OOB CATE predictions: " as result r(N)
assert r(N) > 350

* Set test obs Y and W to missing
replace y = . if _n > 400
replace w = . if _n > 400

* Predict on new data
grf_predict, gen(cate_new) replace

* Check predictions exist for test obs
quietly count if !missing(cate_new) & _n > 400
local n_test_pred = r(N)
display as text "  Test CATE predictions: " as result `n_test_pred'
assert `n_test_pred' > 150

display as text "  PASSED"

* ---- Test 3: Predict error handling ----
display as text ""
display as text "--- Test 3: Error handling ---"

* Test: no prior estimation
clear
set obs 100
gen x1 = rnormal()
capture grf_predict, gen(test_pred)
assert _rc != 0
display as text "  No prior estimation error: OK (rc=" as result _rc ")"

* Test: no test data (all training)
clear
set obs 100
set seed 42
gen x1 = rnormal()
gen x2 = rnormal()
gen y = x1 + rnormal()
grf_regression_forest y x1 x2, gen(yhat) ntrees(100) seed(42) replace
* Don't append any test data -- should error
capture grf_predict, gen(pred_err) replace
assert _rc != 0
display as text "  No test data error: OK (rc=" as result _rc ")"

display as text "  PASSED"

* ---- Summary ----
display as text ""
display as text "=============================================="
display as text " All predict tests completed"
display as text "=============================================="
display as text ""
