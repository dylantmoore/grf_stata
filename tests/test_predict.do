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

* ---- Test 4: Quantile forest predict on new data ----
display as text ""
display as text "--- Test 4: Quantile forest predict on new data ---"

clear
set obs 600
set seed 42

gen x1 = rnormal()
gen x2 = rnormal()
gen x3 = rnormal()
gen y_q = 2 * x1 + x2 + rnormal()

* Train quantile forest on first 400 obs
grf_quantile_forest y_q x1 x2 x3 in 1/400, gen(qpred) quantiles(0.25 0.5 0.75) ntrees(500) seed(42) replace

* Check OOB predictions exist
quietly count if !missing(qpred_q50)
display as text "  OOB median predictions: " as result r(N)
assert r(N) > 350

* Set y missing for test obs
replace y_q = . if _n > 400

* Predict on new data
grf_predict, gen(qpred_new) replace

* Check predictions exist for test obs (one var per quantile)
quietly count if !missing(qpred_new_q25) & _n > 400
local n_q25 = r(N)
display as text "  Test q25 predictions: " as result `n_q25'

quietly count if !missing(qpred_new_q50) & _n > 400
local n_q50 = r(N)
display as text "  Test q50 predictions: " as result `n_q50'

quietly count if !missing(qpred_new_q75) & _n > 400
local n_q75 = r(N)
display as text "  Test q75 predictions: " as result `n_q75'

assert `n_q25' > 150
assert `n_q50' > 150
assert `n_q75' > 150

* Check no predictions for training obs
quietly count if !missing(qpred_new_q50) & _n <= 400
display as text "  Train predictions (should be 0): " as result r(N)

display as text "  PASSED"

* ---- Test 5: Probability forest predict on new data ----
display as text ""
display as text "--- Test 5: Probability forest predict on new data ---"

clear
set obs 600
set seed 42

gen x1 = rnormal()
gen x2 = rnormal()
gen x3 = rnormal()
gen y_p = cond(x1 + x2 > 0, 1, 0)

* Train probability forest on first 400 obs
grf_probability_forest y_p x1 x2 x3 in 1/400, gen(ppred) ntrees(500) seed(42) replace

* Check OOB predictions exist
quietly count if !missing(ppred_c0)
display as text "  OOB class 0 predictions: " as result r(N)
assert r(N) > 350

* Set y missing for test obs
replace y_p = . if _n > 400

* Predict on new data
grf_predict, gen(ppred_new) replace

* Check predictions exist for test obs
quietly count if !missing(ppred_new_c0) & _n > 400
local n_c0 = r(N)
display as text "  Test class 0 predictions: " as result `n_c0'

quietly count if !missing(ppred_new_c1) & _n > 400
local n_c1 = r(N)
display as text "  Test class 1 predictions: " as result `n_c1'

assert `n_c0' > 150
assert `n_c1' > 150

* Check no predictions for training obs
quietly count if !missing(ppred_new_c0) & _n <= 400
display as text "  Train predictions (should be 0): " as result r(N)

display as text "  PASSED"

* ---- Test 6: Instrumental forest predict on new data ----
display as text ""
display as text "--- Test 6: Instrumental forest predict on new data ---"

clear
set obs 600
set seed 42

gen x1 = rnormal()
gen x2 = rnormal()
gen x3 = rnormal()
gen z_iv = (runiform() > 0.5)
gen w_iv = z_iv * (runiform() > 0.3)
gen y_iv = 2 * x1 + w_iv * x1 + rnormal()

* Train instrumental forest on first 400 obs (syntax: Y W Z X1..Xp)
grf_instrumental_forest y_iv w_iv z_iv x1 x2 x3 in 1/400, gen(ipred) ntrees(500) seed(42) replace

* Check OOB predictions exist
quietly count if !missing(ipred)
display as text "  OOB IV predictions: " as result r(N)
assert r(N) > 350

* Set y, w, z missing for test obs
replace y_iv = . if _n > 400
replace w_iv = . if _n > 400
replace z_iv = . if _n > 400

* Predict on new data
grf_predict, gen(ipred_new) replace

* Check predictions exist for test obs
quietly count if !missing(ipred_new) & _n > 400
local n_iv_pred = r(N)
display as text "  Test IV predictions: " as result `n_iv_pred'
assert `n_iv_pred' > 150

* Check no predictions for training obs
quietly count if !missing(ipred_new) & _n <= 400
display as text "  Train predictions (should be 0): " as result r(N)

display as text "  PASSED"

* ---- Test 7: Survival forest predict on new data ----
display as text ""
display as text "--- Test 7: Survival forest predict on new data ---"

clear
set obs 600
set seed 42

gen x1 = rnormal()
gen x2 = rnormal()
gen x3 = rnormal()
gen sv_time = -log(runiform()) * exp(-x1)
replace sv_time = max(sv_time, 0.001)
gen sv_status = (runiform() > 0.3)

* Train survival forest on first 400 obs (syntax: time status X1..Xp)
grf_survival_forest sv_time sv_status x1 x2 x3 in 1/400, gen(spred) ntrees(500) seed(42) replace

* Check OOB predictions exist (survival outputs stub_s1, stub_s2, ...)
quietly count if !missing(spred_s1)
display as text "  OOB survival predictions (s1): " as result r(N)
assert r(N) > 350

* Set time and status missing for test obs
replace sv_time = . if _n > 400
replace sv_status = . if _n > 400

* Predict on new data
grf_predict, gen(spred_new) replace

* Check predictions exist for test obs
quietly count if !missing(spred_new_s1) & _n > 400
local n_sv_pred = r(N)
display as text "  Test survival predictions (s1): " as result `n_sv_pred'
assert `n_sv_pred' > 150

quietly count if !missing(spred_new_s2) & _n > 400
display as text "  Test survival predictions (s2): " as result r(N)

* Check no predictions for training obs
quietly count if !missing(spred_new_s1) & _n <= 400
display as text "  Train predictions (should be 0): " as result r(N)

display as text "  PASSED"

* ---- Test 8: Multi-regression forest predict on new data ----
display as text ""
display as text "--- Test 8: Multi-regression forest predict on new data ---"

clear
set obs 600
set seed 42

gen x1 = rnormal()
gen x2 = rnormal()
gen x3 = rnormal()
gen mr_y1 = x1 + rnormal()
gen mr_y2 = x2 + rnormal()

* Train multi-regression forest on first 400 obs (syntax: y1 y2 x1..xp, ndep(2))
grf_multi_regression_forest mr_y1 mr_y2 x1 x2 x3 in 1/400, gen(mrpred) ndep(2) ntrees(500) seed(42) replace

* Check OOB predictions exist
quietly count if !missing(mrpred_y1)
display as text "  OOB y1 predictions: " as result r(N)
assert r(N) > 350

* Set y1, y2 missing for test obs
replace mr_y1 = . if _n > 400
replace mr_y2 = . if _n > 400

* Predict on new data
grf_predict, gen(mrpred_new) replace

* Check predictions exist for test obs
quietly count if !missing(mrpred_new_y1) & _n > 400
local n_mr_y1 = r(N)
display as text "  Test y1 predictions: " as result `n_mr_y1'

quietly count if !missing(mrpred_new_y2) & _n > 400
local n_mr_y2 = r(N)
display as text "  Test y2 predictions: " as result `n_mr_y2'

assert `n_mr_y1' > 150
assert `n_mr_y2' > 150

* Check no predictions for training obs
quietly count if !missing(mrpred_new_y1) & _n <= 400
display as text "  Train predictions (should be 0): " as result r(N)

display as text "  PASSED"

* ---- Test 9: Local linear regression forest predict on new data ----
display as text ""
display as text "--- Test 9: LL regression forest predict on new data ---"

clear
set obs 600
set seed 42

gen x1 = rnormal()
gen x2 = rnormal()
gen x3 = rnormal()
gen y_ll = sin(x1) + rnormal() * 0.5

* Train local linear regression forest on first 400 obs
grf_ll_regression_forest y_ll x1 x2 x3 in 1/400, gen(llpred) ntrees(500) seed(42) replace

* Check OOB predictions exist
quietly count if !missing(llpred)
display as text "  OOB LL predictions: " as result r(N)
assert r(N) > 350

* Set y missing for test obs
replace y_ll = . if _n > 400

* Predict on new data
grf_predict, gen(llpred_new) replace

* Check predictions exist for test obs
quietly count if !missing(llpred_new) & _n > 400
local n_ll_pred = r(N)
display as text "  Test LL predictions: " as result `n_ll_pred'
assert `n_ll_pred' > 150

* Check no predictions for training obs
quietly count if !missing(llpred_new) & _n <= 400
display as text "  Train predictions (should be 0): " as result r(N)

display as text "  PASSED"

* ---- Test 10: LM forest predict on new data ----
display as text ""
display as text "--- Test 10: LM forest predict on new data ---"

clear
set obs 600
set seed 42

gen x1 = rnormal()
gen x2 = rnormal()
gen x3 = rnormal()
gen lm_w1 = rnormal()
gen y_lm = x1 + lm_w1 * x2 + rnormal()

* Train LM forest on first 400 obs (syntax: y w1 [w2 ...], xvars(x1 x2 x3))
grf_lm_forest y_lm lm_w1 in 1/400, gen(lmpred) xvars(x1 x2 x3) ntrees(500) seed(42) replace

* Check OOB predictions exist (lm_forest outputs stub_1, stub_2, ...)
quietly count if !missing(lmpred_1)
display as text "  OOB LM predictions (coef 1): " as result r(N)
assert r(N) > 350

* Set y, w missing for test obs
replace y_lm = . if _n > 400
replace lm_w1 = . if _n > 400

* Predict on new data
grf_predict, gen(lmpred_new) replace

* Check predictions exist for test obs
quietly count if !missing(lmpred_new_1) & _n > 400
local n_lm_pred = r(N)
display as text "  Test LM predictions (coef 1): " as result `n_lm_pred'
assert `n_lm_pred' > 150

* Check no predictions for training obs
quietly count if !missing(lmpred_new_1) & _n <= 400
display as text "  Train predictions (should be 0): " as result r(N)

display as text "  PASSED"

* ---- Test 11: Multi-arm causal forest predict on new data ----
display as text ""
display as text "--- Test 11: Multi-arm causal forest predict on new data ---"

clear
set obs 600
set seed 42

gen x1 = rnormal()
gen x2 = rnormal()
gen x3 = rnormal()
gen mac_t1 = (runiform() > 0.5)
gen mac_t2 = (runiform() > 0.5) * (1 - mac_t1)
gen y_mac = x1 + mac_t1 * x2 + mac_t2 * x3 + rnormal()

* Train multi-arm causal forest on first 400 obs (syntax: y t1 t2 x1..xp, ntreat(2))
grf_multi_arm_causal_forest y_mac mac_t1 mac_t2 x1 x2 x3 in 1/400, gen(macpred) ntreat(2) ntrees(500) seed(42) replace

* Check OOB predictions exist (multi_causal outputs stub_t1, stub_t2)
quietly count if !missing(macpred_t1)
display as text "  OOB multi-arm t1 predictions: " as result r(N)
assert r(N) > 350

* Set y, t1, t2 missing for test obs
replace y_mac = . if _n > 400
replace mac_t1 = . if _n > 400
replace mac_t2 = . if _n > 400

* Predict on new data
grf_predict, gen(macpred_new) replace

* Check predictions exist for test obs
quietly count if !missing(macpred_new_t1) & _n > 400
local n_mac_t1 = r(N)
display as text "  Test multi-arm t1 predictions: " as result `n_mac_t1'

quietly count if !missing(macpred_new_t2) & _n > 400
local n_mac_t2 = r(N)
display as text "  Test multi-arm t2 predictions: " as result `n_mac_t2'

assert `n_mac_t1' > 150
assert `n_mac_t2' > 150

* Check no predictions for training obs
quietly count if !missing(macpred_new_t1) & _n <= 400
display as text "  Train predictions (should be 0): " as result r(N)

display as text "  PASSED"

* ---- Test 12: Causal survival forest predict on new data ----
display as text ""
display as text "--- Test 12: Causal survival forest predict on new data ---"

clear
set obs 600
set seed 42

gen x1 = rnormal()
gen x2 = rnormal()
gen x3 = rnormal()
gen cs_w = (runiform() > 0.5)
gen cs_time = -log(runiform()) * exp(-x1 - 0.5 * cs_w)
replace cs_time = max(cs_time, 0.001)
gen cs_status = (runiform() > 0.2)

* Train causal survival forest on first 400 obs (syntax: time status w x1..xp)
grf_causal_survival_forest cs_time cs_status cs_w x1 x2 x3 in 1/400, gen(cspred) ntrees(500) seed(42) replace

* Check OOB predictions exist
quietly count if !missing(cspred)
display as text "  OOB causal survival predictions: " as result r(N)
assert r(N) > 350

* Set time, status, w missing for test obs
replace cs_time = . if _n > 400
replace cs_status = . if _n > 400
replace cs_w = . if _n > 400

* Predict on new data
grf_predict, gen(cspred_new) replace

* Check predictions exist for test obs
quietly count if !missing(cspred_new) & _n > 400
local n_cs_pred = r(N)
display as text "  Test causal survival predictions: " as result `n_cs_pred'
assert `n_cs_pred' > 150

* Check no predictions for training obs
quietly count if !missing(cspred_new) & _n <= 400
display as text "  Train predictions (should be 0): " as result r(N)

display as text "  PASSED"

* ---- Test 13: Boosted regression predict (should error) ----
display as text ""
display as text "--- Test 13: Boosted regression predict error ---"

clear
set obs 600
set seed 42

gen x1 = rnormal()
gen x2 = rnormal()
gen x3 = rnormal()
gen y_br = 2 * x1 + x2 + rnormal()

* Train boosted regression forest on first 400 obs
grf_boosted_regression_forest y_br x1 x2 x3 in 1/400, gen(brpred) ntrees(500) seed(42) replace

* Set y missing for test obs
replace y_br = . if _n > 400

* Predict should fail with informative error
capture grf_predict, gen(brpred_new) replace
assert _rc != 0
display as text "  Boosted regression predict error: OK (rc=" as result _rc as text ")"

display as text "  PASSED"

* ---- Summary ----
display as text ""
display as text "=============================================="
display as text " All predict tests completed (Tests 1-13)"
display as text "=============================================="
display as text ""
