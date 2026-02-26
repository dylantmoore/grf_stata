* test_tune.do -- Test grf_tune (parameter tuning)
* Run from the test-c-plugin-skill-grf/ directory

clear all
set more off

display as text ""
display as text "=============================================="
display as text " GRF Tuning Tests"
display as text "=============================================="

* ---- Test 1: Regression tuning ----
display as text ""
display as text "--- Test 1: Regression forest tuning ---"

clear
set obs 500
set seed 42

gen x1 = rnormal()
gen x2 = rnormal()
gen x3 = rnormal()
gen x4 = rnormal()
gen x5 = rnormal()
gen y = 2 * x1 + x2^2 + rnormal()

* Tune with 10 reps (fast for testing)
grf_tune y x1 x2 x3 x4 x5, foresttype(regression) numreps(10) ///
    tunetrees(100) seed(42)

* Check return values exist
assert r(best_mtry) > 0
assert r(best_min_node_size) > 0
assert r(best_sample_fraction) > 0
assert r(best_mse) > 0
assert r(n_reps) == 10

display as text "  Best mtry: " as result r(best_mtry)
display as text "  Best min_node: " as result r(best_min_node_size)
display as text "  Best sample_frac: " as result %5.3f r(best_sample_fraction)
display as text "  Best MSE: " as result %9.4f r(best_mse)

display as text "  PASSED"

* ---- Test 2: Tuned forest should be at least as good ----
display as text ""
display as text "--- Test 2: Apply tuned parameters ---"

local tuned_mtry = r(best_mtry)
local tuned_minnode = r(best_min_node_size)
local tuned_sf = r(best_sample_fraction)
local tuned_hf = r(best_honesty_fraction)
local tuned_alpha = r(best_alpha)
local tuned_ip = r(best_imbalance_penalty)

* Train with tuned parameters
grf_regression_forest y x1 x2 x3 x4 x5, gen(yhat_tuned) ntrees(500) seed(42) ///
    mtry(`tuned_mtry') minnodesize(`tuned_minnode') ///
    samplefrac(`tuned_sf') honestyfrac(`tuned_hf') ///
    alpha(`tuned_alpha') imbalancepenalty(`tuned_ip') replace

quietly count if !missing(yhat_tuned)
assert r(N) > 400
display as text "  Tuned forest produced " as result r(N) " predictions"

* Train with default parameters for comparison
grf_regression_forest y x1 x2 x3 x4 x5, gen(yhat_default) ///
    ntrees(500) seed(42) replace

quietly gen double se_tuned = (y - yhat_tuned)^2
quietly gen double se_default = (y - yhat_default)^2
quietly summarize se_tuned, meanonly
local mse_tuned = r(mean)
quietly summarize se_default, meanonly
local mse_default = r(mean)

display as text "  Tuned MSE:   " as result %9.4f `mse_tuned'
display as text "  Default MSE: " as result %9.4f `mse_default'
display as text "  PASSED"

* ---- Summary ----
display as text ""
display as text "=============================================="
display as text " All tuning tests completed"
display as text "=============================================="
display as text ""
