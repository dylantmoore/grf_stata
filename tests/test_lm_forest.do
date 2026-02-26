* test_lm_forest.do -- Test grf_lm_forest
* Run from the test-c-plugin-skill-grf/ directory

clear all
set more off

display as text ""
display as text "=============================================="
display as text " GRF LM Forest Tests"
display as text "=============================================="

* ---- Test 1: Basic functionality (single regressor) ----
display as text ""
display as text "--- Test 1: Single regressor (treatment effect) ---"

clear
set obs 500
set seed 12345

gen x1 = rnormal()
gen x2 = rnormal()
gen w = (runiform() > 0.5)
gen y = 2*x1 + w*(1 + x2) + rnormal()

grf_lm_forest y w, gen(h) xvars(x1 x2) ntrees(500) seed(42)

* Should create h_1 (coefficient on w)
confirm variable h_1

quietly count if !missing(h_1)
local n_pred = r(N)
display as text "  Predictions: " as result `n_pred' " / 500 non-missing"
assert `n_pred' > 400

summarize h_1
display as text "  Mean h_1(x) [treatment effect]: " as result %9.4f r(mean)
display as text "  SD h_1(x):                       " as result %9.4f r(sd)
* Treatment effect is 1 + x2, mean should be around 1
assert abs(r(mean) - 1) < 1.0
assert r(sd) > 0

assert "`e(cmd)'" == "grf_lm_forest"
assert "`e(forest_type)'" == "lm_forest"
assert e(N) == 500
assert e(n_regressors) == 1

display as text "  PASSED"

* ---- Test 2: Multiple regressors ----
display as text ""
display as text "--- Test 2: Two regressors ---"

gen z = rnormal()
replace y = 2*x1 + w*(1 + x2) + z*(0.5*x1) + rnormal()

grf_lm_forest y w z, gen(coef) xvars(x1 x2) ntrees(500) seed(42) replace

confirm variable coef_1
confirm variable coef_2

quietly count if !missing(coef_1)
local n1 = r(N)
quietly count if !missing(coef_2)
local n2 = r(N)
display as text "  coef_1 non-missing: " as result `n1'
display as text "  coef_2 non-missing: " as result `n2'
assert `n1' > 400
assert `n2' > 400

summarize coef_1
display as text "  Mean coef_1 (w): " as result %9.4f r(mean)
summarize coef_2
display as text "  Mean coef_2 (z): " as result %9.4f r(mean)

assert e(n_regressors) == 2

display as text "  PASSED"

* ---- Test 3: Options (nohonesty, nuisancetrees) ----
display as text ""
display as text "--- Test 3: Options ---"

grf_lm_forest y w, gen(h_opt) xvars(x1 x2) ntrees(200) seed(42) ///
    nohonesty nuisancetrees(200) replace

quietly count if !missing(h_opt_1)
assert r(N) > 400
assert e(honesty) == 0

display as text "  PASSED"

* ---- Test 4: if/in restrictions ----
display as text ""
display as text "--- Test 4: if/in restrictions ---"

grf_lm_forest y w if x1 > 0, gen(h_sub) xvars(x1 x2) ntrees(200) seed(42) ///
    replace

quietly count if !missing(h_sub_1) & x1 > 0
local n_pred = r(N)
assert `n_pred' > 0

quietly count if !missing(h_sub_1) & x1 <= 0
assert r(N) == 0

display as text "  PASSED"

* ---- Summary ----
display as text ""
display as text "=============================================="
display as text " All LM forest tests completed"
display as text "=============================================="
display as text ""
