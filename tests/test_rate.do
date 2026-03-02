* test_rate.do -- Test grf_rate (RATE)
* Run from the test-c-plugin-skill-grf/ directory

clear all
set more off

display as text ""
display as text "=============================================="
display as text " GRF RATE (Rank Average Treatment Effect) Tests"
display as text "=============================================="

* ---- Test 1: AUTOC with CATE as priorities ----
display as text ""
display as text "--- Test 1: AUTOC ---"

clear
set obs 500
set seed 42

* Generate data with heterogeneous treatment effects
gen x1 = rnormal()
gen x2 = rnormal()
gen x3 = rnormal()
gen w = (runiform() > 0.5)
* True CATE = x1 (positive effect for x1 > 0)
gen y = x1 + 0.5 * w * x1 + rnormal()

* Estimate causal forest
grf_causal_forest y w x1 x2 x3, gen(cate) ntrees(500) seed(42) replace

* Compute RATE using CATE as the prioritization score
grf_rate cate, target(AUTOC) bootstrap(50) seed(42)

* Check return values
assert r(estimate) != .
assert r(std_err) > 0
display as text "  AUTOC estimate: " as result %8.4f r(estimate)
display as text "  Std error:      " as result %8.4f r(std_err)
display as text "  z-statistic:    " as result %8.4f r(z_stat)

display as text "  PASSED"

* ---- Test 2: QINI ----
display as text ""
display as text "--- Test 2: QINI ---"

grf_rate cate, target(QINI) bootstrap(50) seed(42)

assert r(estimate) != .
assert r(std_err) > 0
display as text "  QINI estimate:  " as result %8.4f r(estimate)
display as text "  Std error:      " as result %8.4f r(std_err)

display as text "  PASSED"

* ---- Test 3: Custom priorities ----
display as text ""
display as text "--- Test 3: Custom priorities variable ---"

* Use x1 as custom prioritization score
grf_rate x1, target(AUTOC) bootstrap(50) ///
    catevar(cate) seed(42)

assert r(estimate) != .
display as text "  AUTOC with x1 priorities: " as result %8.4f r(estimate)

display as text "  PASSED"

* ---- Test 4: subset() matches if restriction ----
display as text ""
display as text "--- Test 4: subset() alias ---"

gen byte keep_subset = (x1 > 0)
grf_rate cate, target(AUTOC) bootstrap(50) seed(42) subset(keep_subset)
local est_subset = r(estimate)
local se_subset = r(std_err)
assert "`r(subset_var)'" == "keep_subset"

grf_rate cate if keep_subset, target(AUTOC) bootstrap(50) seed(42)
local est_if = r(estimate)
local se_if = r(std_err)

assert reldif(`est_subset', `est_if') < 1e-10
assert reldif(`se_subset', `se_if') < 1e-10

display as text "  subset estimate: " as result %8.4f `est_subset'
display as text "  if-estimate:     " as result %8.4f `est_if'
display as text "  PASSED"

* ---- Summary ----
display as text ""
display as text "=============================================="
display as text " All RATE tests completed"
display as text "=============================================="
display as text ""
