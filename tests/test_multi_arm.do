* test_multi_arm.do -- Test grf_multi_arm_causal_forest
* Run from the test-c-plugin-skill-grf/ directory

clear all
set more off

display as text ""
display as text "=============================================="
display as text " GRF Multi-Arm Causal Forest Tests"
display as text "=============================================="

* ---- Test 1: Basic (2 treatment arms) ----
display as text ""
display as text "--- Test 1: Basic (2 treatment arms) ---"

clear
set obs 500
set seed 42

* Generate covariates
gen x1 = rnormal()
gen x2 = rnormal()
gen x3 = rnormal()

* Generate binary treatment indicators
gen w1 = (x1 > 0)
gen w2 = (x2 > 0)

* Generate outcome with known treatment effects: 1.5 for w1, 0.8 for w2
gen y = 1.5 * w1 + 0.8 * w2 + x1 + 0.5 * rnormal()

* Run multi-arm causal forest
grf_multi_arm_causal_forest y w1 w2 x1 x2 x3, gen(macf) ntreat(2) ///
    ntrees(500) seed(42) replace

* Check that per-arm CATE variables exist
quietly count if !missing(macf_t1)
local n_t1 = r(N)
display as text "  Treatment arm 1 predictions: " as result `n_t1' " / 500"
assert `n_t1' > 400

quietly count if !missing(macf_t2)
local n_t2 = r(N)
display as text "  Treatment arm 2 predictions: " as result `n_t2' " / 500"
assert `n_t2' > 400

* Check predictions have reasonable variation
summarize macf_t1
assert r(sd) > 0.05
display as text "  CATE arm 1 SD: " as result %6.4f r(sd)

summarize macf_t2
assert r(sd) > 0.05
display as text "  CATE arm 2 SD: " as result %6.4f r(sd)

display as text "  PASSED"

* ---- Test 2: Variance estimation ----
display as text ""
display as text "--- Test 2: Variance estimation ---"

grf_multi_arm_causal_forest y w1 w2 x1 x2 x3, gen(macf2) ntreat(2) ///
    ntrees(500) seed(42) estimatevariance replace

* Check variance variables exist for each arm
quietly count if !missing(macf2_t1_var)
local n_var1 = r(N)
display as text "  Arm 1 variance estimates: " as result `n_var1' " / 500"
assert `n_var1' > 400

quietly count if !missing(macf2_t2_var)
local n_var2 = r(N)
display as text "  Arm 2 variance estimates: " as result `n_var2' " / 500"
assert `n_var2' > 400

* Variance should be non-negative
summarize macf2_t1_var
assert r(min) >= 0
display as text "  Arm 1 min variance: " as result %9.6f r(min) " (>= 0 required)"
display as text "  Arm 1 mean variance: " as result %9.6f r(mean)

summarize macf2_t2_var
assert r(min) >= 0
display as text "  Arm 2 min variance: " as result %9.6f r(min) " (>= 0 required)"
display as text "  Arm 2 mean variance: " as result %9.6f r(mean)

display as text "  PASSED"

* ---- Test 3: if/in restrictions ----
display as text ""
display as text "--- Test 3: if/in restrictions ---"

grf_multi_arm_causal_forest y w1 w2 x1 x2 x3 if x3 > 0, gen(macf3) ///
    ntreat(2) ntrees(200) seed(42) replace

* Check predictions exist only for eligible obs
quietly count if !missing(macf3_t1) & x3 > 0
local n_pred = r(N)
quietly count if x3 > 0
local n_eligible = r(N)
display as text "  Arm 1 predictions: " as result `n_pred' " / " as result `n_eligible' " (if x3 > 0)"
assert `n_pred' > 0

* Predictions should not exist for x3 <= 0
quietly count if !missing(macf3_t1) & x3 <= 0
assert r(N) == 0
display as text "  No arm 1 predictions for x3 <= 0: OK"

quietly count if !missing(macf3_t2) & x3 <= 0
assert r(N) == 0
display as text "  No arm 2 predictions for x3 <= 0: OK"

display as text "  PASSED"

* ---- Summary ----
display as text ""
display as text "=============================================="
display as text " All multi-arm causal forest tests completed"
display as text "=============================================="
display as text ""
