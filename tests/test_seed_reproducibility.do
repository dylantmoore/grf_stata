* test_seed_reproducibility.do -- Tests for seed-based reproducibility
*
* Verifies that the same seed produces identical predictions and that
* different seeds produce different predictions.
*
* Run: do tests/test_seed_reproducibility.do

clear all
set more off

local errors = 0

* ============================================================
* Generate common test data
* ============================================================
clear
set obs 300
set seed 12345
gen x1 = rnormal()
gen x2 = rnormal()
gen x3 = rnormal()
gen w = rbinomial(1, 0.5)
gen y = 2*x1 + x1*w + 0.5*rnormal()

* ============================================================
* 1. Regression forest: same seed => identical predictions
* ============================================================
display as text ""
display as text "=== Seed Test 1: Regression forest same seed ==="

capture noisily {
    grf_regression_forest y x1 x2 x3, gen(pred_r1) ntrees(200) seed(42) replace
    grf_regression_forest y x1 x2 x3, gen(pred_r2) ntrees(200) seed(42) replace

    gen diff_r = abs(pred_r1 - pred_r2)
    summarize diff_r
    * Predictions must be identical (zero difference)
    assert r(max) == 0
}
if _rc {
    display as error "FAIL: regression forest same seed"
    local errors = `errors' + 1
}
else {
    display as result "PASS: regression forest same seed => identical predictions"
}

* ============================================================
* 2. Causal forest: same seed => identical CATEs
* ============================================================
display as text ""
display as text "=== Seed Test 2: Causal forest same seed ==="

capture noisily {
    grf_causal_forest y w x1 x2 x3, gen(cate_c1) ntrees(200) seed(42) replace
    grf_causal_forest y w x1 x2 x3, gen(cate_c2) ntrees(200) seed(42) replace

    gen diff_c = abs(cate_c1 - cate_c2)
    summarize diff_c
    assert r(max) == 0
}
if _rc {
    display as error "FAIL: causal forest same seed"
    local errors = `errors' + 1
}
else {
    display as result "PASS: causal forest same seed => identical CATEs"
}

* ============================================================
* 3. Different seeds => different predictions
* ============================================================
display as text ""
display as text "=== Seed Test 3: Different seeds => different predictions ==="

capture noisily {
    grf_regression_forest y x1 x2 x3, gen(pred_s1) ntrees(200) seed(42) replace
    grf_regression_forest y x1 x2 x3, gen(pred_s2) ntrees(200) seed(99) replace

    gen diff_s = abs(pred_s1 - pred_s2)
    summarize diff_s
    * Predictions should differ with different seeds
    assert r(sd) > 0
}
if _rc {
    display as error "FAIL: different seeds produce different predictions"
    local errors = `errors' + 1
}
else {
    display as result "PASS: different seeds => different predictions"
}

* ============================================================
* Summary
* ============================================================
display ""
if `errors' > 0 {
    display as error "FAILED: `errors' seed reproducibility test(s) did not pass"
    exit 1
}
else {
    display as result "ALL SEED REPRODUCIBILITY TESTS PASSED"
}
