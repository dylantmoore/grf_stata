* test_weights.do -- Tests for sample weights support
*
* Tests that the weights() option correctly passes observation-level
* sampling weights to the grf backend.
*
* Run: do tests/test_weights.do

clear all
set more off

local errors = 0

* ============================================================
* 1. Regression forest with weights succeeds
* ============================================================
display as text ""
display as text "=== Weight Test 1: Regression forest with weights ==="

capture noisily {
    clear
    set obs 500
    set seed 42
    gen x1 = rnormal()
    gen x2 = rnormal()
    gen y = 2*x1 + x2 + 0.5*rnormal()

    * Weights proportional to |x1|
    gen wt = abs(x1) + 0.1

    grf_regression_forest y x1 x2, gen(pred_wt) ntrees(500) seed(42) weights(wt)
    assert e(N) == 500
}
if _rc {
    display as error "FAIL: regression forest with weights"
    local errors = `errors' + 1
}
else {
    display as result "PASS: regression forest with weights"
}

* ============================================================
* 2. Weighted predictions differ from unweighted
* ============================================================
display as text ""
display as text "=== Weight Test 2: Weighted vs unweighted predictions differ ==="

capture noisily {
    grf_regression_forest y x1 x2, gen(pred_nowt) ntrees(500) seed(42) replace
    gen diff = pred_wt - pred_nowt
    summarize diff
    assert r(sd) > 0
}
if _rc {
    display as error "FAIL: weighted vs unweighted predictions differ"
    local errors = `errors' + 1
}
else {
    display as result "PASS: weighted vs unweighted predictions differ"
}

* ============================================================
* 3. Causal forest with weights -- ATE should differ
* ============================================================
display as text ""
display as text "=== Weight Test 3: Causal forest with weights ==="

capture noisily {
    clear
    set obs 500
    set seed 42
    gen x1 = rnormal()
    gen x2 = rnormal()
    gen w = rbinomial(1, 0.5)
    gen y = 2*x1 + x1*w + 0.5*rnormal()
    gen wt = abs(x1) + 0.1

    * With weights
    grf_causal_forest y w x1 x2, gen(cate_wt) ntrees(500) seed(42) weights(wt) replace
    grf_ate
    local ate_wt = r(ate)

    * Without weights
    grf_causal_forest y w x1 x2, gen(cate_nowt) ntrees(500) seed(42) replace
    grf_ate
    local ate_nowt = r(ate)

    * ATEs should differ
    assert `ate_wt' != `ate_nowt'
}
if _rc {
    display as error "FAIL: causal forest weighted vs unweighted ATE"
    local errors = `errors' + 1
}
else {
    display as result "PASS: weighted ATE=`ate_wt' vs unweighted ATE=`ate_nowt'"
}

* ============================================================
* 4. Error: negative weights
* ============================================================
display as text ""
display as text "=== Weight Test 4: Error on negative weights ==="

capture {
    clear
    set obs 100
    set seed 42
    gen x1 = rnormal()
    gen y = x1 + rnormal()
    gen wt = -1

    grf_regression_forest y x1, gen(pred_neg) ntrees(100) seed(42) weights(wt)
}
if _rc {
    display as result "PASS: error on negative weights"
}
else {
    display as error "FAIL: should have errored on negative weights"
    local errors = `errors' + 1
}

* ============================================================
* 5. Error: all-zero weights
* ============================================================
display as text ""
display as text "=== Weight Test 5: Error on all-zero weights ==="

capture {
    clear
    set obs 100
    set seed 42
    gen x1 = rnormal()
    gen y = x1 + rnormal()
    gen wt = 0

    grf_regression_forest y x1, gen(pred_zero) ntrees(100) seed(42) weights(wt)
}
if _rc {
    display as result "PASS: error on all-zero weights"
}
else {
    display as error "FAIL: should have errored on all-zero weights"
    local errors = `errors' + 1
}

* ============================================================
* Summary
* ============================================================
display ""
if `errors' > 0 {
    display as error "FAILED: `errors' weight test(s) did not pass"
    exit 1
}
else {
    display as result "ALL WEIGHT TESTS PASSED"
}
