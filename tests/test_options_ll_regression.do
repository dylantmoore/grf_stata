* test_options_ll_regression.do -- Comprehensive option tests for grf_ll_regression_forest
* Tests every option individually and in combination

clear all
set more off

local errors = 0

* ============================================================
* Setup data
* ============================================================
clear
set obs 500
set seed 42
forvalues j = 1/5 {
    gen x`j' = rnormal()
}
gen y = 2*x1 + x2^2 + rnormal()

* ---- Test 1: mtry(3) ----
capture noisily {
    grf_ll_regression_forest y x1-x5, gen(ll1) ntrees(100) seed(42) mtry(3)
    assert !missing(ll1) in 1
    assert e(mtry) == 3
    drop ll1
}
if _rc {
    display as error "FAIL: mtry(3)"
    local errors = `errors' + 1
}
else {
    display as result "PASS: mtry(3)"
}

* ---- Test 2: minnodesize(10) ----
capture noisily {
    grf_ll_regression_forest y x1-x5, gen(ll2) ntrees(100) seed(42) minnodesize(10)
    assert !missing(ll2) in 1
    assert e(min_node) == 10
    drop ll2
}
if _rc {
    display as error "FAIL: minnodesize(10)"
    local errors = `errors' + 1
}
else {
    display as result "PASS: minnodesize(10)"
}

* ---- Test 3: samplefrac(0.7) ----
capture noisily {
    grf_ll_regression_forest y x1-x5, gen(ll3) ntrees(100) seed(42) samplefrac(0.7)
    assert !missing(ll3) in 1
    assert reldif(e(sample_fraction), 0.7) < 1e-6
    drop ll3
}
if _rc {
    display as error "FAIL: samplefrac(0.7)"
    local errors = `errors' + 1
}
else {
    display as result "PASS: samplefrac(0.7)"
}

* ---- Test 4: nohonesty ----
capture noisily {
    grf_ll_regression_forest y x1-x5, gen(ll4) ntrees(100) seed(42) nohonesty
    assert !missing(ll4) in 1
    assert e(honesty) == 0
    drop ll4
}
if _rc {
    display as error "FAIL: nohonesty"
    local errors = `errors' + 1
}
else {
    display as result "PASS: nohonesty"
}

* ---- Test 5: honestyfrac(0.7) ----
capture noisily {
    grf_ll_regression_forest y x1-x5, gen(ll5) ntrees(100) seed(42) honestyfrac(0.7)
    assert !missing(ll5) in 1
    assert reldif(e(honesty_fraction), 0.7) < 1e-6
    drop ll5
}
if _rc {
    display as error "FAIL: honestyfrac(0.7)"
    local errors = `errors' + 1
}
else {
    display as result "PASS: honestyfrac(0.7)"
}

* ---- Test 6: nohonestyprune ----
capture noisily {
    grf_ll_regression_forest y x1-x5, gen(ll6) ntrees(100) seed(42) nohonestyprune
    assert !missing(ll6) in 1
    assert e(honesty_prune) == 0
    drop ll6
}
if _rc {
    display as error "FAIL: nohonestyprune"
    local errors = `errors' + 1
}
else {
    display as result "PASS: nohonestyprune"
}

* ---- Test 7: alpha(0.1) ----
capture noisily {
    grf_ll_regression_forest y x1-x5, gen(ll7) ntrees(100) seed(42) alpha(0.1)
    assert !missing(ll7) in 1
    assert reldif(e(alpha), 0.1) < 1e-6
    drop ll7
}
if _rc {
    display as error "FAIL: alpha(0.1)"
    local errors = `errors' + 1
}
else {
    display as result "PASS: alpha(0.1)"
}

* ---- Test 8: imbalancepenalty(0.5) ----
capture noisily {
    grf_ll_regression_forest y x1-x5, gen(ll8) ntrees(100) seed(42) imbalancepenalty(0.5)
    assert !missing(ll8) in 1
    assert reldif(e(imbalance_penalty), 0.5) < 1e-6
    drop ll8
}
if _rc {
    display as error "FAIL: imbalancepenalty(0.5)"
    local errors = `errors' + 1
}
else {
    display as result "PASS: imbalancepenalty(0.5)"
}

* ---- Test 9: cigroupsize(2) ----
capture noisily {
    grf_ll_regression_forest y x1-x5, gen(ll9) ntrees(100) seed(42) cigroupsize(2) samplefrac(0.4)
    assert !missing(ll9) in 1
    assert e(ci_group_size) == 2
    drop ll9
}
if _rc {
    display as error "FAIL: cigroupsize(2)"
    local errors = `errors' + 1
}
else {
    display as result "PASS: cigroupsize(2)"
}

* ---- Test 10: numthreads(2) ----
capture noisily {
    grf_ll_regression_forest y x1-x5, gen(ll10) ntrees(100) seed(42) numthreads(2)
    assert !missing(ll10) in 1
    drop ll10
}
if _rc {
    display as error "FAIL: numthreads(2)"
    local errors = `errors' + 1
}
else {
    display as result "PASS: numthreads(2)"
}

* ---- Test 11: estimatevariance ----
capture noisily {
    grf_ll_regression_forest y x1-x5, gen(ll11) ntrees(100) seed(42) estimatevariance
    assert !missing(ll11) in 1
    assert !missing(ll11_var) in 1
    assert ll11_var[1] > 0
    assert "`e(variance_var)'" == "ll11_var"
    drop ll11 ll11_var
}
if _rc {
    display as error "FAIL: estimatevariance"
    local errors = `errors' + 1
}
else {
    display as result "PASS: estimatevariance"
}

* ---- Test 12: estimatevariance + vargenerate ----
capture noisily {
    grf_ll_regression_forest y x1-x5, gen(ll12) ntrees(100) seed(42) estimatevariance vargenerate(llvar12)
    assert !missing(ll12) in 1
    assert !missing(llvar12) in 1
    assert "`e(variance_var)'" == "llvar12"
    drop ll12 llvar12
}
if _rc {
    display as error "FAIL: estimatevariance + vargenerate"
    local errors = `errors' + 1
}
else {
    display as result "PASS: estimatevariance + vargenerate"
}

* ---- Test 13: llenable ----
capture noisily {
    grf_ll_regression_forest y x1-x5, gen(ll13) ntrees(100) seed(42) llenable
    assert !missing(ll13) in 1
    assert e(enable_ll_split) == 1
    drop ll13
}
if _rc {
    display as error "FAIL: llenable"
    local errors = `errors' + 1
}
else {
    display as result "PASS: llenable"
}

* ---- Test 14: lllambda(0.5) ----
capture noisily {
    grf_ll_regression_forest y x1-x5, gen(ll14) ntrees(100) seed(42) lllambda(0.5)
    assert !missing(ll14) in 1
    assert reldif(e(ll_lambda), 0.5) < 1e-6
    drop ll14
}
if _rc {
    display as error "FAIL: lllambda(0.5)"
    local errors = `errors' + 1
}
else {
    display as result "PASS: lllambda(0.5)"
}

* ---- Test 15: llweightpenalty ----
capture noisily {
    grf_ll_regression_forest y x1-x5, gen(ll15) ntrees(100) seed(42) llweightpenalty
    assert !missing(ll15) in 1
    assert e(ll_weight_penalty) == 1
    drop ll15
}
if _rc {
    display as error "FAIL: llweightpenalty"
    local errors = `errors' + 1
}
else {
    display as result "PASS: llweightpenalty"
}

* ---- Test 16: llcutoff(10) ----
capture noisily {
    grf_ll_regression_forest y x1-x5, gen(ll16) ntrees(100) seed(42) llcutoff(10)
    assert !missing(ll16) in 1
    assert e(ll_split_cutoff) == 10
    drop ll16
}
if _rc {
    display as error "FAIL: llcutoff(10)"
    local errors = `errors' + 1
}
else {
    display as result "PASS: llcutoff(10)"
}

* ---- Test 17: if restriction ----
capture noisily {
    grf_ll_regression_forest y x1-x5 if x1 > 0, gen(ll17) ntrees(100) seed(42)
    quietly count if x1 > 0
    local n_sub = r(N)
    assert e(N) == `n_sub'
    drop ll17
}
if _rc {
    display as error "FAIL: if restriction"
    local errors = `errors' + 1
}
else {
    display as result "PASS: if restriction"
}

* ---- Test 18: replace ----
capture noisily {
    grf_ll_regression_forest y x1-x5, gen(ll18) ntrees(100) seed(42)
    grf_ll_regression_forest y x1-x5, gen(ll18) ntrees(100) seed(42) replace
    assert !missing(ll18) in 1
    drop ll18
}
if _rc {
    display as error "FAIL: replace"
    local errors = `errors' + 1
}
else {
    display as result "PASS: replace"
}

* ---- Test 19: All non-default options combined ----
capture noisily {
    grf_ll_regression_forest y x1-x5, gen(ll19) ntrees(100) seed(123)  ///
        mtry(3) minnodesize(10) samplefrac(0.4) honestyfrac(0.7)      ///
        nohonestyprune alpha(0.1) imbalancepenalty(0.5) cigroupsize(2) ///
        numthreads(2) estimatevariance vargenerate(llvar19)            ///
        llenable lllambda(0.5) llweightpenalty llcutoff(10)
    assert !missing(ll19) in 1
    assert !missing(llvar19) in 1
    assert e(mtry) == 3
    assert e(min_node) == 10
    assert reldif(e(sample_fraction), 0.4) < 1e-6
    assert reldif(e(honesty_fraction), 0.7) < 1e-6
    assert e(honesty_prune) == 0
    assert reldif(e(alpha), 0.1) < 1e-6
    assert reldif(e(imbalance_penalty), 0.5) < 1e-6
    assert e(ci_group_size) == 2
    assert e(enable_ll_split) == 1
    assert reldif(e(ll_lambda), 0.5) < 1e-6
    assert e(ll_weight_penalty) == 1
    assert e(ll_split_cutoff) == 10
    drop ll19 llvar19
}
if _rc {
    display as error "FAIL: all non-default options combined"
    local errors = `errors' + 1
}
else {
    display as result "PASS: all non-default options combined"
}

* ---- Test 20: Verify default e() values ----
capture noisily {
    grf_ll_regression_forest y x1-x5, gen(ll20) ntrees(100) seed(42)
    assert e(honesty) == 1
    assert e(honesty_prune) == 1
    assert e(mtry) == 0
    assert e(min_node) == 5
    assert reldif(e(sample_fraction), 0.5) < 1e-6
    assert reldif(e(alpha), 0.05) < 1e-6
    assert e(enable_ll_split) == 0
    assert reldif(e(ll_lambda), 0.1) < 1e-6
    assert e(ll_weight_penalty) == 0
    assert e(ll_split_cutoff) == 0
    assert "`e(forest_type)'" == "ll_regression"
    drop ll20
}
if _rc {
    display as error "FAIL: default e() values"
    local errors = `errors' + 1
}
else {
    display as result "PASS: default e() values"
}

* ============================================================
* Summary
* ============================================================
display ""
if `errors' > 0 {
    display as error "FAILED: `errors' errors in LL regression forest option tests"
    exit 1
}
else {
    display as result "ALL LL REGRESSION FOREST OPTION TESTS PASSED"
}
