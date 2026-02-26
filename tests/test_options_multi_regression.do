* test_options_multi_regression.do -- Comprehensive option tests for grf_multi_regression_forest
* Tests every option individually and in combination

clear all
set more off

local errors = 0

* ============================================================
* Setup data
* Syntax: grf_multi_regression_forest y1 y2 x1..xp, gen(...) ndep(2)
* ============================================================
clear
set obs 500
set seed 42
forvalues j = 1/5 {
    gen x`j' = rnormal()
}
gen y1 = 2*x1 + x2^2 + rnormal()
gen y2 = x1 - x3 + rnormal()

* ---- Test 1: mtry(3) ----
capture noisily {
    grf_multi_regression_forest y1 y2 x1-x5, gen(mrp1) ndep(2) ntrees(100) seed(42) mtry(3)
    assert !missing(mrp1_y1) in 1
    assert !missing(mrp1_y2) in 1
    assert e(mtry) == 3
    drop mrp1_y1 mrp1_y2
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
    grf_multi_regression_forest y1 y2 x1-x5, gen(mrp2) ndep(2) ntrees(100) seed(42) minnodesize(10)
    assert !missing(mrp2_y1) in 1
    assert e(min_node) == 10
    drop mrp2_y1 mrp2_y2
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
    grf_multi_regression_forest y1 y2 x1-x5, gen(mrp3) ndep(2) ntrees(100) seed(42) samplefrac(0.7)
    assert !missing(mrp3_y1) in 1
    assert reldif(e(sample_fraction), 0.7) < 1e-6
    drop mrp3_y1 mrp3_y2
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
    grf_multi_regression_forest y1 y2 x1-x5, gen(mrp4) ndep(2) ntrees(100) seed(42) nohonesty
    assert !missing(mrp4_y1) in 1
    assert e(honesty) == 0
    drop mrp4_y1 mrp4_y2
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
    grf_multi_regression_forest y1 y2 x1-x5, gen(mrp5) ndep(2) ntrees(100) seed(42) honestyfrac(0.7)
    assert !missing(mrp5_y1) in 1
    assert reldif(e(honesty_fraction), 0.7) < 1e-6
    drop mrp5_y1 mrp5_y2
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
    grf_multi_regression_forest y1 y2 x1-x5, gen(mrp6) ndep(2) ntrees(100) seed(42) nohonestyprune
    assert !missing(mrp6_y1) in 1
    assert e(honesty_prune) == 0
    drop mrp6_y1 mrp6_y2
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
    grf_multi_regression_forest y1 y2 x1-x5, gen(mrp7) ndep(2) ntrees(100) seed(42) alpha(0.1)
    assert !missing(mrp7_y1) in 1
    assert reldif(e(alpha), 0.1) < 1e-6
    drop mrp7_y1 mrp7_y2
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
    grf_multi_regression_forest y1 y2 x1-x5, gen(mrp8) ndep(2) ntrees(100) seed(42) imbalancepenalty(0.5)
    assert !missing(mrp8_y1) in 1
    assert reldif(e(imbalance_penalty), 0.5) < 1e-6
    drop mrp8_y1 mrp8_y2
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
    grf_multi_regression_forest y1 y2 x1-x5, gen(mrp9) ndep(2) ntrees(100) seed(42) cigroupsize(2) samplefrac(0.4)
    assert !missing(mrp9_y1) in 1
    assert e(ci_group_size) == 2
    drop mrp9_y1 mrp9_y2
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
    grf_multi_regression_forest y1 y2 x1-x5, gen(mrp10) ndep(2) ntrees(100) seed(42) numthreads(2)
    assert !missing(mrp10_y1) in 1
    drop mrp10_y1 mrp10_y2
}
if _rc {
    display as error "FAIL: numthreads(2)"
    local errors = `errors' + 1
}
else {
    display as result "PASS: numthreads(2)"
}

* ---- Test 11: estimatevariance (should warn and ignore for multi_regression) ----
capture noisily {
    grf_multi_regression_forest y1 y2 x1-x5, gen(mrp11) ndep(2) ntrees(100) seed(42) estimatevariance
    assert !missing(mrp11_y1) in 1
    assert !missing(mrp11_y2) in 1
    drop mrp11_y1 mrp11_y2
}
if _rc {
    display as error "FAIL: estimatevariance (should warn but succeed)"
    local errors = `errors' + 1
}
else {
    display as result "PASS: estimatevariance (warns but succeeds)"
}

* ---- Test 12: ndep(3) with 3 outcomes ----
capture noisily {
    gen y3 = x4 + rnormal()
    grf_multi_regression_forest y1 y2 y3 x1-x5, gen(mrp12) ndep(3) ntrees(100) seed(42)
    assert !missing(mrp12_y1) in 1
    assert !missing(mrp12_y2) in 1
    assert !missing(mrp12_y3) in 1
    assert e(n_outcomes) == 3
    drop mrp12_y1 mrp12_y2 mrp12_y3 y3
}
if _rc {
    display as error "FAIL: ndep(3)"
    local errors = `errors' + 1
}
else {
    display as result "PASS: ndep(3)"
}

* ---- Test 13: if restriction ----
capture noisily {
    grf_multi_regression_forest y1 y2 x1-x5 if x1 > 0, gen(mrp13) ndep(2) ntrees(100) seed(42)
    quietly count if x1 > 0
    local n_sub = r(N)
    assert e(N) == `n_sub'
    drop mrp13_y1 mrp13_y2
}
if _rc {
    display as error "FAIL: if restriction"
    local errors = `errors' + 1
}
else {
    display as result "PASS: if restriction"
}

* ---- Test 14: replace ----
capture noisily {
    grf_multi_regression_forest y1 y2 x1-x5, gen(mrp14) ndep(2) ntrees(100) seed(42)
    grf_multi_regression_forest y1 y2 x1-x5, gen(mrp14) ndep(2) ntrees(100) seed(42) replace
    assert !missing(mrp14_y1) in 1
    drop mrp14_y1 mrp14_y2
}
if _rc {
    display as error "FAIL: replace"
    local errors = `errors' + 1
}
else {
    display as result "PASS: replace"
}

* ---- Test 15: All non-default options combined ----
capture noisily {
    grf_multi_regression_forest y1 y2 x1-x5, gen(mrp15) ndep(2) ntrees(100) seed(123) ///
        mtry(3) minnodesize(10) samplefrac(0.4) honestyfrac(0.7)                       ///
        nohonestyprune alpha(0.1) imbalancepenalty(0.5) cigroupsize(2)                  ///
        numthreads(2) nohonesty
    assert !missing(mrp15_y1) in 1
    assert !missing(mrp15_y2) in 1
    assert e(mtry) == 3
    assert e(min_node) == 10
    assert reldif(e(sample_fraction), 0.4) < 1e-6
    assert e(honesty) == 0
    assert e(honesty_prune) == 0
    assert reldif(e(alpha), 0.1) < 1e-6
    assert reldif(e(imbalance_penalty), 0.5) < 1e-6
    assert e(ci_group_size) == 2
    assert e(n_outcomes) == 2
    drop mrp15_y1 mrp15_y2
}
if _rc {
    display as error "FAIL: all non-default options combined"
    local errors = `errors' + 1
}
else {
    display as result "PASS: all non-default options combined"
}

* ---- Test 16: Verify default e() values ----
capture noisily {
    grf_multi_regression_forest y1 y2 x1-x5, gen(mrp16) ndep(2) ntrees(100) seed(42)
    assert e(honesty) == 1
    assert e(honesty_prune) == 1
    assert e(mtry) == 0
    assert e(min_node) == 5
    assert reldif(e(sample_fraction), 0.5) < 1e-6
    assert reldif(e(alpha), 0.05) < 1e-6
    assert e(n_outcomes) == 2
    assert "`e(forest_type)'" == "multi_regression"
    drop mrp16_y1 mrp16_y2
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
    display as error "FAILED: `errors' errors in multi-regression forest option tests"
    exit 1
}
else {
    display as result "ALL MULTI-REGRESSION FOREST OPTION TESTS PASSED"
}
