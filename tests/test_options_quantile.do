* test_options_quantile.do -- Comprehensive option tests for grf_quantile_forest
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
    grf_quantile_forest y x1-x5, gen(qpred1) ntrees(100) seed(42) mtry(3)
    assert !missing(qpred1_q10) in 1
    assert !missing(qpred1_q50) in 1
    assert !missing(qpred1_q90) in 1
    assert e(mtry) == 3
    drop qpred1_q10 qpred1_q50 qpred1_q90
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
    grf_quantile_forest y x1-x5, gen(qpred2) ntrees(100) seed(42) minnodesize(10)
    assert !missing(qpred2_q10) in 1
    assert e(min_node) == 10
    drop qpred2_q10 qpred2_q50 qpred2_q90
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
    grf_quantile_forest y x1-x5, gen(qpred3) ntrees(100) seed(42) samplefrac(0.7)
    assert !missing(qpred3_q10) in 1
    assert reldif(e(sample_fraction), 0.7) < 1e-6
    drop qpred3_q10 qpred3_q50 qpred3_q90
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
    grf_quantile_forest y x1-x5, gen(qpred4) ntrees(100) seed(42) nohonesty
    assert !missing(qpred4_q10) in 1
    assert e(honesty) == 0
    drop qpred4_q10 qpred4_q50 qpred4_q90
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
    grf_quantile_forest y x1-x5, gen(qpred5) ntrees(100) seed(42) honestyfrac(0.7)
    assert !missing(qpred5_q10) in 1
    assert reldif(e(honesty_fraction), 0.7) < 1e-6
    drop qpred5_q10 qpred5_q50 qpred5_q90
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
    grf_quantile_forest y x1-x5, gen(qpred6) ntrees(100) seed(42) nohonestyprune
    assert !missing(qpred6_q10) in 1
    assert e(honesty_prune) == 0
    drop qpred6_q10 qpred6_q50 qpred6_q90
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
    grf_quantile_forest y x1-x5, gen(qpred7) ntrees(100) seed(42) alpha(0.1)
    assert !missing(qpred7_q10) in 1
    assert reldif(e(alpha), 0.1) < 1e-6
    drop qpred7_q10 qpred7_q50 qpred7_q90
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
    grf_quantile_forest y x1-x5, gen(qpred8) ntrees(100) seed(42) imbalancepenalty(0.5)
    assert !missing(qpred8_q10) in 1
    assert reldif(e(imbalance_penalty), 0.5) < 1e-6
    drop qpred8_q10 qpred8_q50 qpred8_q90
}
if _rc {
    display as error "FAIL: imbalancepenalty(0.5)"
    local errors = `errors' + 1
}
else {
    display as result "PASS: imbalancepenalty(0.5)"
}

* ---- Test 9: numthreads(2) ----
capture noisily {
    grf_quantile_forest y x1-x5, gen(qpred9) ntrees(100) seed(42) numthreads(2)
    assert !missing(qpred9_q10) in 1
    drop qpred9_q10 qpred9_q50 qpred9_q90
}
if _rc {
    display as error "FAIL: numthreads(2)"
    local errors = `errors' + 1
}
else {
    display as result "PASS: numthreads(2)"
}

* ---- Test 10: custom quantiles ----
capture noisily {
    grf_quantile_forest y x1-x5, gen(qpred10) ntrees(100) seed(42) quantiles(0.25 0.5 0.75)
    assert !missing(qpred10_q25) in 1
    assert !missing(qpred10_q50) in 1
    assert !missing(qpred10_q75) in 1
    assert e(n_quantiles) == 3
    drop qpred10_q25 qpred10_q50 qpred10_q75
}
if _rc {
    display as error "FAIL: custom quantiles(0.25 0.5 0.75)"
    local errors = `errors' + 1
}
else {
    display as result "PASS: custom quantiles(0.25 0.5 0.75)"
}

* ---- Test 11: single quantile ----
capture noisily {
    grf_quantile_forest y x1-x5, gen(qpred11) ntrees(100) seed(42) quantiles(0.5)
    assert !missing(qpred11_q50) in 1
    assert e(n_quantiles) == 1
    drop qpred11_q50
}
if _rc {
    display as error "FAIL: single quantile(0.5)"
    local errors = `errors' + 1
}
else {
    display as result "PASS: single quantile(0.5)"
}

* ---- Test 12: if restriction ----
capture noisily {
    grf_quantile_forest y x1-x5 if x1 > 0, gen(qpred12) ntrees(100) seed(42)
    quietly count if x1 > 0
    local n_sub = r(N)
    assert e(N) == `n_sub'
    drop qpred12_q10 qpred12_q50 qpred12_q90
}
if _rc {
    display as error "FAIL: if restriction"
    local errors = `errors' + 1
}
else {
    display as result "PASS: if restriction"
}

* ---- Test 13: replace ----
capture noisily {
    grf_quantile_forest y x1-x5, gen(qpred13) ntrees(100) seed(42)
    grf_quantile_forest y x1-x5, gen(qpred13) ntrees(100) seed(42) replace
    assert !missing(qpred13_q10) in 1
    drop qpred13_q10 qpred13_q50 qpred13_q90
}
if _rc {
    display as error "FAIL: replace"
    local errors = `errors' + 1
}
else {
    display as result "PASS: replace"
}

* ---- Test 14: All non-default options combined ----
capture noisily {
    grf_quantile_forest y x1-x5, gen(qpred14) ntrees(100) seed(123) ///
        mtry(3) minnodesize(10) samplefrac(0.7) honestyfrac(0.7)    ///
        nohonestyprune alpha(0.1) imbalancepenalty(0.5)              ///
        numthreads(2) quantiles(0.25 0.75) nohonesty
    assert !missing(qpred14_q25) in 1
    assert !missing(qpred14_q75) in 1
    assert e(mtry) == 3
    assert e(min_node) == 10
    assert reldif(e(sample_fraction), 0.7) < 1e-6
    assert e(honesty) == 0
    assert e(honesty_prune) == 0
    assert reldif(e(alpha), 0.1) < 1e-6
    assert reldif(e(imbalance_penalty), 0.5) < 1e-6
    assert e(n_quantiles) == 2
    drop qpred14_q25 qpred14_q75
}
if _rc {
    display as error "FAIL: all non-default options combined"
    local errors = `errors' + 1
}
else {
    display as result "PASS: all non-default options combined"
}

* ---- Test 15: Verify default e() values ----
capture noisily {
    grf_quantile_forest y x1-x5, gen(qpred15) ntrees(100) seed(42)
    assert e(honesty) == 1
    assert e(honesty_prune) == 1
    assert e(mtry) == 0
    assert e(min_node) == 5
    assert reldif(e(sample_fraction), 0.5) < 1e-6
    assert reldif(e(honesty_fraction), 0.5) < 1e-6
    assert reldif(e(alpha), 0.05) < 1e-6
    assert reldif(e(imbalance_penalty), 0.0) < 1e-6
    assert e(n_quantiles) == 3
    assert "`e(forest_type)'" == "quantile"
    drop qpred15_q10 qpred15_q50 qpred15_q90
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
    display as error "FAILED: `errors' errors in quantile forest option tests"
    exit 1
}
else {
    display as result "ALL QUANTILE FOREST OPTION TESTS PASSED"
}
