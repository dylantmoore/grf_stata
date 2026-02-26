* test_options_regression.do -- Comprehensive option tests for grf_regression_forest
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

* ---- Test 1: nohonestyprune ----
capture noisily {
    grf_regression_forest y x1-x5, gen(pred1) ntrees(100) seed(42) nohonestyprune
    assert !missing(pred1) in 1
    assert e(honesty_prune) == 0
    assert "`e(cmd)'" == "grf_regression_forest"
    drop pred1
}
if _rc {
    display as error "FAIL: nohonestyprune"
    local errors = `errors' + 1
}
else {
    display as result "PASS: nohonestyprune"
}

* ---- Test 2: cigroupsize(2) ----
capture noisily {
    grf_regression_forest y x1-x5, gen(pred2) ntrees(100) seed(42) cigroupsize(2) samplefrac(0.4)
    assert !missing(pred2) in 1
    assert e(ci_group_size) == 2
    drop pred2
}
if _rc {
    display as error "FAIL: cigroupsize(2)"
    local errors = `errors' + 1
}
else {
    display as result "PASS: cigroupsize(2)"
}

* ---- Test 3: numthreads(2) ----
capture noisily {
    grf_regression_forest y x1-x5, gen(pred3) ntrees(100) seed(42) numthreads(2)
    assert !missing(pred3) in 1
    assert "`e(cmd)'" == "grf_regression_forest"
    drop pred3
}
if _rc {
    display as error "FAIL: numthreads(2)"
    local errors = `errors' + 1
}
else {
    display as result "PASS: numthreads(2)"
}

* ---- Test 4: mtry(3) ----
capture noisily {
    grf_regression_forest y x1-x5, gen(pred4) ntrees(100) seed(42) mtry(3)
    assert !missing(pred4) in 1
    assert e(mtry) == 3
    drop pred4
}
if _rc {
    display as error "FAIL: mtry(3)"
    local errors = `errors' + 1
}
else {
    display as result "PASS: mtry(3)"
}

* ---- Test 5: minnodesize(10) ----
capture noisily {
    grf_regression_forest y x1-x5, gen(pred5) ntrees(100) seed(42) minnodesize(10)
    assert !missing(pred5) in 1
    assert e(min_node) == 10
    drop pred5
}
if _rc {
    display as error "FAIL: minnodesize(10)"
    local errors = `errors' + 1
}
else {
    display as result "PASS: minnodesize(10)"
}

* ---- Test 6: samplefrac(0.7) ----
capture noisily {
    grf_regression_forest y x1-x5, gen(pred6) ntrees(100) seed(42) samplefrac(0.7)
    assert !missing(pred6) in 1
    assert reldif(e(sample_fraction), 0.7) < 1e-6
    drop pred6
}
if _rc {
    display as error "FAIL: samplefrac(0.7)"
    local errors = `errors' + 1
}
else {
    display as result "PASS: samplefrac(0.7)"
}

* ---- Test 7: honestyfrac(0.7) ----
capture noisily {
    grf_regression_forest y x1-x5, gen(pred7) ntrees(100) seed(42) honestyfrac(0.7)
    assert !missing(pred7) in 1
    assert reldif(e(honesty_fraction), 0.7) < 1e-6
    drop pred7
}
if _rc {
    display as error "FAIL: honestyfrac(0.7)"
    local errors = `errors' + 1
}
else {
    display as result "PASS: honestyfrac(0.7)"
}

* ---- Test 8: alpha(0.1) ----
capture noisily {
    grf_regression_forest y x1-x5, gen(pred8) ntrees(100) seed(42) alpha(0.1)
    assert !missing(pred8) in 1
    assert reldif(e(alpha), 0.1) < 1e-6
    drop pred8
}
if _rc {
    display as error "FAIL: alpha(0.1)"
    local errors = `errors' + 1
}
else {
    display as result "PASS: alpha(0.1)"
}

* ---- Test 9: imbalancepenalty(0.5) ----
capture noisily {
    grf_regression_forest y x1-x5, gen(pred9) ntrees(100) seed(42) imbalancepenalty(0.5)
    assert !missing(pred9) in 1
    assert reldif(e(imbalance_penalty), 0.5) < 1e-6
    drop pred9
}
if _rc {
    display as error "FAIL: imbalancepenalty(0.5)"
    local errors = `errors' + 1
}
else {
    display as result "PASS: imbalancepenalty(0.5)"
}

* ---- Test 10: nohonesty ----
capture noisily {
    grf_regression_forest y x1-x5, gen(pred10) ntrees(100) seed(42) nohonesty
    assert !missing(pred10) in 1
    assert e(honesty) == 0
    drop pred10
}
if _rc {
    display as error "FAIL: nohonesty"
    local errors = `errors' + 1
}
else {
    display as result "PASS: nohonesty"
}

* ---- Test 11: estimatevariance ----
capture noisily {
    grf_regression_forest y x1-x5, gen(pred11) ntrees(100) seed(42) estimatevariance
    assert !missing(pred11) in 1
    assert !missing(pred11_var) in 1
    assert pred11_var[1] > 0
    * estimatevariance auto-bumps cigroupsize to 2
    assert e(ci_group_size) >= 2
    assert "`e(variance_var)'" == "pred11_var"
    drop pred11 pred11_var
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
    grf_regression_forest y x1-x5, gen(pred12) ntrees(100) seed(42) estimatevariance vargenerate(myvar12)
    assert !missing(pred12) in 1
    assert !missing(myvar12) in 1
    assert myvar12[1] > 0
    assert "`e(variance_var)'" == "myvar12"
    drop pred12 myvar12
}
if _rc {
    display as error "FAIL: estimatevariance + vargenerate"
    local errors = `errors' + 1
}
else {
    display as result "PASS: estimatevariance + vargenerate"
}

* ---- Test 13: if/in restriction ----
capture noisily {
    grf_regression_forest y x1-x5 if x1 > 0, gen(pred13) ntrees(100) seed(42)
    quietly count if x1 > 0
    local n_subset = r(N)
    assert e(N) == `n_subset'
    * Predictions should be non-missing for included obs, missing for excluded
    quietly count if !missing(pred13) & x1 > 0
    assert r(N) == `n_subset'
    drop pred13
}
if _rc {
    display as error "FAIL: if restriction"
    local errors = `errors' + 1
}
else {
    display as result "PASS: if restriction"
}

* ---- Test 14: replace option ----
capture noisily {
    grf_regression_forest y x1-x5, gen(pred14) ntrees(100) seed(42)
    assert !missing(pred14) in 1
    grf_regression_forest y x1-x5, gen(pred14) ntrees(100) seed(42) replace
    assert !missing(pred14) in 1
    drop pred14
}
if _rc {
    display as error "FAIL: replace"
    local errors = `errors' + 1
}
else {
    display as result "PASS: replace"
}

* ---- Test 15: All non-default options together ----
capture noisily {
    grf_regression_forest y x1-x5, gen(pred15) ntrees(100) seed(123) ///
        mtry(3) minnodesize(10) samplefrac(0.4) honestyfrac(0.7)    ///
        nohonestyprune alpha(0.1) imbalancepenalty(0.5)              ///
        cigroupsize(2) numthreads(2) estimatevariance                ///
        vargenerate(var15)
    assert !missing(pred15) in 1
    assert !missing(var15) in 1
    assert e(mtry) == 3
    assert e(min_node) == 10
    assert reldif(e(sample_fraction), 0.4) < 1e-6
    assert reldif(e(honesty_fraction), 0.7) < 1e-6
    assert e(honesty_prune) == 0
    assert reldif(e(alpha), 0.1) < 1e-6
    assert reldif(e(imbalance_penalty), 0.5) < 1e-6
    assert e(ci_group_size) == 2
    assert e(seed) == 123
    assert e(n_trees) == 100
    assert "`e(variance_var)'" == "var15"
    drop pred15 var15
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
    grf_regression_forest y x1-x5, gen(pred16) ntrees(100) seed(42)
    assert e(mtry) == 0
    assert e(min_node) == 5
    assert reldif(e(sample_fraction), 0.5) < 1e-6
    assert reldif(e(honesty_fraction), 0.5) < 1e-6
    assert e(honesty) == 1
    assert e(honesty_prune) == 1
    assert reldif(e(alpha), 0.05) < 1e-6
    assert reldif(e(imbalance_penalty), 0.0) < 1e-6
    assert e(ci_group_size) == 1
    assert "`e(forest_type)'" == "regression"
    assert "`e(depvar)'" == "y"
    assert "`e(predict_var)'" == "pred16"
    drop pred16
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
    display as error "FAILED: `errors' errors in regression forest option tests"
    exit 1
}
else {
    display as result "ALL REGRESSION FOREST OPTION TESTS PASSED"
}
