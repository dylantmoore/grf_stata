* test_options_instrumental.do -- Comprehensive option tests for grf_instrumental_forest
* Tests every option individually and in combination

clear all
set more off

local errors = 0

* ============================================================
* Setup data
* Syntax: grf_instrumental_forest Y W Z X1..Xp, gen(...)
* ============================================================
clear
set obs 500
set seed 42
forvalues j = 1/5 {
    gen x`j' = rnormal()
}
* Instrument z -> treatment w -> outcome y
gen z = (x1 + rnormal() > 0)
gen w = 0.5*z + 0.3*x2 + (rnormal() > 0)
gen y = 2*w + x1 + x2^2 + rnormal()

* ---- Test 1: mtry(3) ----
capture noisily {
    grf_instrumental_forest y w z x1-x5, gen(iv1) ntrees(100) seed(42) mtry(3)
    assert !missing(iv1) in 1
    assert e(mtry) == 3
    drop iv1
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
    grf_instrumental_forest y w z x1-x5, gen(iv2) ntrees(100) seed(42) minnodesize(10)
    assert !missing(iv2) in 1
    assert e(min_node) == 10
    drop iv2
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
    grf_instrumental_forest y w z x1-x5, gen(iv3) ntrees(100) seed(42) samplefrac(0.7)
    assert !missing(iv3) in 1
    assert reldif(e(sample_fraction), 0.7) < 1e-6
    drop iv3
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
    grf_instrumental_forest y w z x1-x5, gen(iv4) ntrees(100) seed(42) nohonesty
    assert !missing(iv4) in 1
    assert e(honesty) == 0
    drop iv4
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
    grf_instrumental_forest y w z x1-x5, gen(iv5) ntrees(100) seed(42) honestyfrac(0.7)
    assert !missing(iv5) in 1
    assert reldif(e(honesty_fraction), 0.7) < 1e-6
    drop iv5
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
    grf_instrumental_forest y w z x1-x5, gen(iv6) ntrees(100) seed(42) nohonestyprune
    assert !missing(iv6) in 1
    assert e(honesty_prune) == 0
    drop iv6
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
    grf_instrumental_forest y w z x1-x5, gen(iv7) ntrees(100) seed(42) alpha(0.1)
    assert !missing(iv7) in 1
    assert reldif(e(alpha), 0.1) < 1e-6
    drop iv7
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
    grf_instrumental_forest y w z x1-x5, gen(iv8) ntrees(100) seed(42) imbalancepenalty(0.5)
    assert !missing(iv8) in 1
    assert reldif(e(imbalance_penalty), 0.5) < 1e-6
    drop iv8
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
    grf_instrumental_forest y w z x1-x5, gen(iv9) ntrees(100) seed(42) cigroupsize(2) samplefrac(0.4)
    assert !missing(iv9) in 1
    assert e(ci_group_size) == 2
    drop iv9
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
    grf_instrumental_forest y w z x1-x5, gen(iv10) ntrees(100) seed(42) numthreads(2)
    assert !missing(iv10) in 1
    drop iv10
}
if _rc {
    display as error "FAIL: numthreads(2)"
    local errors = `errors' + 1
}
else {
    display as result "PASS: numthreads(2)"
}

* ---- Test 11: reducedformweight(0.5) ----
capture noisily {
    grf_instrumental_forest y w z x1-x5, gen(iv11) ntrees(100) seed(42) reducedformweight(0.5)
    assert !missing(iv11) in 1
    assert reldif(e(reduced_form_wt), 0.5) < 1e-6
    drop iv11
}
if _rc {
    display as error "FAIL: reducedformweight(0.5)"
    local errors = `errors' + 1
}
else {
    display as result "PASS: reducedformweight(0.5)"
}

* ---- Test 12: nostabilizesplits (opt-out, default is ON for instrumental) ----
capture noisily {
    grf_instrumental_forest y w z x1-x5, gen(iv12) ntrees(100) seed(42) nostabilizesplits
    assert !missing(iv12) in 1
    assert e(stabilize_splits) == 0
    drop iv12
}
if _rc {
    display as error "FAIL: nostabilizesplits"
    local errors = `errors' + 1
}
else {
    display as result "PASS: nostabilizesplits"
}

* ---- Test 13: nuisancetrees(100) ----
capture noisily {
    grf_instrumental_forest y w z x1-x5, gen(iv13) ntrees(100) seed(42) nuisancetrees(100)
    assert !missing(iv13) in 1
    drop iv13
}
if _rc {
    display as error "FAIL: nuisancetrees(100)"
    local errors = `errors' + 1
}
else {
    display as result "PASS: nuisancetrees(100)"
}

* ---- Test 14: estimatevariance ----
capture noisily {
    grf_instrumental_forest y w z x1-x5, gen(iv14) ntrees(100) seed(42) estimatevariance
    assert !missing(iv14) in 1
    assert !missing(iv14_var) in 1
    assert iv14_var[1] > 0
    assert "`e(variance_var)'" == "iv14_var"
    drop iv14 iv14_var
}
if _rc {
    display as error "FAIL: estimatevariance"
    local errors = `errors' + 1
}
else {
    display as result "PASS: estimatevariance"
}

* ---- Test 15: estimatevariance + vargenerate ----
capture noisily {
    grf_instrumental_forest y w z x1-x5, gen(iv15) ntrees(100) seed(42) ///
        estimatevariance vargenerate(myivvar)
    assert !missing(iv15) in 1
    assert !missing(myivvar) in 1
    assert "`e(variance_var)'" == "myivvar"
    drop iv15 myivvar
}
if _rc {
    display as error "FAIL: estimatevariance + vargenerate"
    local errors = `errors' + 1
}
else {
    display as result "PASS: estimatevariance + vargenerate"
}

* ---- Test 16: if restriction ----
capture noisily {
    grf_instrumental_forest y w z x1-x5 if x1 > 0, gen(iv16) ntrees(100) seed(42)
    quietly count if x1 > 0
    local n_sub = r(N)
    assert e(N) == `n_sub'
    drop iv16
}
if _rc {
    display as error "FAIL: if restriction"
    local errors = `errors' + 1
}
else {
    display as result "PASS: if restriction"
}

* ---- Test 17: replace ----
capture noisily {
    grf_instrumental_forest y w z x1-x5, gen(iv17) ntrees(100) seed(42)
    grf_instrumental_forest y w z x1-x5, gen(iv17) ntrees(100) seed(42) replace
    assert !missing(iv17) in 1
    drop iv17
}
if _rc {
    display as error "FAIL: replace"
    local errors = `errors' + 1
}
else {
    display as result "PASS: replace"
}

* ---- Test 18: All non-default options combined ----
capture noisily {
    grf_instrumental_forest y w z x1-x5, gen(iv18) ntrees(100) seed(123) ///
        mtry(3) minnodesize(10) samplefrac(0.4) honestyfrac(0.7)         ///
        nohonestyprune alpha(0.1) imbalancepenalty(0.5)                   ///
        cigroupsize(2) numthreads(2) reducedformweight(0.5)               ///
        nostabilizesplits nuisancetrees(100) estimatevariance               ///
        vargenerate(var18)
    assert !missing(iv18) in 1
    assert !missing(var18) in 1
    assert e(mtry) == 3
    assert e(min_node) == 10
    assert reldif(e(sample_fraction), 0.4) < 1e-6
    assert reldif(e(honesty_fraction), 0.7) < 1e-6
    assert e(honesty_prune) == 0
    assert reldif(e(alpha), 0.1) < 1e-6
    assert reldif(e(imbalance_penalty), 0.5) < 1e-6
    assert e(ci_group_size) == 2
    assert reldif(e(reduced_form_wt), 0.5) < 1e-6
    assert e(stabilize_splits) == 0
    drop iv18 var18
}
if _rc {
    display as error "FAIL: all non-default options combined"
    local errors = `errors' + 1
}
else {
    display as result "PASS: all non-default options combined"
}

* ---- Test 19: Verify default e() values ----
capture noisily {
    grf_instrumental_forest y w z x1-x5, gen(iv19) ntrees(100) seed(42)
    assert e(honesty) == 1
    assert e(honesty_prune) == 1
    assert e(stabilize_splits) == 1
    assert reldif(e(reduced_form_wt), 0.0) < 1e-6
    assert e(mtry) == 0
    assert e(min_node) == 5
    assert reldif(e(sample_fraction), 0.5) < 1e-6
    assert reldif(e(alpha), 0.05) < 1e-6
    assert "`e(forest_type)'" == "instrumental"
    assert "`e(depvar)'" == "y"
    assert "`e(treatment)'" == "w"
    assert "`e(instrument)'" == "z"
    drop iv19
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
    display as error "FAILED: `errors' errors in instrumental forest option tests"
    exit 1
}
else {
    display as result "ALL INSTRUMENTAL FOREST OPTION TESTS PASSED"
}
