* test_options_causal.do -- Comprehensive option tests for grf_causal_forest
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
gen w = (x1 + rnormal() > 0)

* ---- Test 1: mtry(3) ----
capture noisily {
    grf_causal_forest y w x1-x5, gen(tau1) ntrees(100) seed(42) mtry(3)
    assert !missing(tau1) in 1
    assert e(mtry) == 3
    drop tau1 _grf_yhat _grf_what
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
    grf_causal_forest y w x1-x5, gen(tau2) ntrees(100) seed(42) minnodesize(10)
    assert !missing(tau2) in 1
    assert e(min_node) == 10
    drop tau2 _grf_yhat _grf_what
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
    grf_causal_forest y w x1-x5, gen(tau3) ntrees(100) seed(42) samplefrac(0.7)
    assert !missing(tau3) in 1
    assert reldif(e(sample_fraction), 0.7) < 1e-6
    drop tau3 _grf_yhat _grf_what
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
    grf_causal_forest y w x1-x5, gen(tau4) ntrees(100) seed(42) nohonesty
    assert !missing(tau4) in 1
    assert e(honesty) == 0
    drop tau4 _grf_yhat _grf_what
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
    grf_causal_forest y w x1-x5, gen(tau5) ntrees(100) seed(42) honestyfrac(0.7)
    assert !missing(tau5) in 1
    assert reldif(e(honesty_fraction), 0.7) < 1e-6
    drop tau5 _grf_yhat _grf_what
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
    grf_causal_forest y w x1-x5, gen(tau6) ntrees(100) seed(42) nohonestyprune
    assert !missing(tau6) in 1
    assert e(honesty_prune) == 0
    drop tau6 _grf_yhat _grf_what
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
    grf_causal_forest y w x1-x5, gen(tau7) ntrees(100) seed(42) alpha(0.1)
    assert !missing(tau7) in 1
    assert reldif(e(alpha), 0.1) < 1e-6
    drop tau7 _grf_yhat _grf_what
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
    grf_causal_forest y w x1-x5, gen(tau8) ntrees(100) seed(42) imbalancepenalty(0.5)
    assert !missing(tau8) in 1
    assert reldif(e(imbalance_penalty), 0.5) < 1e-6
    drop tau8 _grf_yhat _grf_what
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
    grf_causal_forest y w x1-x5, gen(tau9) ntrees(100) seed(42) cigroupsize(2) samplefrac(0.4)
    assert !missing(tau9) in 1
    assert e(ci_group_size) == 2
    drop tau9 _grf_yhat _grf_what
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
    grf_causal_forest y w x1-x5, gen(tau10) ntrees(100) seed(42) numthreads(2)
    assert !missing(tau10) in 1
    drop tau10 _grf_yhat _grf_what
}
if _rc {
    display as error "FAIL: numthreads(2)"
    local errors = `errors' + 1
}
else {
    display as result "PASS: numthreads(2)"
}

* ---- Test 11: nostabilizesplits ----
capture noisily {
    grf_causal_forest y w x1-x5, gen(tau11) ntrees(100) seed(42) nostabilizesplits
    assert !missing(tau11) in 1
    assert e(stabilize) == 0
    drop tau11 _grf_yhat _grf_what
}
if _rc {
    display as error "FAIL: nostabilizesplits"
    local errors = `errors' + 1
}
else {
    display as result "PASS: nostabilizesplits"
}

* ---- Test 12: estimatevariance ----
capture noisily {
    grf_causal_forest y w x1-x5, gen(tau12) ntrees(100) seed(42) estimatevariance
    assert !missing(tau12) in 1
    assert !missing(tau12_var) in 1
    assert tau12_var[1] > 0
    assert "`e(variance_var)'" == "tau12_var"
    drop tau12 tau12_var _grf_yhat _grf_what
}
if _rc {
    display as error "FAIL: estimatevariance"
    local errors = `errors' + 1
}
else {
    display as result "PASS: estimatevariance"
}

* ---- Test 13: estimatevariance + vargenerate ----
capture noisily {
    grf_causal_forest y w x1-x5, gen(tau13) ntrees(100) seed(42) estimatevariance vargenerate(myvar13)
    assert !missing(tau13) in 1
    assert !missing(myvar13) in 1
    assert myvar13[1] > 0
    assert "`e(variance_var)'" == "myvar13"
    drop tau13 myvar13 _grf_yhat _grf_what
}
if _rc {
    display as error "FAIL: estimatevariance + vargenerate"
    local errors = `errors' + 1
}
else {
    display as result "PASS: estimatevariance + vargenerate"
}

* ---- Test 14: yhatgenerate + whatgenerate ----
capture noisily {
    grf_causal_forest y w x1-x5, gen(tau14) ntrees(100) seed(42) ///
        yhatgenerate(my_yhat) whatgenerate(my_what)
    assert !missing(tau14) in 1
    assert !missing(my_yhat) in 1
    assert !missing(my_what) in 1
    drop tau14 my_yhat my_what _grf_yhat _grf_what
}
if _rc {
    display as error "FAIL: yhatgenerate + whatgenerate"
    local errors = `errors' + 1
}
else {
    display as result "PASS: yhatgenerate + whatgenerate"
}

* ---- Test 15: nuisancetrees(100) ----
capture noisily {
    grf_causal_forest y w x1-x5, gen(tau15) ntrees(100) seed(42) nuisancetrees(100)
    assert !missing(tau15) in 1
    drop tau15 _grf_yhat _grf_what
}
if _rc {
    display as error "FAIL: nuisancetrees(100)"
    local errors = `errors' + 1
}
else {
    display as result "PASS: nuisancetrees(100)"
}

* ---- Test 16: if restriction ----
capture noisily {
    grf_causal_forest y w x1-x5 if x1 > 0, gen(tau16) ntrees(100) seed(42)
    quietly count if x1 > 0
    local n_sub = r(N)
    assert e(N) == `n_sub'
    drop tau16 _grf_yhat _grf_what
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
    grf_causal_forest y w x1-x5, gen(tau17) ntrees(100) seed(42)
    grf_causal_forest y w x1-x5, gen(tau17) ntrees(100) seed(42) replace
    assert !missing(tau17) in 1
    drop tau17 _grf_yhat _grf_what
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
    grf_causal_forest y w x1-x5, gen(tau18) ntrees(100) seed(123) ///
        mtry(3) minnodesize(10) samplefrac(0.4) honestyfrac(0.7)  ///
        nohonestyprune alpha(0.1) imbalancepenalty(0.5)            ///
        cigroupsize(2) numthreads(2) nostabilizesplits             ///
        estimatevariance vargenerate(var18)                        ///
        yhatgenerate(yh18) whatgenerate(wh18)                      ///
        nuisancetrees(100)
    assert !missing(tau18) in 1
    assert !missing(var18) in 1
    assert !missing(yh18) in 1
    assert !missing(wh18) in 1
    assert e(mtry) == 3
    assert e(min_node) == 10
    assert reldif(e(sample_fraction), 0.4) < 1e-6
    assert reldif(e(honesty_fraction), 0.7) < 1e-6
    assert e(honesty_prune) == 0
    assert reldif(e(alpha), 0.1) < 1e-6
    assert reldif(e(imbalance_penalty), 0.5) < 1e-6
    assert e(ci_group_size) == 2
    assert e(stabilize) == 0
    assert e(seed) == 123
    drop tau18 var18 yh18 wh18 _grf_yhat _grf_what
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
    grf_causal_forest y w x1-x5, gen(tau19) ntrees(100) seed(42)
    assert e(honesty) == 1
    assert e(honesty_prune) == 1
    assert e(stabilize) == 1
    assert e(mtry) == 0
    assert e(min_node) == 5
    assert reldif(e(sample_fraction), 0.5) < 1e-6
    assert reldif(e(honesty_fraction), 0.5) < 1e-6
    assert reldif(e(alpha), 0.05) < 1e-6
    assert reldif(e(imbalance_penalty), 0.0) < 1e-6
    assert e(ci_group_size) == 1
    assert "`e(forest_type)'" == "causal"
    assert "`e(depvar)'" == "y"
    assert "`e(treatvar)'" == "w"
    assert !missing(e(ate))
    assert !missing(e(ate_se))
    drop tau19 _grf_yhat _grf_what
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
    display as error "FAILED: `errors' errors in causal forest option tests"
    exit 1
}
else {
    display as result "ALL CAUSAL FOREST OPTION TESTS PASSED"
}
