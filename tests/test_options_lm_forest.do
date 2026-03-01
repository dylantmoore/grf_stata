* test_options_lm_forest.do -- Comprehensive option tests for grf_lm_forest
* Tests every option individually and in combination

clear all
set more off

local errors = 0

* ============================================================
* Setup data
* Syntax: grf_lm_forest depvar W1 [W2 ...], xvars(X1 X2 ...) gen(...)
* Model: Y = c(X) + h_1(X)*W1 + ... + h_K(X)*Wk
* ============================================================
clear
set obs 500
set seed 42
forvalues j = 1/5 {
    gen x`j' = rnormal()
}
gen w1 = rnormal()
gen w2 = rnormal()
gen y = 2*x1 + 0.5*w1*x1 + 0.3*w2*x2 + rnormal()

* ---- Test 1: mtry(3) ----
capture noisily {
    grf_lm_forest y w1 w2, xvars(x1 x2 x3 x4 x5) gen(lm1) ntrees(100) seed(42) mtry(3)
    assert !missing(lm1_1) in 1
    assert !missing(lm1_2) in 1
    assert e(mtry) == 3
    drop lm1_1 lm1_2 _grf_lm_yhat _grf_lm_what1 _grf_lm_what2
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
    grf_lm_forest y w1 w2, xvars(x1 x2 x3 x4 x5) gen(lm2) ntrees(100) seed(42) minnodesize(10)
    assert !missing(lm2_1) in 1
    assert e(min_node) == 10
    drop lm2_1 lm2_2 _grf_lm_yhat _grf_lm_what1 _grf_lm_what2
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
    grf_lm_forest y w1 w2, xvars(x1 x2 x3 x4 x5) gen(lm3) ntrees(100) seed(42) samplefrac(0.7)
    assert !missing(lm3_1) in 1
    assert reldif(e(sample_fraction), 0.7) < 1e-6
    drop lm3_1 lm3_2 _grf_lm_yhat _grf_lm_what1 _grf_lm_what2
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
    grf_lm_forest y w1 w2, xvars(x1 x2 x3 x4 x5) gen(lm4) ntrees(100) seed(42) nohonesty
    assert !missing(lm4_1) in 1
    assert e(honesty) == 0
    drop lm4_1 lm4_2 _grf_lm_yhat _grf_lm_what1 _grf_lm_what2
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
    grf_lm_forest y w1 w2, xvars(x1 x2 x3 x4 x5) gen(lm5) ntrees(100) seed(42) honestyfrac(0.7)
    assert !missing(lm5_1) in 1
    assert reldif(e(honesty_fraction), 0.7) < 1e-6
    drop lm5_1 lm5_2 _grf_lm_yhat _grf_lm_what1 _grf_lm_what2
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
    grf_lm_forest y w1 w2, xvars(x1 x2 x3 x4 x5) gen(lm6) ntrees(100) seed(42) nohonestyprune
    assert !missing(lm6_1) in 1
    assert e(honesty_prune) == 0
    drop lm6_1 lm6_2 _grf_lm_yhat _grf_lm_what1 _grf_lm_what2
}
if _rc {
    display as error "FAIL: nohonestyprune"
    local errors = `errors' + 1
}
else {
    display as result "PASS: nohonestyprune"
}

* ---- Test 7: stabilizesplits (opt-in, default OFF) ----
capture noisily {
    grf_lm_forest y w1 w2, xvars(x1 x2 x3 x4 x5) gen(lm7) ntrees(100) seed(42) stabilizesplits
    assert !missing(lm7_1) in 1
    assert e(stabilize) == 1
    drop lm7_1 lm7_2 _grf_lm_yhat _grf_lm_what1 _grf_lm_what2
}
if _rc {
    display as error "FAIL: stabilizesplits"
    local errors = `errors' + 1
}
else {
    display as result "PASS: stabilizesplits"
}

* ---- Test 8: alpha(0.1) ----
capture noisily {
    grf_lm_forest y w1 w2, xvars(x1 x2 x3 x4 x5) gen(lm8) ntrees(100) seed(42) alpha(0.1)
    assert !missing(lm8_1) in 1
    assert reldif(e(alpha), 0.1) < 1e-6
    drop lm8_1 lm8_2 _grf_lm_yhat _grf_lm_what1 _grf_lm_what2
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
    grf_lm_forest y w1 w2, xvars(x1 x2 x3 x4 x5) gen(lm9) ntrees(100) seed(42) imbalancepenalty(0.5)
    assert !missing(lm9_1) in 1
    assert reldif(e(imbalance_penalty), 0.5) < 1e-6
    drop lm9_1 lm9_2 _grf_lm_yhat _grf_lm_what1 _grf_lm_what2
}
if _rc {
    display as error "FAIL: imbalancepenalty(0.5)"
    local errors = `errors' + 1
}
else {
    display as result "PASS: imbalancepenalty(0.5)"
}

* ---- Test 10: cigroupsize(2) ----
capture noisily {
    grf_lm_forest y w1 w2, xvars(x1 x2 x3 x4 x5) gen(lm10) ntrees(100) seed(42) cigroupsize(2) samplefrac(0.4)
    assert !missing(lm10_1) in 1
    assert e(ci_group_size) == 2
    drop lm10_1 lm10_2 _grf_lm_yhat _grf_lm_what1 _grf_lm_what2
}
if _rc {
    display as error "FAIL: cigroupsize(2)"
    local errors = `errors' + 1
}
else {
    display as result "PASS: cigroupsize(2)"
}

* ---- Test 11: numthreads(2) ----
capture noisily {
    grf_lm_forest y w1 w2, xvars(x1 x2 x3 x4 x5) gen(lm11) ntrees(100) seed(42) numthreads(2)
    assert !missing(lm11_1) in 1
    drop lm11_1 lm11_2 _grf_lm_yhat _grf_lm_what1 _grf_lm_what2
}
if _rc {
    display as error "FAIL: numthreads(2)"
    local errors = `errors' + 1
}
else {
    display as result "PASS: numthreads(2)"
}

* ---- Test 12: estimatevariance ----
capture noisily {
    grf_lm_forest y w1 w2, xvars(x1 x2 x3 x4 x5) gen(lm12) ntrees(100) seed(42) estimatevariance
    assert !missing(lm12_1) in 1
    assert !missing(lm12_2) in 1
    assert !missing(lm12_1_var) in 1
    assert !missing(lm12_2_var) in 1
    assert lm12_1_var[1] > 0
    assert lm12_2_var[1] > 0
    drop lm12_1 lm12_2 lm12_1_var lm12_2_var _grf_lm_yhat _grf_lm_what1 _grf_lm_what2
}
if _rc {
    display as error "FAIL: estimatevariance"
    local errors = `errors' + 1
}
else {
    display as result "PASS: estimatevariance"
}

* ---- Test 13: nuisancetrees(100) ----
capture noisily {
    grf_lm_forest y w1 w2, xvars(x1 x2 x3 x4 x5) gen(lm13) ntrees(100) seed(42) nuisancetrees(100)
    assert !missing(lm13_1) in 1
    drop lm13_1 lm13_2 _grf_lm_yhat _grf_lm_what1 _grf_lm_what2
}
if _rc {
    display as error "FAIL: nuisancetrees(100)"
    local errors = `errors' + 1
}
else {
    display as result "PASS: nuisancetrees(100)"
}

* ---- Test 14: Single regressor ----
capture noisily {
    grf_lm_forest y w1, xvars(x1 x2 x3 x4 x5) gen(lm14) ntrees(100) seed(42)
    assert !missing(lm14_1) in 1
    assert e(n_regressors) == 1
    drop lm14_1 _grf_lm_yhat _grf_lm_what1
}
if _rc {
    display as error "FAIL: single regressor"
    local errors = `errors' + 1
}
else {
    display as result "PASS: single regressor"
}

* ---- Test 15: if restriction ----
capture noisily {
    grf_lm_forest y w1 w2 if x1 > 0, xvars(x1 x2 x3 x4 x5) gen(lm15) ntrees(100) seed(42)
    quietly count if x1 > 0
    local n_sub = r(N)
    assert e(N) == `n_sub'
    drop lm15_1 lm15_2 _grf_lm_yhat _grf_lm_what1 _grf_lm_what2
}
if _rc {
    display as error "FAIL: if restriction"
    local errors = `errors' + 1
}
else {
    display as result "PASS: if restriction"
}

* ---- Test 16: replace ----
capture noisily {
    grf_lm_forest y w1 w2, xvars(x1 x2 x3 x4 x5) gen(lm16) ntrees(100) seed(42)
    grf_lm_forest y w1 w2, xvars(x1 x2 x3 x4 x5) gen(lm16) ntrees(100) seed(42) replace
    assert !missing(lm16_1) in 1
    drop lm16_1 lm16_2 _grf_lm_yhat _grf_lm_what1 _grf_lm_what2
}
if _rc {
    display as error "FAIL: replace"
    local errors = `errors' + 1
}
else {
    display as result "PASS: replace"
}

* ---- Test 17: All non-default options combined ----
capture noisily {
    grf_lm_forest y w1 w2, xvars(x1 x2 x3 x4 x5) gen(lm17) ntrees(100) seed(123) ///
        mtry(3) minnodesize(10) samplefrac(0.4) honestyfrac(0.7)                    ///
        nohonestyprune alpha(0.1) imbalancepenalty(0.5) cigroupsize(2)               ///
        numthreads(2) stabilizesplits estimatevariance nuisancetrees(100)
    assert !missing(lm17_1) in 1
    assert !missing(lm17_2) in 1
    assert !missing(lm17_1_var) in 1
    assert !missing(lm17_2_var) in 1
    assert e(mtry) == 3
    assert e(min_node) == 10
    assert reldif(e(sample_fraction), 0.4) < 1e-6
    assert reldif(e(honesty_fraction), 0.7) < 1e-6
    assert e(honesty_prune) == 0
    assert reldif(e(alpha), 0.1) < 1e-6
    assert reldif(e(imbalance_penalty), 0.5) < 1e-6
    assert e(ci_group_size) == 2
    assert e(stabilize) == 1
    assert e(n_regressors) == 2
    drop lm17_1 lm17_2 lm17_1_var lm17_2_var _grf_lm_yhat _grf_lm_what1 _grf_lm_what2
}
if _rc {
    display as error "FAIL: all non-default options combined"
    local errors = `errors' + 1
}
else {
    display as result "PASS: all non-default options combined"
}

* ---- Test 18: Verify default e() values ----
capture noisily {
    grf_lm_forest y w1 w2, xvars(x1 x2 x3 x4 x5) gen(lm18) ntrees(100) seed(42)
    assert e(honesty) == 1
    assert e(honesty_prune) == 1
    assert e(stabilize) == 0
    assert e(mtry) == 0
    assert e(min_node) == 5
    assert reldif(e(sample_fraction), 0.5) < 1e-6
    assert reldif(e(alpha), 0.05) < 1e-6
    assert e(n_regressors) == 2
    assert "`e(forest_type)'" == "lm_forest"
    assert "`e(depvar)'" == "y"
    assert !missing(e(coef_mean_1))
    assert !missing(e(coef_sd_1))
    assert !missing(e(coef_mean_2))
    assert !missing(e(coef_sd_2))
    drop lm18_1 lm18_2 _grf_lm_yhat _grf_lm_what1 _grf_lm_what2
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
    display as error "FAILED: `errors' errors in LM forest option tests"
    exit 1
}
else {
    display as result "ALL LM FOREST OPTION TESTS PASSED"
}
