* test_options_probability.do -- Comprehensive option tests for grf_probability_forest
* Tests every option individually and in combination

clear all
set more off

local errors = 0

* ============================================================
* Setup data (binary classification)
* ============================================================
clear
set obs 500
set seed 42
forvalues j = 1/5 {
    gen x`j' = rnormal()
}
gen yclass = (2*x1 + x2 + rnormal() > 0)

* ---- Test 1: mtry(3) ----
capture noisily {
    grf_probability_forest yclass x1-x5, gen(pp1) ntrees(100) seed(42) mtry(3)
    assert !missing(pp1_c0) in 1
    assert !missing(pp1_c1) in 1
    assert e(mtry) == 3
    drop pp1_c0 pp1_c1
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
    grf_probability_forest yclass x1-x5, gen(pp2) ntrees(100) seed(42) minnodesize(10)
    assert !missing(pp2_c0) in 1
    assert e(min_node) == 10
    drop pp2_c0 pp2_c1
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
    grf_probability_forest yclass x1-x5, gen(pp3) ntrees(100) seed(42) samplefrac(0.7)
    assert !missing(pp3_c0) in 1
    assert reldif(e(sample_fraction), 0.7) < 1e-6
    drop pp3_c0 pp3_c1
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
    grf_probability_forest yclass x1-x5, gen(pp4) ntrees(100) seed(42) nohonesty
    assert !missing(pp4_c0) in 1
    assert e(honesty) == 0
    drop pp4_c0 pp4_c1
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
    grf_probability_forest yclass x1-x5, gen(pp5) ntrees(100) seed(42) honestyfrac(0.7)
    assert !missing(pp5_c0) in 1
    assert reldif(e(honesty_fraction), 0.7) < 1e-6
    drop pp5_c0 pp5_c1
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
    grf_probability_forest yclass x1-x5, gen(pp6) ntrees(100) seed(42) nohonestyprune
    assert !missing(pp6_c0) in 1
    assert e(honesty_prune) == 0
    drop pp6_c0 pp6_c1
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
    grf_probability_forest yclass x1-x5, gen(pp7) ntrees(100) seed(42) alpha(0.1)
    assert !missing(pp7_c0) in 1
    assert reldif(e(alpha), 0.1) < 1e-6
    drop pp7_c0 pp7_c1
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
    grf_probability_forest yclass x1-x5, gen(pp8) ntrees(100) seed(42) imbalancepenalty(0.5)
    assert !missing(pp8_c0) in 1
    assert reldif(e(imbalance_penalty), 0.5) < 1e-6
    drop pp8_c0 pp8_c1
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
    grf_probability_forest yclass x1-x5, gen(pp9) ntrees(100) seed(42) numthreads(2)
    assert !missing(pp9_c0) in 1
    drop pp9_c0 pp9_c1
}
if _rc {
    display as error "FAIL: numthreads(2)"
    local errors = `errors' + 1
}
else {
    display as result "PASS: numthreads(2)"
}

* ---- Test 10: nclasses(2) explicit ----
capture noisily {
    grf_probability_forest yclass x1-x5, gen(pp10) ntrees(100) seed(42) nclasses(2)
    assert !missing(pp10_c0) in 1
    assert !missing(pp10_c1) in 1
    assert e(n_classes) == 2
    drop pp10_c0 pp10_c1
}
if _rc {
    display as error "FAIL: nclasses(2)"
    local errors = `errors' + 1
}
else {
    display as result "PASS: nclasses(2)"
}

* ---- Test 11: nclasses(3) with only 2 actual classes ----
* This should create 3 output columns even though data only has 0/1
capture noisily {
    grf_probability_forest yclass x1-x5, gen(pp11) ntrees(100) seed(42) nclasses(3)
    assert !missing(pp11_c0) in 1
    assert !missing(pp11_c1) in 1
    assert !missing(pp11_c2) in 1
    assert e(n_classes) == 3
    drop pp11_c0 pp11_c1 pp11_c2
}
if _rc {
    display as error "FAIL: nclasses(3) with binary data"
    local errors = `errors' + 1
}
else {
    display as result "PASS: nclasses(3) with binary data"
}

* ---- Test 12: 3-class data with auto-detection ----
capture noisily {
    gen yclass3 = floor(3*uniform())
    grf_probability_forest yclass3 x1-x5, gen(pp12) ntrees(100) seed(42)
    assert !missing(pp12_c0) in 1
    assert !missing(pp12_c1) in 1
    assert !missing(pp12_c2) in 1
    assert e(n_classes) == 3
    drop pp12_c0 pp12_c1 pp12_c2 yclass3
}
if _rc {
    display as error "FAIL: 3-class auto-detection"
    local errors = `errors' + 1
}
else {
    display as result "PASS: 3-class auto-detection"
}

* ---- Test 13: if restriction ----
capture noisily {
    grf_probability_forest yclass x1-x5 if x1 > 0, gen(pp13) ntrees(100) seed(42)
    quietly count if x1 > 0
    local n_sub = r(N)
    assert e(N) == `n_sub'
    drop pp13_c0 pp13_c1
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
    grf_probability_forest yclass x1-x5, gen(pp14) ntrees(100) seed(42)
    grf_probability_forest yclass x1-x5, gen(pp14) ntrees(100) seed(42) replace
    assert !missing(pp14_c0) in 1
    drop pp14_c0 pp14_c1
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
    grf_probability_forest yclass x1-x5, gen(pp15) ntrees(100) seed(123) ///
        mtry(3) minnodesize(10) samplefrac(0.7) honestyfrac(0.7)        ///
        nohonestyprune alpha(0.1) imbalancepenalty(0.5) numthreads(2)    ///
        nclasses(2) nohonesty
    assert !missing(pp15_c0) in 1
    assert !missing(pp15_c1) in 1
    assert e(mtry) == 3
    assert e(min_node) == 10
    assert reldif(e(sample_fraction), 0.7) < 1e-6
    assert e(honesty) == 0
    assert e(honesty_prune) == 0
    assert reldif(e(alpha), 0.1) < 1e-6
    assert reldif(e(imbalance_penalty), 0.5) < 1e-6
    assert e(n_classes) == 2
    drop pp15_c0 pp15_c1
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
    grf_probability_forest yclass x1-x5, gen(pp16) ntrees(100) seed(42)
    assert e(honesty) == 1
    assert e(honesty_prune) == 1
    assert e(mtry) == 0
    assert e(min_node) == 5
    assert reldif(e(sample_fraction), 0.5) < 1e-6
    assert reldif(e(alpha), 0.05) < 1e-6
    assert "`e(forest_type)'" == "probability"
    assert e(n_classes) == 2
    drop pp16_c0 pp16_c1
}
if _rc {
    display as error "FAIL: default e() values"
    local errors = `errors' + 1
}
else {
    display as result "PASS: default e() values"
}

* ---- Test 17: probabilities sum to approximately 1 ----
capture noisily {
    grf_probability_forest yclass x1-x5, gen(pp17) ntrees(100) seed(42)
    gen pp17_sum = pp17_c0 + pp17_c1
    quietly summarize pp17_sum
    * Probabilities should sum close to 1
    assert reldif(r(mean), 1.0) < 0.05
    drop pp17_c0 pp17_c1 pp17_sum
}
if _rc {
    display as error "FAIL: probabilities sum to ~1"
    local errors = `errors' + 1
}
else {
    display as result "PASS: probabilities sum to ~1"
}

* ============================================================
* Summary
* ============================================================
display ""
if `errors' > 0 {
    display as error "FAILED: `errors' errors in probability forest option tests"
    exit 1
}
else {
    display as result "ALL PROBABILITY FOREST OPTION TESTS PASSED"
}
