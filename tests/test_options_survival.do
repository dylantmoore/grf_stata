* test_options_survival.do -- Comprehensive option tests for grf_survival_forest
* Tests every option individually and in combination

clear all
set more off

local errors = 0

* ============================================================
* Setup survival data
* Syntax: grf_survival_forest time status X1..Xp, gen(...)
* ============================================================
clear
set obs 500
set seed 42
forvalues j = 1/5 {
    gen x`j' = rnormal()
}
* Generate survival time (positive) and censoring indicator
gen time = exp(x1 + 0.5*x2 + rnormal())
gen status = (uniform() > 0.3)

* ---- Test 1: mtry(3) ----
capture noisily {
    grf_survival_forest time status x1-x5, gen(sp1) ntrees(100) seed(42) mtry(3)
    assert !missing(sp1_s1) in 1
    assert e(mtry) == 3
    forvalues j = 1/20 {
        capture drop sp1_s`j'
    }
}
if _rc {
    display as error "FAIL: mtry(3)"
    local errors = `errors' + 1
}
else {
    display as result "PASS: mtry(3)"
}

* ---- Test 2: minnodesize(20) ----
capture noisily {
    grf_survival_forest time status x1-x5, gen(sp2) ntrees(100) seed(42) minnodesize(20)
    assert !missing(sp2_s1) in 1
    assert e(min_node) == 20
    forvalues j = 1/20 {
        capture drop sp2_s`j'
    }
}
if _rc {
    display as error "FAIL: minnodesize(20)"
    local errors = `errors' + 1
}
else {
    display as result "PASS: minnodesize(20)"
}

* ---- Test 3: samplefrac(0.7) ----
capture noisily {
    grf_survival_forest time status x1-x5, gen(sp3) ntrees(100) seed(42) samplefrac(0.7)
    assert !missing(sp3_s1) in 1
    assert reldif(e(sample_fraction), 0.7) < 1e-6
    forvalues j = 1/20 {
        capture drop sp3_s`j'
    }
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
    grf_survival_forest time status x1-x5, gen(sp4) ntrees(100) seed(42) nohonesty
    assert !missing(sp4_s1) in 1
    assert e(honesty) == 0
    forvalues j = 1/20 {
        capture drop sp4_s`j'
    }
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
    grf_survival_forest time status x1-x5, gen(sp5) ntrees(100) seed(42) honestyfrac(0.7)
    assert !missing(sp5_s1) in 1
    assert reldif(e(honesty_fraction), 0.7) < 1e-6
    forvalues j = 1/20 {
        capture drop sp5_s`j'
    }
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
    grf_survival_forest time status x1-x5, gen(sp6) ntrees(100) seed(42) nohonestyprune
    assert !missing(sp6_s1) in 1
    assert e(honesty_prune) == 0
    forvalues j = 1/20 {
        capture drop sp6_s`j'
    }
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
    grf_survival_forest time status x1-x5, gen(sp7) ntrees(100) seed(42) alpha(0.1)
    assert !missing(sp7_s1) in 1
    assert reldif(e(alpha), 0.1) < 1e-6
    forvalues j = 1/20 {
        capture drop sp7_s`j'
    }
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
    grf_survival_forest time status x1-x5, gen(sp8) ntrees(100) seed(42) imbalancepenalty(0.5)
    assert !missing(sp8_s1) in 1
    assert reldif(e(imbalance_penalty), 0.5) < 1e-6
    forvalues j = 1/20 {
        capture drop sp8_s`j'
    }
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
    grf_survival_forest time status x1-x5, gen(sp9) ntrees(100) seed(42) cigroupsize(2) samplefrac(0.4)
    assert !missing(sp9_s1) in 1
    assert e(ci_group_size) == 2
    forvalues j = 1/20 {
        capture drop sp9_s`j'
    }
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
    grf_survival_forest time status x1-x5, gen(sp10) ntrees(100) seed(42) numthreads(2)
    assert !missing(sp10_s1) in 1
    forvalues j = 1/20 {
        capture drop sp10_s`j'
    }
}
if _rc {
    display as error "FAIL: numthreads(2)"
    local errors = `errors' + 1
}
else {
    display as result "PASS: numthreads(2)"
}

* ---- Test 11: numfailures(5) ----
capture noisily {
    grf_survival_forest time status x1-x5, gen(sp11) ntrees(100) seed(42) numfailures(5)
    assert !missing(sp11_s1) in 1
    forvalues j = 1/20 {
        capture drop sp11_s`j'
    }
}
if _rc {
    display as error "FAIL: numfailures(5)"
    local errors = `errors' + 1
}
else {
    display as result "PASS: numfailures(5)"
}

* ---- Test 12: predtype(0) (Nelson-Aalen) ----
capture noisily {
    grf_survival_forest time status x1-x5, gen(sp12) ntrees(100) seed(42) predtype(0)
    assert !missing(sp12_s1) in 1
    assert e(pred_type) == 0
    forvalues j = 1/20 {
        capture drop sp12_s`j'
    }
}
if _rc {
    display as error "FAIL: predtype(0) Nelson-Aalen"
    local errors = `errors' + 1
}
else {
    display as result "PASS: predtype(0) Nelson-Aalen"
}

* ---- Test 13: predtype(1) (Kaplan-Meier, default) ----
capture noisily {
    grf_survival_forest time status x1-x5, gen(sp13) ntrees(100) seed(42) predtype(1)
    assert !missing(sp13_s1) in 1
    assert e(pred_type) == 1
    forvalues j = 1/20 {
        capture drop sp13_s`j'
    }
}
if _rc {
    display as error "FAIL: predtype(1) Kaplan-Meier"
    local errors = `errors' + 1
}
else {
    display as result "PASS: predtype(1) Kaplan-Meier"
}

* ---- Test 14: noutput(10) ----
capture noisily {
    grf_survival_forest time status x1-x5, gen(sp14) ntrees(100) seed(42) noutput(10)
    assert !missing(sp14_s1) in 1
    assert !missing(sp14_s10) in 1
    assert e(n_output) == 10
    forvalues j = 1/10 {
        capture drop sp14_s`j'
    }
}
if _rc {
    display as error "FAIL: noutput(10)"
    local errors = `errors' + 1
}
else {
    display as result "PASS: noutput(10)"
}

* ---- Test 15: if restriction ----
capture noisily {
    grf_survival_forest time status x1-x5 if x1 > 0, gen(sp15) ntrees(100) seed(42)
    quietly count if x1 > 0
    local n_sub = r(N)
    assert e(N) == `n_sub'
    forvalues j = 1/20 {
        capture drop sp15_s`j'
    }
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
    grf_survival_forest time status x1-x5, gen(sp16) ntrees(100) seed(42)
    grf_survival_forest time status x1-x5, gen(sp16) ntrees(100) seed(42) replace
    assert !missing(sp16_s1) in 1
    forvalues j = 1/20 {
        capture drop sp16_s`j'
    }
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
    grf_survival_forest time status x1-x5, gen(sp17) ntrees(100) seed(123) ///
        mtry(3) minnodesize(20) samplefrac(0.4) honestyfrac(0.7)          ///
        nohonestyprune alpha(0.1) imbalancepenalty(0.5) cigroupsize(2)     ///
        numthreads(2) noutput(10) numfailures(5) predtype(0)
    assert !missing(sp17_s1) in 1
    assert !missing(sp17_s5) in 1
    assert e(mtry) == 3
    assert e(min_node) == 20
    assert reldif(e(sample_fraction), 0.4) < 1e-6
    assert e(honesty_prune) == 0
    assert reldif(e(alpha), 0.1) < 1e-6
    assert reldif(e(imbalance_penalty), 0.5) < 1e-6
    assert e(ci_group_size) == 2
    assert e(n_output) == 10
    assert e(pred_type) == 0
    forvalues j = 1/10 {
        capture drop sp17_s`j'
    }
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
    grf_survival_forest time status x1-x5, gen(sp18) ntrees(100) seed(42)
    assert e(honesty) == 1
    assert e(honesty_prune) == 1
    assert e(mtry) == 0
    assert e(min_node) == 15
    assert reldif(e(sample_fraction), 0.5) < 1e-6
    assert reldif(e(alpha), 0.05) < 1e-6
    assert e(n_output) == 20
    assert e(pred_type) == 1
    assert "`e(forest_type)'" == "survival"
    assert "`e(timevar)'" == "time"
    assert "`e(statusvar)'" == "status"
    assert !missing(e(n_events))
    assert !missing(e(n_censored))
    forvalues j = 1/20 {
        capture drop sp18_s`j'
    }
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
    display as error "FAILED: `errors' errors in survival forest option tests"
    exit 1
}
else {
    display as result "ALL SURVIVAL FOREST OPTION TESTS PASSED"
}
