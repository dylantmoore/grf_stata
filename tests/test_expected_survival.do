* test_expected_survival.do -- Tests for grf_expected_survival
clear all
set more off

local errors = 0

* ============================================================
* Setup survival data
* ============================================================
clear
set obs 500
set seed 42
forvalues j = 1/5 {
    gen x`j' = rnormal()
}
gen double T_latent = -ln(runiform()) / exp(0.5 * x1)
gen double C = -ln(runiform()) / 0.5
gen double time = min(T_latent, C)
gen byte status = (T_latent <= C)
drop T_latent C

* ---- Test 1: basic expected survival time ----
capture noisily {
    grf_survival_forest time status x1-x5, gen(surv) ntrees(100) seed(42)
    grf_expected_survival, gen(etime)
    assert r(N) > 0
    assert r(mean) > 0
    quietly count if !missing(etime)
    assert r(N) > 0
    drop etime
    * Clean up survival curve variables
    local ncols = e(n_output)
    forvalues j = 1/`ncols' {
        drop surv_s`j'
    }
}
if _rc {
    display as error "FAIL: grf_expected_survival basic"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_expected_survival basic"
}

* ---- Test 2: custom predictions stub ----
capture noisily {
    grf_survival_forest time status x1-x5, gen(mysurv) ntrees(100) seed(42)
    grf_expected_survival, gen(etime2) predictions(mysurv)
    assert r(N) > 0
    assert !missing(r(mean))
    drop etime2
    local ncols = e(n_output)
    forvalues j = 1/`ncols' {
        drop mysurv_s`j'
    }
}
if _rc {
    display as error "FAIL: grf_expected_survival custom predictions stub"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_expected_survival custom predictions stub"
}

* ---- Test 3: custom grid ----
capture noisily {
    * Fit survival forest with noutput(4) so we get exactly 4 columns
    grf_survival_forest time status x1-x5, gen(sg) ntrees(100) seed(42) noutput(4)
    grf_expected_survival, gen(etime3) grid(0.5 1.0 1.5 2.0)
    assert r(N) > 0
    assert r(mean) > 0
    assert r(n_grid) == 4
    drop etime3
    forvalues j = 1/4 {
        drop sg_s`j'
    }
}
if _rc {
    display as error "FAIL: grf_expected_survival custom grid"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_expected_survival custom grid"
}

* ---- Test 4: replace option ----
capture noisily {
    grf_survival_forest time status x1-x5, gen(sr) ntrees(100) seed(42)
    grf_expected_survival, gen(etime4)
    * Second call without replace should fail
    capture grf_expected_survival, gen(etime4)
    assert _rc != 0
    * With replace should succeed
    grf_expected_survival, gen(etime4) replace
    quietly count if !missing(etime4)
    assert r(N) > 0
    drop etime4
    local ncols = e(n_output)
    forvalues j = 1/`ncols' {
        drop sr_s`j'
    }
}
if _rc {
    display as error "FAIL: grf_expected_survival replace option"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_expected_survival replace option"
}

* ---- Test 5: error without survival forest ----
capture noisily {
    * Fit a regression forest (wrong type)
    grf_regression_forest time x1-x5, gen(rpred) ntrees(100) seed(42)
    capture grf_expected_survival, gen(etime5)
    assert _rc == 301
    drop rpred
}
if _rc {
    display as error "FAIL: grf_expected_survival error without survival forest"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_expected_survival error without survival forest"
}

* ---- Test 6: all expected survival times are positive ----
capture noisily {
    grf_survival_forest time status x1-x5, gen(sp) ntrees(100) seed(42)
    grf_expected_survival, gen(etime6)
    quietly summarize etime6
    assert r(min) > 0
    drop etime6
    local ncols = e(n_output)
    forvalues j = 1/`ncols' {
        drop sp_s`j'
    }
}
if _rc {
    display as error "FAIL: grf_expected_survival positive values"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_expected_survival positive values"
}

* ---- Test 7: heterogeneity -- variance of E[T|X] > 0 ----
capture noisily {
    grf_survival_forest time status x1-x5, gen(sh) ntrees(100) seed(42)
    grf_expected_survival, gen(etime7)
    assert r(sd) > 0
    drop etime7
    local ncols = e(n_output)
    forvalues j = 1/`ncols' {
        drop sh_s`j'
    }
}
if _rc {
    display as error "FAIL: grf_expected_survival heterogeneity (variance > 0)"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_expected_survival heterogeneity (variance > 0)"
}

* ============================================================
* Summary
* ============================================================
display ""
if `errors' > 0 {
    display as error "FAILED: `errors' errors in expected survival tests"
    exit 1
}
else {
    display as result "ALL EXPECTED SURVIVAL TESTS PASSED"
}
