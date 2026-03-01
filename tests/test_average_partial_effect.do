* test_average_partial_effect.do -- Tests for grf_average_partial_effect
* Tests: basic continuous W, if/in restrictions, debiasing weights,
*        nocalibrate, wrong forest error, binary W vs ATE, CI properties

clear all
set more off

local errors = 0

* ============================================================
* Setup data (continuous treatment)
* ============================================================
clear
set obs 500
set seed 42
forvalues j = 1/5 {
    gen x`j' = rnormal()
}
gen w = rnormal()
gen y = x1 + w * x2 + rnormal()

* ============================================================
* Test 1: Basic test with continuous W
* ============================================================
capture noisily {
    grf_causal_forest y w x1-x5, gen(tau_ape1) ntrees(100) seed(42)
    grf_average_partial_effect
    assert !missing(r(estimate))
    assert !missing(r(std_err))
    assert r(std_err) > 0
    assert !missing(r(ci_lower))
    assert !missing(r(ci_upper))
    assert !missing(r(pvalue))
    assert r(N) == 500
    drop tau_ape1 _grf_yhat _grf_what
}
if _rc {
    display as error "FAIL: APE basic continuous W"
    local errors = `errors' + 1
}
else {
    display as result "PASS: APE basic continuous W"
}

* ============================================================
* Test 2: if/in restrictions
* ============================================================
capture noisily {
    grf_causal_forest y w x1-x5, gen(tau_ape2) ntrees(100) seed(42)
    grf_average_partial_effect if x1 > 0
    local saved_N = r(N)
    assert !missing(r(estimate))
    assert `saved_N' < 500
    assert `saved_N' > 0
    drop tau_ape2 _grf_yhat _grf_what
}
if _rc {
    display as error "FAIL: APE with if restriction"
    local errors = `errors' + 1
}
else {
    display as result "PASS: APE with if restriction"
}

* ============================================================
* Test 3: debiasing weights option
* ============================================================
capture noisily {
    grf_causal_forest y w x1-x5, gen(tau_ape3) ntrees(100) seed(42)
    gen dbw = abs(x1) + 0.1
    grf_average_partial_effect, debiasweights(dbw)
    assert !missing(r(estimate))
    assert !missing(r(std_err))
    assert r(std_err) > 0
    assert r(N) == 500
    drop tau_ape3 _grf_yhat _grf_what dbw
}
if _rc {
    display as error "FAIL: APE with debiasing weights"
    local errors = `errors' + 1
}
else {
    display as result "PASS: APE with debiasing weights"
}

* ============================================================
* Test 4: nocalibrate option
* ============================================================
capture noisily {
    grf_causal_forest y w x1-x5, gen(tau_ape4) ntrees(100) seed(42)
    grf_average_partial_effect, nocalibrate
    assert !missing(r(estimate))
    assert !missing(r(std_err))
    assert r(std_err) > 0
    drop tau_ape4 _grf_yhat _grf_what
}
if _rc {
    display as error "FAIL: APE nocalibrate"
    local errors = `errors' + 1
}
else {
    display as result "PASS: APE nocalibrate"
}

* ============================================================
* Test 5: Error with wrong forest type (regression)
* ============================================================
capture noisily {
    grf_regression_forest y x1-x5, gen(rpred) ntrees(100) seed(42)
    capture grf_average_partial_effect
    assert _rc == 301
    drop rpred
}
if _rc {
    display as error "FAIL: APE error with regression forest"
    local errors = `errors' + 1
}
else {
    display as result "PASS: APE error with regression forest"
}

* ============================================================
* Test 6: Binary W -- APE should approximate ATE
* ============================================================
capture noisily {
    * Regenerate data with binary treatment
    drop y w
    gen w = (x1 + rnormal() > 0)
    gen y = 2*x1 + x2^2 + 0.5*w + rnormal()

    grf_causal_forest y w x1-x5, gen(tau_ape6) ntrees(100) seed(42)

    * Get ATE
    grf_ate
    local ate_est = r(ate)

    * Get APE
    grf_average_partial_effect
    local ape_est = r(estimate)

    * They should be in the same ballpark (within 1.0 of each other)
    assert abs(`ape_est' - `ate_est') < 1.0

    drop tau_ape6 _grf_yhat _grf_what y w
}
if _rc {
    display as error "FAIL: APE vs ATE for binary W"
    local errors = `errors' + 1
}
else {
    display as result "PASS: APE vs ATE for binary W"
}

* ============================================================
* Test 7: Confidence interval properties
* ============================================================
capture noisily {
    * Restore continuous treatment data
    gen w = rnormal()
    gen y = x1 + w * x2 + rnormal()

    grf_causal_forest y w x1-x5, gen(tau_ape7) ntrees(100) seed(42)
    grf_average_partial_effect
    local saved_est = r(estimate)
    local saved_lo  = r(ci_lower)
    local saved_hi  = r(ci_upper)
    local saved_p   = r(pvalue)

    assert `saved_lo' < `saved_est'
    assert `saved_est' < `saved_hi'
    assert `saved_p' >= 0 & `saved_p' <= 1

    drop tau_ape7 _grf_yhat _grf_what y w
}
if _rc {
    display as error "FAIL: APE confidence interval properties"
    local errors = `errors' + 1
}
else {
    display as result "PASS: APE confidence interval properties"
}

* ============================================================
* Test 8: numtreesvariance() changes SE
* ============================================================
capture noisily {
    * Restore continuous treatment data
    capture drop w y
    gen w = rnormal()
    gen y = x1 + w * x2 + rnormal()

    grf_causal_forest y w x1-x5, gen(tau_ape8) ntrees(100) seed(42)

    * Default numtreesvariance (500)
    grf_average_partial_effect
    local se_default = r(std_err)

    * Custom numtreesvariance (100)
    grf_average_partial_effect, numtreesvariance(100)
    local se_100 = r(std_err)

    * Both should be positive
    assert `se_default' > 0
    assert `se_100' > 0

    drop tau_ape8 _grf_yhat _grf_what y w
}
if _rc {
    display as error "FAIL: APE numtreesvariance changes SE"
    local errors = `errors' + 1
}
else {
    display as result "PASS: APE numtreesvariance changes SE"
}

* ============================================================
* Test 9: APE with in restriction
* ============================================================
capture noisily {
    capture drop w y
    gen w = rnormal()
    gen y = x1 + w * x2 + rnormal()

    grf_causal_forest y w x1-x5, gen(tau_ape9) ntrees(100) seed(42)
    grf_average_partial_effect in 1/400
    local saved_N = r(N)
    assert !missing(r(estimate))
    assert `saved_N' <= 400
    assert `saved_N' > 0
    drop tau_ape9 _grf_yhat _grf_what y w
}
if _rc {
    display as error "FAIL: APE with in restriction"
    local errors = `errors' + 1
}
else {
    display as result "PASS: APE with in restriction"
}

* ============================================================
* Summary
* ============================================================
display ""
if `errors' > 0 {
    display as error "FAILED: `errors' errors in average partial effect tests"
    exit 1
}
else {
    display as result "ALL AVERAGE PARTIAL EFFECT TESTS PASSED"
}
