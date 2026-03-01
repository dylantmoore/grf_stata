* Minimal TMLE verification test
clear all
set more off
local errors = 0

* Setup: fit one causal forest
clear
set obs 200
set seed 42
gen x1 = rnormal()
gen x2 = rnormal()
gen w = rbinomial(1, normal(0.3*x1))
gen y = 2*x1 + x1*w + rnormal()

grf_causal_forest y w x1 x2, gen(tau_e5) ntrees(100) seed(42)

* E5.1: TMLE basic
capture noisily {
    grf_ate, method(TMLE)
    assert !missing(r(ate))
    assert !missing(r(se))
    assert r(se) > 0
    assert "`r(method)'" == "TMLE"
}
if _rc {
    display as error "FAIL: E5.1 TMLE basic"
    local errors = `errors' + 1
}
else {
    display as result "PASS: E5.1 TMLE basic"
}

* E5.2: TMLE vs AIPW differ
capture noisily {
    grf_ate, method(AIPW)
    local ate_aipw = r(ate)
    grf_ate, method(TMLE)
    local ate_tmle = r(ate)
    assert reldif(`ate_aipw', `ate_tmle') > 1e-10
}
if _rc {
    display as error "FAIL: E5.2 TMLE vs AIPW differ"
    local errors = `errors' + 1
}
else {
    display as result "PASS: E5.2 TMLE vs AIPW differ"
}

* E5.3: ATT
capture noisily {
    grf_ate, method(TMLE) targetsample(treated)
    assert !missing(r(ate))
    assert !missing(r(se))
    assert r(se) > 0
}
if _rc {
    display as error "FAIL: E5.3 TMLE ATT"
    local errors = `errors' + 1
}
else {
    display as result "PASS: E5.3 TMLE ATT"
}

* E5.4: ATC
capture noisily {
    grf_ate, method(TMLE) targetsample(control)
    assert !missing(r(ate))
    assert !missing(r(se))
    assert r(se) > 0
}
if _rc {
    display as error "FAIL: E5.4 TMLE ATC"
    local errors = `errors' + 1
}
else {
    display as result "PASS: E5.4 TMLE ATC"
}

* E5.5: overlap fallback
capture noisily {
    grf_ate, method(TMLE) targetsample(overlap)
    assert !missing(r(ate))
    assert !missing(r(se))
    assert "`r(method)'" == "AIPW"
}
if _rc {
    display as error "FAIL: E5.5 TMLE overlap fallback"
    local errors = `errors' + 1
}
else {
    display as result "PASS: E5.5 TMLE overlap fallback"
}

* E5.7b: TMLE + debiasingweights errors
capture {
    gen dbwt = abs(x1) + 0.5
    grf_ate, method(TMLE) debiasingweights(dbwt)
}
if _rc {
    display as result "PASS: E5.7b TMLE+debiasingweights errors"
}
else {
    display as error "FAIL: E5.7b should error"
    local errors = `errors' + 1
}

* Summary
display as text ""
if `errors' > 0 {
    display as error "FAILED: `errors' error(s)"
    exit 1
}
else {
    display as result "ALL TMLE TESTS PASSED"
}
