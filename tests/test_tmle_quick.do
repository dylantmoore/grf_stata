* test_tmle_quick.do -- Quick verification of TMLE fixes and equalizeclusterweights
* Runs only E4.1-E4.3 (equalize) and E5 tests

clear all
set more off

local errors = 0

* ============================================================
* E4: equalizeclusterweights (subset)
* ============================================================

display as text ""
display as text "E4: equalizeclusterweights (subset)"
display as text "============================================================"

* ---- E4.1: equalize.cluster.weights on regression forest ----
capture noisily {
    clear
    set obs 300
    set seed 42
    gen cluster_id = cond(_n <= 100, 1, 2)
    gen x1 = rnormal()
    gen x2 = rnormal()
    gen y = 2*x1 + x2 + rnormal()

    grf_regression_forest y x1 x2, gen(pred_eq) ntrees(100) seed(42) ///
        cluster(cluster_id) equalizeclusterweights
    assert !missing(pred_eq) in 1
    assert "`e(cluster_var)'" == "cluster_id"

    grf_regression_forest y x1 x2, gen(pred_neq) ntrees(100) seed(42) ///
        cluster(cluster_id) replace
    gen diff_eq = pred_eq - pred_neq
    quietly summarize diff_eq
    assert r(sd) > 0
    drop pred_eq pred_neq diff_eq
}
if _rc {
    display as error "FAIL: E4.1 equalize.cluster.weights regression forest"
    local errors = `errors' + 1
}
else {
    display as result "PASS: E4.1 equalize.cluster.weights regression forest"
}

* ---- E4.3: equalize.cluster.weights on causal forest ----
capture noisily {
    clear
    set obs 300
    set seed 42
    gen cluster_id = ceil(_n / 10)
    gen x1 = rnormal()
    gen x2 = rnormal()
    gen w = rbinomial(1, 0.5)
    gen y = 2*x1 + w + rnormal()

    grf_causal_forest y w x1 x2, gen(cate_eq) ntrees(100) seed(42) ///
        cluster(cluster_id) equalizeclusterweights
    assert !missing(cate_eq) in 1
    assert "`e(cluster_var)'" == "cluster_id"
    drop cate_eq _grf_yhat _grf_what
}
if _rc {
    display as error "FAIL: E4.3 equalize.cluster.weights causal forest"
    local errors = `errors' + 1
}
else {
    display as result "PASS: E4.3 equalize.cluster.weights causal forest"
}


* ============================================================
* E5: TMLE tests
* ============================================================

display as text ""
display as text "E5: TMLE tests"
display as text "============================================================"

* ---- Setup ----
clear
set obs 200
set seed 42
gen x1 = rnormal()
gen x2 = rnormal()
gen w = rbinomial(1, normal(0.3*x1))
gen y = 2*x1 + x1*w + rnormal()

grf_causal_forest y w x1 x2, gen(tau_e5) ntrees(100) seed(42)

* ---- E5.1: TMLE basic ----
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

* ---- E5.2: TMLE vs AIPW differ ----
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

* ---- E5.3: TMLE + target.sample(treated) ----
capture noisily {
    grf_ate, method(TMLE) targetsample(treated)
    assert !missing(r(ate))
    assert !missing(r(se))
    assert r(se) > 0
}
if _rc {
    display as error "FAIL: E5.3 TMLE target.sample(treated)"
    local errors = `errors' + 1
}
else {
    display as result "PASS: E5.3 TMLE target.sample(treated)"
}

* ---- E5.4: TMLE + target.sample(control) ----
capture noisily {
    grf_ate, method(TMLE) targetsample(control)
    assert !missing(r(ate))
    assert !missing(r(se))
    assert r(se) > 0
}
if _rc {
    display as error "FAIL: E5.4 TMLE target.sample(control)"
    local errors = `errors' + 1
}
else {
    display as result "PASS: E5.4 TMLE target.sample(control)"
}

* ---- E5.5: TMLE + overlap falls back to AIPW ----
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

* ---- E5.6: TMLE with clustered SEs ----
capture noisily {
    drop tau_e5 _grf_yhat _grf_what
    clear
    set obs 200
    set seed 42
    gen cluster_id = ceil(_n / 20)
    gen x1 = rnormal()
    gen x2 = rnormal()
    gen w = rbinomial(1, normal(0.3*x1))
    gen y = 2*x1 + x1*w + rnormal()

    grf_causal_forest y w x1 x2, gen(tau_cl) ntrees(100) seed(42) ///
        cluster(cluster_id)
    grf_ate, method(TMLE)
    assert !missing(r(ate))
    assert !missing(r(se))
    assert r(se) > 0
    drop tau_cl _grf_yhat _grf_what
}
if _rc {
    display as error "FAIL: E5.6 TMLE clustered SEs"
    local errors = `errors' + 1
}
else {
    display as result "PASS: E5.6 TMLE clustered SEs"
}

* ---- E5.7: TMLE error on invalid method ----
capture {
    clear
    set obs 100
    set seed 42
    gen x1 = rnormal()
    gen w = rbinomial(1, 0.5)
    gen y = x1 + w + rnormal()

    grf_causal_forest y w x1, gen(tau_err) ntrees(50) seed(42)
    grf_ate, method(INVALID)
}
if _rc {
    display as result "PASS: E5.7 error on invalid method"
}
else {
    display as error "FAIL: E5.7 should error on invalid method"
    local errors = `errors' + 1
}

* ---- E5.7b: TMLE + debiasingweights errors out ----
capture {
    clear
    set obs 200
    set seed 42
    gen x1 = rnormal()
    gen x2 = rnormal()
    gen w = rbinomial(1, 0.5)
    gen y = x1 + w + rnormal()
    gen dbwt = abs(x1) + 0.5

    grf_causal_forest y w x1 x2, gen(tau_db_err) ntrees(50) seed(42)
    grf_ate, method(TMLE) debiasingweights(dbwt)
}
if _rc {
    display as result "PASS: E5.7b TMLE + debiasingweights errors"
}
else {
    display as error "FAIL: E5.7b should error with method(TMLE) debiasingweights()"
    local errors = `errors' + 1
}

* ============================================================
* Summary
* ============================================================
display as text ""
display as text "============================================================"
display as text "Summary"
display as text "============================================================"

if `errors' > 0 {
    display as error "FAILED: `errors' error(s) in TMLE verification tests"
    exit 1
}
else {
    display as result "ALL TMLE VERIFICATION TESTS PASSED (`errors' errors)"
}
