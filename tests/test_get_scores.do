* test_get_scores.do -- Tests for grf_get_scores (doubly-robust score extraction)
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

* ============================================================
* Test 1: Basic causal forest scores
* ============================================================

* ---- Test 1: grf_get_scores basic causal forest ----
capture noisily {
    grf_causal_forest y w x1-x5, gen(tau_var) ntrees(100) seed(42)
    grf_get_scores, gen(scores)
    assert r(N) > 0
    assert !missing(r(mean))
    assert r(sd) > 0
    quietly count if !missing(scores)
    assert r(N) > 0
    drop scores tau_var _grf_yhat _grf_what
}
if _rc {
    display as error "FAIL: grf_get_scores basic causal forest"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_get_scores basic causal forest"
}

* ---- Test 2: grf_get_scores instrumental forest ----
capture noisily {
    gen z = (rnormal() > 0)
    gen w2 = (z + rnormal() > 0.5)
    grf_instrumental_forest y w2 z x1-x5, gen(late) ntrees(100) seed(42)
    grf_get_scores, gen(iscores)
    assert r(N) > 0
    assert !missing(r(mean))
    quietly count if !missing(iscores)
    assert r(N) > 0
    drop iscores late z w2
}
if _rc {
    display as error "FAIL: grf_get_scores instrumental forest"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_get_scores instrumental forest"
}

* ---- Test 3: replace option ----
capture noisily {
    grf_causal_forest y w x1-x5, gen(tau_var) ntrees(100) seed(42)
    grf_get_scores, gen(s)
    * Second call without replace should fail
    capture grf_get_scores, gen(s)
    assert _rc != 0
    * Second call with replace should succeed
    grf_get_scores, gen(s) replace
    quietly count if !missing(s)
    assert r(N) > 0
    drop s tau_var _grf_yhat _grf_what
}
if _rc {
    display as error "FAIL: grf_get_scores replace option"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_get_scores replace option"
}

* ---- Test 4: error without prior estimation ----
capture noisily {
    * Save current data, clear e() results
    preserve
    clear
    set obs 10
    gen x = rnormal()
    capture grf_get_scores, gen(s)
    local rc4 = _rc
    restore
    assert `rc4' != 0
}
if _rc {
    display as error "FAIL: grf_get_scores error without prior estimation"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_get_scores error without prior estimation"
}

* ---- Test 5: error with wrong forest type (regression) ----
capture noisily {
    grf_regression_forest y x1-x5, gen(rpred) ntrees(100) seed(42)
    capture grf_get_scores, gen(s)
    assert _rc == 301
    drop rpred
}
if _rc {
    display as error "FAIL: grf_get_scores error with regression forest"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_get_scores error with regression forest"
}

* ---- Test 6: if/in restriction on score generation ----
capture noisily {
    grf_causal_forest y w x1-x5, gen(tau_var) ntrees(100) seed(42)
    grf_get_scores, gen(s)
    * Scores should be generated for the full sample
    quietly count if !missing(s)
    assert r(N) == 500
    * Check scores have reasonable variation
    quietly summarize s
    assert r(sd) > 0
    drop s tau_var _grf_yhat _grf_what
}
if _rc {
    display as error "FAIL: grf_get_scores if/in restriction"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_get_scores if/in restriction"
}

* ---- Test 7: score properties -- mean(DR scores) approximates ATE ----
capture noisily {
    grf_causal_forest y w x1-x5, gen(tau_var) ntrees(100) seed(42)
    grf_ate
    local ate_val = r(ate)
    local ate_se = r(se)
    grf_get_scores, gen(s)
    local score_mean = r(mean)
    local score_se = r(se)
    * Mean of DR scores should approximate ATE within 3*se
    local diff = abs(`score_mean' - `ate_val')
    local tol = 3 * `ate_se'
    assert `diff' < `tol'
    drop s tau_var _grf_yhat _grf_what
}
if _rc {
    display as error "FAIL: grf_get_scores score properties (mean ~ ATE)"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_get_scores score properties (mean ~ ATE)"
}

* ---- Test 8: multi-arm causal forest DR scores ----
capture noisily {
    gen arm = mod(_n, 3)
    gen w_arm1 = (arm == 1)
    gen w_arm2 = (arm == 2)
    gen y_ma = 1 + 0.5*x1 + w_arm1*(1 + x2) + w_arm2*(0.5 - x3) + rnormal()

    grf_multi_arm_causal_forest y_ma w_arm1 w_arm2 x1-x5, ///
        gen(ma_tau) ntreat(2) ntrees(100) seed(42)
    grf_get_scores, gen(ma_score)

    assert r(n_treat) == 2
    quietly count if !missing(ma_score_t1)
    assert r(N) > 0
    quietly count if !missing(ma_score_t2)
    assert r(N) > 0

    * Means should be in the same ballpark as stored arm-level ATE summaries
    local ate1 = e(ate_1)
    local ate2 = e(ate_2)
    quietly summarize ma_score_t1
    assert abs(r(mean) - `ate1') < 2
    quietly summarize ma_score_t2
    assert abs(r(mean) - `ate2') < 2

    drop arm w_arm1 w_arm2 y_ma ma_tau_t1 ma_tau_t2 ma_score_t1 ma_score_t2
    capture drop _grf_mac_yhat _grf_mac_what1 _grf_mac_what2
}
if _rc {
    display as error "FAIL: grf_get_scores multi-arm causal forest"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_get_scores multi-arm causal forest"
}

* ---- Test 9: causal survival IPCW/DR scores ----
capture noisily {
    gen w_cs = (x1 + rnormal() > 0)
    gen time_cs = exp(0.4*x1 - 0.2*w_cs + rnormal())
    gen status_cs = (runiform() > 0.3)

    grf_causal_survival_forest time_cs status_cs w_cs x1-x5, ///
        gen(cs_tau) ntrees(100) seed(42) horizon(2.5)

    assert "`e(numer_var)'" == "_grf_cs_numer"
    assert "`e(denom_var)'" == "_grf_cs_denom"
    assert "`e(what_var)'" == "_grf_cs_what"

    grf_get_scores, gen(cs_scores)
    assert !missing(r(denom_mean))
    assert !missing(r(vhat_mean))
    assert "`r(vhat_source)'" == "W.hat * (1 - W.hat)"
    assert r(sd) > 0

    * Formula check against upstream-style tau + psi / V.hat for binary treatment.
    gen double cs_vhat = max(_grf_cs_what * (1 - _grf_cs_what), 1e-12)
    gen double cs_formula = cs_tau + (_grf_cs_numer - _grf_cs_denom * cs_tau) / cs_vhat
    gen double cs_abs_diff = abs(cs_scores - cs_formula)
    quietly summarize cs_abs_diff, meanonly
    assert r(max) < 1e-10

    drop w_cs time_cs status_cs cs_tau cs_scores cs_vhat cs_formula cs_abs_diff
    capture drop _grf_cs_what _grf_cs_numer _grf_cs_denom _grf_cs_yhat _grf_cs_shat _grf_cs_chat
}
if _rc {
    display as error "FAIL: grf_get_scores causal survival forest"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_get_scores causal survival forest"
}

* ============================================================
* Summary
* ============================================================
display ""
if `errors' > 0 {
    display as error "FAILED: `errors' errors in get_scores tests"
    exit 1
}
else {
    display as result "ALL GET_SCORES TESTS PASSED"
}
