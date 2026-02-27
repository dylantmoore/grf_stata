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
