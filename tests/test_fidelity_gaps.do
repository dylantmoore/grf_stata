* test_fidelity_gaps.do -- Tests for all fidelity gap fixes (Groups A-D)
*
* E1: Changed defaults (survival ntrees, lm stabilize, instrumental stabilize)
* E2: Wired parameters (fast.logrank, gradient_weights, ll.split.variables)
* E3: User-supplied nuisance estimates
* E4: New features (equalize.cluster.weights, debiasing.weights, decay.exponent)
*
* Run: do tests/test_fidelity_gaps.do

clear all
set more off

local errors = 0

* ============================================================
* E1: Changed Defaults
* ============================================================

display as text ""
display as text "============================================================"
display as text "E1: Changed Defaults"
display as text "============================================================"

* ---- E1.1: Survival forest default ntrees = 1000 ----
capture noisily {
    clear
    set obs 200
    set seed 42
    gen x1 = rnormal()
    gen x2 = rnormal()
    gen time = exp(x1 + 0.5*x2 + rnormal())
    gen status = (uniform() > 0.3)

    grf_survival_forest time status x1 x2, gen(sv_def) seed(42)
    assert e(n_trees) == 1000
    forvalues j = 1/20 {
        capture drop sv_def_s`j'
    }
}
if _rc {
    display as error "FAIL: E1.1 survival default ntrees=1000"
    local errors = `errors' + 1
}
else {
    display as result "PASS: E1.1 survival default ntrees=1000"
}

* ---- E1.2: LM forest default stabilize = OFF (0) ----
capture noisily {
    clear
    set obs 300
    set seed 42
    gen x1 = rnormal()
    gen x2 = rnormal()
    gen w1 = rnormal()
    gen y = 2*x1 + 0.5*w1*x1 + rnormal()

    grf_lm_forest y w1, xvars(x1 x2) gen(lm_def) ntrees(100) seed(42)
    assert e(stabilize) == 0
    drop lm_def_1 _grf_lm_yhat _grf_lm_what1
}
if _rc {
    display as error "FAIL: E1.2 lm_forest default stabilize=0"
    local errors = `errors' + 1
}
else {
    display as result "PASS: E1.2 lm_forest default stabilize=0"
}

* ---- E1.3: Instrumental forest default stabilize = ON (1) ----
capture noisily {
    clear
    set obs 300
    set seed 42
    gen x1 = rnormal()
    gen x2 = rnormal()
    gen z = (x1 + rnormal() > 0)
    gen w = 0.5*z + 0.3*x2 + (rnormal() > 0)
    gen y = 2*w + x1 + rnormal()

    grf_instrumental_forest y w z x1 x2, gen(iv_def) ntrees(100) seed(42)
    assert e(stabilize_splits) == 1
    drop iv_def
}
if _rc {
    display as error "FAIL: E1.3 instrumental default stabilize=1"
    local errors = `errors' + 1
}
else {
    display as result "PASS: E1.3 instrumental default stabilize=1"
}


* ============================================================
* E2: Wired Parameters
* ============================================================

display as text ""
display as text "============================================================"
display as text "E2: Wired Parameters"
display as text "============================================================"

* ---- E2.1: fast.logrank default ON for survival ----
capture noisily {
    clear
    set obs 200
    set seed 42
    gen x1 = rnormal()
    gen x2 = rnormal()
    gen time = exp(x1 + 0.5*x2 + rnormal())
    gen status = (uniform() > 0.3)

    * Default: fast.logrank ON -- should run without error
    grf_survival_forest time status x1 x2, gen(sv_fl1) ntrees(100) seed(42)
    assert !missing(sv_fl1_s1) in 1
    forvalues j = 1/20 {
        capture drop sv_fl1_s`j'
    }
}
if _rc {
    display as error "FAIL: E2.1 survival fast.logrank default ON"
    local errors = `errors' + 1
}
else {
    display as result "PASS: E2.1 survival fast.logrank default ON"
}

* ---- E2.2: nofastlogrank disables fast logrank ----
capture noisily {
    grf_survival_forest time status x1 x2, gen(sv_fl2) ntrees(100) seed(42) nofastlogrank
    assert !missing(sv_fl2_s1) in 1
    forvalues j = 1/20 {
        capture drop sv_fl2_s`j'
    }
}
if _rc {
    display as error "FAIL: E2.2 survival nofastlogrank"
    local errors = `errors' + 1
}
else {
    display as result "PASS: E2.2 survival nofastlogrank"
}

* ---- E2.3: gradient_weights for lm_forest ----
capture noisily {
    clear
    set obs 300
    set seed 42
    gen x1 = rnormal()
    gen x2 = rnormal()
    gen w1 = rnormal()
    gen w2 = rnormal()
    gen y = 2*x1 + 0.5*w1*x1 + 0.3*w2*x2 + rnormal()

    * Specify gradient weights (one per regressor)
    grf_lm_forest y w1 w2, xvars(x1 x2) gen(lm_gw) ntrees(100) seed(42) ///
        gradientweights(0.7 0.3)
    assert !missing(lm_gw_1) in 1
    assert !missing(lm_gw_2) in 1
    drop lm_gw_1 lm_gw_2 _grf_lm_yhat _grf_lm_what1 _grf_lm_what2
}
if _rc {
    display as error "FAIL: E2.3 lm_forest gradient_weights"
    local errors = `errors' + 1
}
else {
    display as result "PASS: E2.3 lm_forest gradient_weights"
}

* ---- E2.4: gradient_weights changes predictions ----
capture noisily {
    * Fit without gradient_weights
    grf_lm_forest y w1 w2, xvars(x1 x2) gen(lm_gw0) ntrees(100) seed(42)

    * Fit with gradient_weights
    grf_lm_forest y w1 w2, xvars(x1 x2) gen(lm_gw1) ntrees(100) seed(42) ///
        gradientweights(0.9 0.1) replace
    gen diff_gw = lm_gw0_1 - lm_gw1_1
    quietly summarize diff_gw
    * Predictions should differ (gradient weights change the splitting criterion)
    assert r(sd) > 0
    drop lm_gw0_1 lm_gw0_2 lm_gw1_1 lm_gw1_2 diff_gw
    capture drop _grf_lm_yhat _grf_lm_what1 _grf_lm_what2
}
if _rc {
    display as error "FAIL: E2.4 gradient_weights changes predictions"
    local errors = `errors' + 1
}
else {
    display as result "PASS: E2.4 gradient_weights changes predictions"
}

* ---- E2.5: ll.split.variables varlist for ll_regression_forest ----
capture noisily {
    clear
    set obs 300
    set seed 42
    gen x1 = rnormal()
    gen x2 = rnormal()
    gen x3 = rnormal()
    gen y = 2*x1 + x2^2 + rnormal()

    * Specify only x1 and x2 for LL splitting
    grf_ll_regression_forest y x1 x2 x3, gen(llsv1) ntrees(100) seed(42) ///
        llsplitvars(x1 x2)
    assert !missing(llsv1) in 1
    * llsplitvars implies llsplit
    assert e(enable_ll_split) == 1
    drop llsv1
}
if _rc {
    display as error "FAIL: E2.5 ll.split.variables varlist"
    local errors = `errors' + 1
}
else {
    display as result "PASS: E2.5 ll.split.variables varlist"
}

* ---- E2.6: llsplitvars with invalid variable errors ----
capture {
    grf_ll_regression_forest y x1 x2 x3, gen(llsv_err) ntrees(100) seed(42) ///
        llsplitvars(x1 nonexistent)
}
if _rc {
    display as result "PASS: E2.6 llsplitvars error on invalid variable"
}
else {
    display as error "FAIL: E2.6 should error on invalid variable in llsplitvars"
    local errors = `errors' + 1
}

* ---- E2.7: llsplitvars with explicit llsplit ----
capture noisily {
    grf_ll_regression_forest y x1 x2 x3, gen(llsv3) ntrees(100) seed(42) ///
        llsplit llsplitvars(x1)
    assert !missing(llsv3) in 1
    assert e(enable_ll_split) == 1
    drop llsv3
}
if _rc {
    display as error "FAIL: E2.7 llsplitvars with explicit llsplit"
    local errors = `errors' + 1
}
else {
    display as result "PASS: E2.7 llsplitvars with explicit llsplit"
}


* ============================================================
* E3: User-Supplied Nuisance Estimates
* ============================================================

display as text ""
display as text "============================================================"
display as text "E3: User-Supplied Nuisance Estimates"
display as text "============================================================"

* ---- E3.1: causal_forest with yhatinput/whatinput ----
capture noisily {
    clear
    set obs 400
    set seed 42
    gen x1 = rnormal()
    gen x2 = rnormal()
    gen w = rbinomial(1, normal(x1))
    gen y = 2*x1 + x1*w + rnormal()

    * Step 1: Fit causal forest (auto-computes nuisance)
    grf_causal_forest y w x1 x2, gen(tau1) ntrees(200) seed(42)
    * Save the nuisance variables
    gen my_yhat = _grf_yhat
    gen my_what = _grf_what

    * Step 2: Refit with user-supplied nuisance
    grf_causal_forest y w x1 x2, gen(tau2) ntrees(200) seed(42) ///
        yhatinput(my_yhat) whatinput(my_what) replace
    assert !missing(tau2) in 1
    * Predictions should be similar since nuisance is the same
    gen tau_diff = abs(tau1 - tau2)
    quietly summarize tau_diff
    * Allow some tolerance for tree randomness
    assert r(mean) < 1.0

    drop tau1 tau2 my_yhat my_what tau_diff _grf_yhat _grf_what
}
if _rc {
    display as error "FAIL: E3.1 causal_forest yhatinput/whatinput"
    local errors = `errors' + 1
}
else {
    display as result "PASS: E3.1 causal_forest yhatinput/whatinput"
}

* ---- E3.2: causal_forest error if only yhatinput specified ----
capture {
    clear
    set obs 200
    set seed 42
    gen x1 = rnormal()
    gen w = rbinomial(1, 0.5)
    gen y = x1 + w + rnormal()
    gen fake_yhat = x1

    grf_causal_forest y w x1, gen(tau_err) ntrees(100) seed(42) ///
        yhatinput(fake_yhat)
}
if _rc {
    display as result "PASS: E3.2 error when only yhatinput specified"
}
else {
    display as error "FAIL: E3.2 should error when only yhatinput specified"
    local errors = `errors' + 1
}

* ---- E3.3: causal_forest with yhatgenerate/whatgenerate ----
capture noisily {
    clear
    set obs 300
    set seed 42
    gen x1 = rnormal()
    gen x2 = rnormal()
    gen w = rbinomial(1, 0.5)
    gen y = 2*x1 + w + rnormal()

    grf_causal_forest y w x1 x2, gen(tau3) ntrees(100) seed(42) ///
        yhatgenerate(custom_yhat) whatgenerate(custom_what)
    assert !missing(custom_yhat) in 1
    assert !missing(custom_what) in 1
    drop tau3 custom_yhat custom_what _grf_yhat _grf_what
}
if _rc {
    display as error "FAIL: E3.3 causal_forest yhatgenerate/whatgenerate"
    local errors = `errors' + 1
}
else {
    display as result "PASS: E3.3 causal_forest yhatgenerate/whatgenerate"
}

* ---- E3.4: instrumental_forest with nuisance inputs ----
capture noisily {
    clear
    set obs 400
    set seed 42
    gen x1 = rnormal()
    gen x2 = rnormal()
    gen z = (x1 + rnormal() > 0)
    gen w = 0.5*z + 0.3*x2 + (rnormal() > 0)
    gen y = 2*w + x1 + rnormal()

    * Step 1: Fit with auto-nuisance
    grf_instrumental_forest y w z x1 x2, gen(iv_n1) ntrees(100) seed(42) ///
        yhatgenerate(iv_yhat) whatgenerate(iv_what) zhatgenerate(iv_zhat)
    assert !missing(iv_yhat) in 1
    assert !missing(iv_what) in 1
    assert !missing(iv_zhat) in 1

    * Step 2: Refit with user-supplied nuisance
    grf_instrumental_forest y w z x1 x2, gen(iv_n2) ntrees(100) seed(42) ///
        yhatinput(iv_yhat) whatinput(iv_what) zhatinput(iv_zhat) replace
    assert !missing(iv_n2) in 1

    drop iv_n1 iv_n2 iv_yhat iv_what iv_zhat
}
if _rc {
    display as error "FAIL: E3.4 instrumental_forest nuisance inputs"
    local errors = `errors' + 1
}
else {
    display as result "PASS: E3.4 instrumental_forest nuisance inputs"
}

* ---- E3.5: lm_forest with nuisance inputs ----
capture noisily {
    clear
    set obs 300
    set seed 42
    gen x1 = rnormal()
    gen x2 = rnormal()
    gen w1 = rnormal()
    gen y = 2*x1 + 0.5*w1*x1 + rnormal()

    * Step 1: Fit with auto-nuisance, custom generate names
    grf_lm_forest y w1, xvars(x1 x2) gen(lmn1) ntrees(100) seed(42) ///
        yhatgenerate(lm_yhat) whatgenerate(lm_what1)
    assert !missing(lm_yhat) in 1
    assert !missing(lm_what1) in 1

    * Step 2: Refit with user-supplied nuisance
    grf_lm_forest y w1, xvars(x1 x2) gen(lmn2) ntrees(100) seed(42) ///
        yhatinput(lm_yhat) whatinput(lm_what1) replace
    assert !missing(lmn2_1) in 1

    drop lmn1_1 lmn2_1 lm_yhat lm_what1
    capture drop _grf_lm_yhat _grf_lm_what1
}
if _rc {
    display as error "FAIL: E3.5 lm_forest nuisance inputs"
    local errors = `errors' + 1
}
else {
    display as result "PASS: E3.5 lm_forest nuisance inputs"
}

* ---- E3.6: causal_survival_forest with whatinput ----
capture noisily {
    clear
    set obs 300
    set seed 42
    gen x1 = rnormal()
    gen x2 = rnormal()
    gen w = rbinomial(1, normal(0.3*x1))
    gen t = exp(0.5*x1 + w + rnormal())
    gen c = exp(2 + rnormal())
    gen time = min(t, c)
    gen status = (t <= c)

    * Step 1: Compute propensity score externally
    quietly logit w x1 x2
    quietly predict my_what_cs, pr

    * Step 2: Fit with user-supplied propensity
    grf_causal_survival_forest time status w x1 x2, gen(csf_n1) ///
        ntrees(100) seed(42) whatinput(my_what_cs)
    assert !missing(csf_n1) in 1

    drop csf_n1 my_what_cs
}
if _rc {
    display as error "FAIL: E3.6 causal_survival_forest whatinput"
    local errors = `errors' + 1
}
else {
    display as result "PASS: E3.6 causal_survival_forest whatinput"
}


* ============================================================
* E4: New Features
* ============================================================

display as text ""
display as text "============================================================"
display as text "E4: New Features"
display as text "============================================================"

* ---- E4.1: equalize.cluster.weights on regression forest ----
capture noisily {
    clear
    set obs 300
    set seed 42
    * Unequal clusters: cluster 1 has 100, cluster 2 has 200
    gen cluster_id = cond(_n <= 100, 1, 2)
    gen x1 = rnormal()
    gen x2 = rnormal()
    gen y = 2*x1 + x2 + rnormal()

    grf_regression_forest y x1 x2, gen(pred_eq) ntrees(100) seed(42) ///
        cluster(cluster_id) equalizeclusterweights
    assert !missing(pred_eq) in 1
    assert "`e(cluster_var)'" == "cluster_id"

    * Compare with non-equalized
    grf_regression_forest y x1 x2, gen(pred_neq) ntrees(100) seed(42) ///
        cluster(cluster_id) replace
    gen diff_eq = pred_eq - pred_neq
    quietly summarize diff_eq
    * Predictions should differ with equalized weights
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

* ---- E4.2: equalize.cluster.weights error without cluster ----
capture {
    grf_regression_forest y x1 x2, gen(pred_eqerr) ntrees(100) seed(42) ///
        equalizeclusterweights
}
if _rc {
    display as result "PASS: E4.2 error: equalizeclusterweights without cluster()"
}
else {
    display as error "FAIL: E4.2 should error without cluster()"
    local errors = `errors' + 1
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

* ---- E4.4: debiasing.weights on grf_ate ----
capture noisily {
    clear
    set obs 400
    set seed 42
    gen x1 = rnormal()
    gen x2 = rnormal()
    gen w = rbinomial(1, normal(0.3*x1))
    gen y = 2*x1 + x1*w + rnormal()
    gen dbwt = abs(x1) + 0.5

    * Fit causal forest
    grf_causal_forest y w x1 x2, gen(tau_db) ntrees(200) seed(42)

    * ATE without debiasing weights
    grf_ate
    local ate_nodb = r(ate)

    * ATE with debiasing weights
    grf_ate, debiasingweights(dbwt)
    local ate_db = r(ate)

    * Should differ
    assert `ate_nodb' != `ate_db'
    drop tau_db _grf_yhat _grf_what dbwt
}
if _rc {
    display as error "FAIL: E4.4 debiasing.weights grf_ate"
    local errors = `errors' + 1
}
else {
    display as result "PASS: E4.4 debiasing.weights grf_ate"
}

* ---- E4.5: debiasing.weights on grf_best_linear_projection ----
capture noisily {
    clear
    set obs 400
    set seed 42
    gen x1 = rnormal()
    gen x2 = rnormal()
    gen w = rbinomial(1, normal(0.3*x1))
    gen y = 2*x1 + x1*w + rnormal()
    gen dbwt = abs(x1) + 0.5

    grf_causal_forest y w x1 x2, gen(tau_blp) ntrees(200) seed(42)

    * BLP without debiasing weights
    grf_best_linear_projection x1 x2
    matrix b_nodb = e(b)

    * Re-fit causal forest (BLP replaces e() results)
    grf_causal_forest y w x1 x2, gen(tau_blp2) ntrees(200) seed(42) replace

    * BLP with debiasing weights
    grf_best_linear_projection x1 x2, debiasingweights(dbwt)
    matrix b_db = e(b)

    * Coefficients should differ
    assert b_nodb[1,1] != b_db[1,1]
    drop tau_blp tau_blp2 _grf_yhat _grf_what dbwt
}
if _rc {
    display as error "FAIL: E4.5 debiasing.weights grf_best_linear_projection"
    local errors = `errors' + 1
}
else {
    display as result "PASS: E4.5 debiasing.weights grf_best_linear_projection"
}

* ---- E4.6: debiasing.weights on grf_rate ----
capture noisily {
    clear
    set obs 400
    set seed 42
    gen x1 = rnormal()
    gen x2 = rnormal()
    gen w = rbinomial(1, normal(0.3*x1))
    gen y = 2*x1 + x1*w + rnormal()
    gen dbwt = abs(x1) + 0.5

    grf_causal_forest y w x1 x2, gen(tau_rate) ntrees(200) seed(42)

    * RATE without debiasing weights
    grf_rate tau_rate, target(AUTOC) bootstrap(50) seed(42)
    local rate_nodb = r(estimate)

    * RATE with debiasing weights
    grf_rate tau_rate, target(AUTOC) bootstrap(50) seed(42) ///
        debiasingweights(dbwt)
    local rate_db = r(estimate)

    * Should differ
    assert `rate_nodb' != `rate_db'
    drop tau_rate _grf_yhat _grf_what dbwt
}
if _rc {
    display as error "FAIL: E4.6 debiasing.weights grf_rate"
    local errors = `errors' + 1
}
else {
    display as result "PASS: E4.6 debiasing.weights grf_rate"
}

* ---- E4.7: decay.exponent on grf_variable_importance ----
capture noisily {
    clear
    set obs 300
    set seed 42
    gen x1 = rnormal()
    gen x2 = rnormal()
    gen x3 = rnormal()
    gen y = 2*x1 + x2^2 + rnormal()

    * Default decay exponent (2.0)
    grf_variable_importance y x1 x2 x3, ntrees(200) seed(42)
    matrix vi_d2 = r(importance)

    * Different decay exponent
    grf_variable_importance y x1 x2 x3, ntrees(200) seed(42) decayexponent(0.5)
    matrix vi_d05 = r(importance)

    * Importance rankings should differ with different decay
    * (absolute values will differ at minimum)
    local diff = abs(vi_d2[1,1] - vi_d05[1,1]) + ///
                 abs(vi_d2[1,2] - vi_d05[1,2]) + ///
                 abs(vi_d2[1,3] - vi_d05[1,3])
    assert `diff' > 0
}
if _rc {
    display as error "FAIL: E4.7 decay.exponent variable_importance"
    local errors = `errors' + 1
}
else {
    display as result "PASS: E4.7 decay.exponent variable_importance"
}

* ---- E4.8: decay.exponent default is 2.0 ----
capture noisily {
    * With explicit 2.0, results should match default
    grf_variable_importance y x1 x2 x3, ntrees(200) seed(42) decayexponent(2.0)
    matrix vi_explicit = r(importance)

    grf_variable_importance y x1 x2 x3, ntrees(200) seed(42)
    matrix vi_default = r(importance)

    * Should be identical
    forvalues j = 1/3 {
        assert reldif(vi_explicit[1,`j'], vi_default[1,`j']) < 1e-10
    }
}
if _rc {
    display as error "FAIL: E4.8 decay.exponent default=2.0"
    local errors = `errors' + 1
}
else {
    display as result "PASS: E4.8 decay.exponent default=2.0"
}

* ---- E4.9: equalize.cluster.weights on survival forest ----
capture noisily {
    clear
    set obs 300
    set seed 42
    gen cluster_id = ceil(_n / 15)
    gen x1 = rnormal()
    gen x2 = rnormal()
    gen time = exp(x1 + rnormal())
    gen status = (uniform() > 0.3)

    grf_survival_forest time status x1 x2, gen(sv_eq) ntrees(100) seed(42) ///
        cluster(cluster_id) equalizeclusterweights
    assert !missing(sv_eq_s1) in 1
    forvalues j = 1/20 {
        capture drop sv_eq_s`j'
    }
}
if _rc {
    display as error "FAIL: E4.9 equalize.cluster.weights survival forest"
    local errors = `errors' + 1
}
else {
    display as result "PASS: E4.9 equalize.cluster.weights survival forest"
}

* ---- E4.10: equalize.cluster.weights combined with weights ----
capture noisily {
    clear
    set obs 300
    set seed 42
    gen cluster_id = ceil(_n / 10)
    gen x1 = rnormal()
    gen x2 = rnormal()
    gen y = 2*x1 + x2 + rnormal()
    gen user_wt = abs(x1) + 0.5

    * Should combine equalized weights with user weights
    grf_regression_forest y x1 x2, gen(pred_eqwt) ntrees(100) seed(42) ///
        cluster(cluster_id) weights(user_wt) equalizeclusterweights
    assert !missing(pred_eqwt) in 1
    drop pred_eqwt
}
if _rc {
    display as error "FAIL: E4.10 equalize.cluster.weights + weights"
    local errors = `errors' + 1
}
else {
    display as result "PASS: E4.10 equalize.cluster.weights + weights"
}

* ---- E4.11: ATE with target.sample + debiasing.weights ----
capture noisily {
    clear
    set obs 400
    set seed 42
    gen x1 = rnormal()
    gen x2 = rnormal()
    gen w = rbinomial(1, normal(0.3*x1))
    gen y = 2*x1 + x1*w + rnormal()
    gen dbwt = abs(x1) + 0.5

    grf_causal_forest y w x1 x2, gen(tau_ts) ntrees(200) seed(42)

    * ATE overlap + debiasing
    grf_ate, targetsample(overlap) debiasingweights(dbwt)
    assert !missing(r(ate))
    assert !missing(r(se))
    assert r(se) > 0

    * ATE treated + debiasing
    grf_ate, targetsample(treated) debiasingweights(dbwt)
    assert !missing(r(ate))

    drop tau_ts _grf_yhat _grf_what dbwt
}
if _rc {
    display as error "FAIL: E4.11 ATE target.sample + debiasing.weights"
    local errors = `errors' + 1
}
else {
    display as result "PASS: E4.11 ATE target.sample + debiasing.weights"
}


* ============================================================
* E5: TMLE and regression.splitting
* ============================================================

display as text ""
display as text "============================================================"
display as text "E5: TMLE and regression.splitting"
display as text "============================================================"

* ---- Setup: shared causal forest for TMLE tests ----
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

* ---- E5.2: TMLE vs AIPW estimates differ ----
capture noisily {
    grf_ate, method(AIPW)
    local ate_aipw = r(ate)
    grf_ate, method(TMLE)
    local ate_tmle = r(ate)
    * Estimates should differ (TMLE applies bias correction)
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

* ---- E5.5: TMLE + target.sample(overlap) falls back to AIPW ----
capture noisily {
    grf_ate, method(TMLE) targetsample(overlap)
    assert !missing(r(ate))
    assert !missing(r(se))
    * Should fall back to AIPW
    assert "`r(method)'" == "AIPW"
}
if _rc {
    display as error "FAIL: E5.5 TMLE overlap fallback to AIPW"
    local errors = `errors' + 1
}
else {
    display as result "PASS: E5.5 TMLE overlap fallback to AIPW"
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

* ---- E5.8: regression.splitting quantile forest ----
capture noisily {
    clear
    set obs 200
    set seed 42
    gen x1 = rnormal()
    gen x2 = rnormal()
    gen y = 2*x1 + x2^2 + rnormal()

    grf_quantile_forest y x1 x2, gen(qf_rs) ntrees(100) seed(42) ///
        quantiles(0.1 0.5 0.9) regressionsplitting
    assert !missing(qf_rs_q10) in 1
    assert !missing(qf_rs_q50) in 1
    assert !missing(qf_rs_q90) in 1
    assert e(regression_splitting) == 1
    drop qf_rs_q10 qf_rs_q50 qf_rs_q90
}
if _rc {
    display as error "FAIL: E5.8 regression.splitting quantile forest"
    local errors = `errors' + 1
}
else {
    display as result "PASS: E5.8 regression.splitting quantile forest"
}

* ---- E5.9: regression.splitting changes predictions ----
capture noisily {
    clear
    set obs 200
    set seed 42
    gen x1 = rnormal()
    gen x2 = rnormal()
    gen y = 2*x1 + x2^2 + rnormal()

    * Without regression splitting
    grf_quantile_forest y x1 x2, gen(qf_def) ntrees(100) seed(42) ///
        quantiles(0.5)
    assert e(regression_splitting) == 0

    * With regression splitting
    grf_quantile_forest y x1 x2, gen(qf_rs2) ntrees(100) seed(42) ///
        quantiles(0.5) regressionsplitting

    gen diff_rs = qf_def_q50 - qf_rs2_q50
    quietly summarize diff_rs
    * Predictions should differ (different splitting rule)
    assert r(sd) > 0
    drop qf_def_q50 qf_rs2_q50 diff_rs
}
if _rc {
    display as error "FAIL: E5.9 regression.splitting changes predictions"
    local errors = `errors' + 1
}
else {
    display as result "PASS: E5.9 regression.splitting changes predictions"
}


* ============================================================
* Summary
* ============================================================
display as text ""
display as text "============================================================"
display as text "Summary"
display as text "============================================================"

if `errors' > 0 {
    display as error "FAILED: `errors' error(s) in fidelity gap tests"
    exit 1
}
else {
    display as result "ALL FIDELITY GAP TESTS PASSED (`errors' errors)"
}
