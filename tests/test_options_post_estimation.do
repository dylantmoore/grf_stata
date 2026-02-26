* test_options_post_estimation.do -- Comprehensive option tests for post-estimation commands
* Tests: grf_ate, grf_best_linear_projection, grf_variable_importance,
*        grf_test_calibration, grf_rate, grf_predict, grf_tune

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
* grf_ate tests
* ============================================================

* ---- Test 1: grf_ate basic ----
capture noisily {
    grf_causal_forest y w x1-x5, gen(tau_ate) ntrees(100) seed(42)
    grf_ate
    assert !missing(r(ate))
    assert !missing(r(se))
    assert r(se) > 0
    assert !missing(r(ci_lower))
    assert !missing(r(ci_upper))
    assert !missing(r(pvalue))
    assert r(N) == 500
    drop tau_ate _grf_yhat _grf_what
}
if _rc {
    display as error "FAIL: grf_ate basic"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_ate basic"
}

* ---- Test 2: grf_ate with if restriction ----
capture noisily {
    grf_causal_forest y w x1-x5, gen(tau_ate2) ntrees(100) seed(42)
    grf_ate if x1 > 0
    * Save r() results before any other r-class command overwrites them
    local saved_ate = r(ate)
    local saved_N = r(N)
    assert !missing(`saved_ate')
    drop tau_ate2 _grf_yhat _grf_what
}
if _rc {
    display as error "FAIL: grf_ate with if restriction"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_ate with if restriction"
}

* ---- Test 3: grf_ate with in restriction ----
capture noisily {
    grf_causal_forest y w x1-x5, gen(tau_ate3) ntrees(100) seed(42)
    grf_ate in 1/250
    assert r(N) == 250
    assert !missing(r(ate))
    drop tau_ate3 _grf_yhat _grf_what
}
if _rc {
    display as error "FAIL: grf_ate with in restriction"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_ate with in restriction"
}

* ---- Test 4: grf_ate error without prior causal forest ----
capture noisily {
    * Run a regression forest first (not causal)
    grf_regression_forest y x1-x5, gen(rpred) ntrees(100) seed(42)
    capture grf_ate
    assert _rc == 301
    drop rpred
}
if _rc {
    display as error "FAIL: grf_ate error without causal forest"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_ate error without causal forest"
}

* ============================================================
* grf_best_linear_projection tests
* ============================================================

* ---- Test 5: grf_best_linear_projection default (no varlist, uses e(indepvars)) ----
capture noisily {
    grf_causal_forest y w x1-x5, gen(tau_blp1) ntrees(100) seed(42)
    grf_best_linear_projection
    * After regress, e() is from the BLP regression
    assert e(N) > 0
    drop tau_blp1 _grf_yhat _grf_what
}
if _rc {
    display as error "FAIL: grf_best_linear_projection default"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_best_linear_projection default"
}

* ---- Test 6: grf_best_linear_projection with explicit varlist ----
capture noisily {
    grf_causal_forest y w x1-x5, gen(tau_blp2) ntrees(100) seed(42)
    grf_best_linear_projection x1 x2
    assert e(N) > 0
    drop tau_blp2 _grf_yhat _grf_what
}
if _rc {
    display as error "FAIL: grf_best_linear_projection with varlist"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_best_linear_projection with varlist"
}

* ---- Test 7: grf_best_linear_projection with if ----
capture noisily {
    grf_causal_forest y w x1-x5, gen(tau_blp3) ntrees(100) seed(42)
    grf_best_linear_projection x1 x2 if x1 > 0
    assert e(N) > 0
    drop tau_blp3 _grf_yhat _grf_what
}
if _rc {
    display as error "FAIL: grf_best_linear_projection with if"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_best_linear_projection with if"
}

* ---- Test 8: grf_best_linear_projection with in ----
capture noisily {
    grf_causal_forest y w x1-x5, gen(tau_blp4) ntrees(100) seed(42)
    grf_best_linear_projection x1 in 1/250
    assert e(N) == 250
    drop tau_blp4 _grf_yhat _grf_what
}
if _rc {
    display as error "FAIL: grf_best_linear_projection with in"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_best_linear_projection with in"
}

* ============================================================
* grf_variable_importance tests
* ============================================================

* ---- Test 9: grf_variable_importance basic ----
capture noisily {
    grf_variable_importance y x1-x5, ntrees(100) seed(42)
    assert r(N) == 500
    matrix vi = r(importance)
    assert colsof(vi) == 5
    * All importance values should be non-negative
    forvalues j = 1/5 {
        assert vi[1,`j'] >= 0
    }
    * Sum should be close to 1 (importance is normalized)
    local vi_sum = 0
    forvalues j = 1/5 {
        local vi_sum = `vi_sum' + vi[1,`j']
    }
    assert reldif(`vi_sum', 1.0) < 0.01
}
if _rc {
    display as error "FAIL: grf_variable_importance basic"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_variable_importance basic"
}

* ---- Test 10: grf_variable_importance maxdepth(6) ----
capture noisily {
    grf_variable_importance y x1-x5, ntrees(100) seed(42) maxdepth(6)
    assert r(N) == 500
    matrix vi6 = r(importance)
    assert colsof(vi6) == 5
}
if _rc {
    display as error "FAIL: grf_variable_importance maxdepth(6)"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_variable_importance maxdepth(6)"
}

* ---- Test 11: grf_variable_importance with non-default ntrees ----
capture noisily {
    grf_variable_importance y x1-x5, ntrees(50) seed(123)
    assert r(N) == 500
    matrix vi_50 = r(importance)
    assert colsof(vi_50) == 5
}
if _rc {
    display as error "FAIL: grf_variable_importance ntrees(50)"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_variable_importance ntrees(50)"
}

* ---- Test 12: grf_variable_importance with if ----
capture noisily {
    grf_variable_importance y x1-x5 if x1 > 0, ntrees(100) seed(42)
    * Save r() results immediately before any other r-class command
    local saved_vi_N = r(N)
    matrix vi_if = r(importance)
    assert `saved_vi_N' > 0
    assert `saved_vi_N' < 500
    assert colsof(vi_if) == 5
}
if _rc {
    display as error "FAIL: grf_variable_importance with if"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_variable_importance with if"
}

* ---- Test 13: grf_variable_importance with in ----
capture noisily {
    grf_variable_importance y x1-x5 in 1/250, ntrees(100) seed(42)
    assert r(N) == 250
}
if _rc {
    display as error "FAIL: grf_variable_importance with in"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_variable_importance with in"
}

* ============================================================
* grf_test_calibration tests
* ============================================================

* ---- Test 14: grf_test_calibration basic ----
capture noisily {
    grf_causal_forest y w x1-x5, gen(tau_cal1) ntrees(100) seed(42)
    grf_test_calibration
    assert !missing(r(b_mean))
    assert !missing(r(se_mean))
    assert !missing(r(t_mean))
    assert !missing(r(p_mean))
    assert !missing(r(b_diff))
    assert !missing(r(se_diff))
    assert !missing(r(t_diff))
    assert !missing(r(p_diff))
    assert r(N) > 0
    assert r(se_mean) > 0
    assert r(se_diff) > 0
    drop tau_cal1 _grf_yhat _grf_what
}
if _rc {
    display as error "FAIL: grf_test_calibration basic"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_test_calibration basic"
}

* ---- Test 15: grf_test_calibration with if ----
capture noisily {
    grf_causal_forest y w x1-x5, gen(tau_cal2) ntrees(100) seed(42)
    grf_test_calibration if x1 > 0
    assert !missing(r(b_mean))
    assert r(N) > 0
    drop tau_cal2 _grf_yhat _grf_what
}
if _rc {
    display as error "FAIL: grf_test_calibration with if"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_test_calibration with if"
}

* ---- Test 16: grf_test_calibration with in ----
capture noisily {
    grf_causal_forest y w x1-x5, gen(tau_cal3) ntrees(100) seed(42)
    grf_test_calibration in 1/400
    assert r(N) == 400
    drop tau_cal3 _grf_yhat _grf_what
}
if _rc {
    display as error "FAIL: grf_test_calibration with in"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_test_calibration with in"
}

* ============================================================
* grf_rate tests
* ============================================================

* ---- Test 17: grf_rate default (AUTOC) ----
capture noisily {
    grf_causal_forest y w x1-x5, gen(tau_rate1) ntrees(100) seed(42)
    grf_rate tau_rate1, bootstrap(50) seed(42)
    assert !missing(r(estimate))
    assert !missing(r(std_err))
    assert r(std_err) > 0
    assert !missing(r(z_stat))
    assert !missing(r(p_value))
    assert r(n) > 0
    assert "`r(target)'" == "AUTOC"
    drop tau_rate1 _grf_yhat _grf_what
}
if _rc {
    display as error "FAIL: grf_rate AUTOC"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_rate AUTOC"
}

* ---- Test 18: grf_rate QINI ----
capture noisily {
    grf_causal_forest y w x1-x5, gen(tau_rate2) ntrees(100) seed(42)
    grf_rate tau_rate2, target(QINI) bootstrap(50) seed(42)
    assert !missing(r(estimate))
    assert "`r(target)'" == "QINI"
    drop tau_rate2 _grf_yhat _grf_what
}
if _rc {
    display as error "FAIL: grf_rate QINI"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_rate QINI"
}

* ---- Test 19: grf_rate custom quantiles ----
capture noisily {
    grf_causal_forest y w x1-x5, gen(tau_rate3) ntrees(100) seed(42)
    grf_rate tau_rate3, quantiles(0.25 0.5 0.75 1.0) bootstrap(50) seed(42)
    assert !missing(r(estimate))
    drop tau_rate3 _grf_yhat _grf_what
}
if _rc {
    display as error "FAIL: grf_rate custom quantiles"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_rate custom quantiles"
}

* ---- Test 20: grf_rate with if ----
capture noisily {
    grf_causal_forest y w x1-x5, gen(tau_rate4) ntrees(100) seed(42)
    grf_rate tau_rate4 if x1 > 0, bootstrap(50) seed(42)
    assert !missing(r(estimate))
    assert r(n) < 500
    drop tau_rate4 _grf_yhat _grf_what
}
if _rc {
    display as error "FAIL: grf_rate with if"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_rate with if"
}

* ---- Test 21: grf_rate with catevar ----
capture noisily {
    grf_causal_forest y w x1-x5, gen(tau_rate5) ntrees(100) seed(42)
    * Use the forest's predictions directly via catevar
    grf_rate tau_rate5, catevar(tau_rate5) bootstrap(50) seed(42)
    assert !missing(r(estimate))
    drop tau_rate5 _grf_yhat _grf_what
}
if _rc {
    display as error "FAIL: grf_rate with catevar"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_rate with catevar"
}

* ============================================================
* grf_predict tests
* ============================================================

* ---- Test 22: grf_predict regression forest with numthreads ----
capture noisily {
    * Fit on first 400 obs
    grf_regression_forest y x1-x5 in 1/400, gen(reg_pred) ntrees(100) seed(42)
    local n_train = e(N)

    * Append 100 more obs (simulate test data)
    preserve
    set obs 600
    forvalues j = 1/5 {
        quietly replace x`j' = rnormal() if _n > 500
    }
    quietly replace y = rnormal() if _n > 500

    * Predict on test obs using numthreads
    grf_predict, gen(pred_test) numthreads(2)

    * Save r() results immediately before count overwrites them
    local saved_forest_type "`r(forest_type)'"
    quietly count if !missing(pred_test) & _n > `n_train'
    assert r(N) > 0
    assert "`saved_forest_type'" == "regression"
    restore
    drop reg_pred
}
if _rc {
    display as error "FAIL: grf_predict with numthreads(2)"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_predict with numthreads(2)"
}

* ---- Test 23: grf_predict causal forest ----
capture noisily {
    * Fit causal forest on first 400 obs
    grf_causal_forest y w x1-x5 in 1/400, gen(tau_pred) ntrees(100) seed(42)
    local n_train = e(N)

    * Append test data
    preserve
    set obs 600
    forvalues j = 1/5 {
        quietly replace x`j' = rnormal() if _n > 500
    }
    quietly replace y = rnormal() if _n > 500
    quietly replace w = (rnormal() > 0) if _n > 500

    grf_predict, gen(cate_test)
    * Save r() results immediately before count overwrites them
    local saved_forest_type "`r(forest_type)'"
    quietly count if !missing(cate_test) & _n > `n_train'
    assert r(N) > 0
    assert "`saved_forest_type'" == "causal"
    restore
    drop tau_pred _grf_yhat _grf_what
}
if _rc {
    display as error "FAIL: grf_predict causal forest"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_predict causal forest"
}

* ---- Test 24: grf_predict replace ----
capture noisily {
    grf_regression_forest y x1-x5 in 1/400, gen(reg_pred2) ntrees(100) seed(42)
    preserve
    set obs 600
    forvalues j = 1/5 {
        quietly replace x`j' = rnormal() if _n > 500
    }
    quietly replace y = rnormal() if _n > 500

    grf_predict, gen(pred_r)
    grf_predict, gen(pred_r) replace
    quietly count if !missing(pred_r)
    assert r(N) > 0
    restore
    drop reg_pred2
}
if _rc {
    display as error "FAIL: grf_predict replace"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_predict replace"
}

* ============================================================
* grf_tune tests
* ============================================================

* ---- Test 25: grf_tune foresttype("regression") ----
capture noisily {
    grf_tune y x1-x5, foresttype("regression") numreps(10) tunetrees(50) seed(42)
    assert !missing(r(best_mtry))
    assert r(best_mtry) >= 1
    assert !missing(r(best_min_node_size))
    assert r(best_min_node_size) >= 1
    assert !missing(r(best_sample_fraction))
    assert r(best_sample_fraction) > 0 & r(best_sample_fraction) < 1
    assert !missing(r(best_honesty_fraction))
    assert !missing(r(best_alpha))
    assert !missing(r(best_imbalance_penalty))
    assert !missing(r(best_mse))
    assert r(best_mse) > 0
    assert r(n_reps) == 10
    assert r(N) == 500
    assert "`r(forest_type)'" == "regression"
}
if _rc {
    display as error "FAIL: grf_tune regression"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_tune regression"
}

* ---- Test 26: grf_tune foresttype("causal") ----
capture noisily {
    grf_tune y w x1-x5, foresttype("causal") numreps(10) tunetrees(50) seed(42)
    assert !missing(r(best_mtry))
    assert !missing(r(best_mse))
    assert "`r(forest_type)'" == "causal"
}
if _rc {
    display as error "FAIL: grf_tune causal"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_tune causal"
}

* ---- Test 27: grf_tune foresttype("quantile") ----
capture noisily {
    grf_tune y x1-x5, foresttype("quantile") numreps(10) tunetrees(50) seed(42)
    assert !missing(r(best_mtry))
    assert !missing(r(best_mse))
    assert "`r(forest_type)'" == "quantile"
}
if _rc {
    display as error "FAIL: grf_tune quantile"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_tune quantile"
}

* ---- Test 28: grf_tune nohonesty ----
capture noisily {
    grf_tune y x1-x5, foresttype("regression") numreps(10) tunetrees(50) seed(42) nohonesty
    assert !missing(r(best_mtry))
    assert !missing(r(best_mse))
}
if _rc {
    display as error "FAIL: grf_tune nohonesty"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_tune nohonesty"
}

* ---- Test 29: grf_tune nohonestyprune ----
capture noisily {
    grf_tune y x1-x5, foresttype("regression") numreps(10) tunetrees(50) seed(42) nohonestyprune
    assert !missing(r(best_mtry))
    assert !missing(r(best_mse))
}
if _rc {
    display as error "FAIL: grf_tune nohonestyprune"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_tune nohonestyprune"
}

* ---- Test 30: grf_tune nostabilizesplits ----
capture noisily {
    grf_tune y w x1-x5, foresttype("causal") numreps(10) tunetrees(50) seed(42) nostabilizesplits
    assert !missing(r(best_mtry))
    assert !missing(r(best_mse))
}
if _rc {
    display as error "FAIL: grf_tune nostabilizesplits"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_tune nostabilizesplits"
}

* ---- Test 31: grf_tune numthreads(2) ----
capture noisily {
    grf_tune y x1-x5, foresttype("regression") numreps(10) tunetrees(50) seed(42) numthreads(2)
    assert !missing(r(best_mtry))
    assert !missing(r(best_mse))
}
if _rc {
    display as error "FAIL: grf_tune numthreads(2)"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_tune numthreads(2)"
}

* ---- Test 32: grf_tune with if ----
capture noisily {
    grf_tune y x1-x5 if x1 > 0, foresttype("regression") numreps(10) tunetrees(50) seed(42)
    * Save r() results immediately before any other r-class command
    local saved_tune_N = r(N)
    local saved_best_mse = r(best_mse)
    assert `saved_tune_N' > 0
    assert `saved_tune_N' < 500
    assert !missing(`saved_best_mse')
}
if _rc {
    display as error "FAIL: grf_tune with if"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_tune with if"
}

* ---- Test 33: grf_tune with in ----
capture noisily {
    grf_tune y x1-x5 in 1/300, foresttype("regression") numreps(10) tunetrees(50) seed(42)
    assert r(N) == 300
    assert !missing(r(best_mse))
}
if _rc {
    display as error "FAIL: grf_tune with in"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_tune with in"
}

* ============================================================
* Summary
* ============================================================
display ""
if `errors' > 0 {
    display as error "FAILED: `errors' errors in post-estimation option tests"
    exit 1
}
else {
    display as result "ALL POST-ESTIMATION OPTION TESTS PASSED"
}
