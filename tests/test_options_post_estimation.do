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

* ---- Test 21b: grf_rate with in restriction ----
capture noisily {
    grf_causal_forest y w x1-x5, gen(tau_rate_in) ntrees(100) seed(42)
    grf_rate tau_rate_in in 1/400, bootstrap(50) seed(42)
    local saved_N = r(N)
    assert !missing(r(estimate))
    assert `saved_N' <= 400
    assert `saved_N' > 0
    drop tau_rate_in _grf_yhat _grf_what
}
if _rc {
    display as error "FAIL: grf_rate with in restriction"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_rate with in restriction"
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
* grf_ate target.sample tests
* ============================================================

* ---- Test 34: grf_ate targetsample(all) matches default ----
capture noisily {
    grf_causal_forest y w x1-x5, gen(tau_ts1) ntrees(100) seed(42)
    grf_ate
    local ate_default = r(ate)
    local se_default = r(se)
    drop tau_ts1 _grf_yhat _grf_what

    grf_causal_forest y w x1-x5, gen(tau_ts1b) ntrees(100) seed(42)
    grf_ate, targetsample(all)
    local ate_all = r(ate)
    local se_all = r(se)
    drop tau_ts1b _grf_yhat _grf_what

    assert reldif(`ate_default', `ate_all') < 1e-10
    assert reldif(`se_default', `se_all') < 1e-10
}
if _rc {
    display as error "FAIL: grf_ate targetsample(all) matches default"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_ate targetsample(all) matches default"
}

* ---- Test 35: grf_ate targetsample(all) returns valid results ----
capture noisily {
    grf_causal_forest y w x1-x5, gen(tau_ts2) ntrees(100) seed(42)
    grf_ate, targetsample(all)
    assert !missing(r(ate))
    assert !missing(r(se))
    assert r(se) > 0
    assert !missing(r(ci_lower))
    assert !missing(r(ci_upper))
    assert r(ci_lower) < r(ci_upper)
    assert r(N) == 500
    drop tau_ts2 _grf_yhat _grf_what
}
if _rc {
    display as error "FAIL: grf_ate targetsample(all) valid results"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_ate targetsample(all) valid results"
}

* ---- Test 36: grf_ate targetsample(treated) returns valid results ----
capture noisily {
    grf_causal_forest y w x1-x5, gen(tau_ts3) ntrees(100) seed(42)
    grf_ate, targetsample(treated)
    assert !missing(r(ate))
    assert !missing(r(se))
    assert r(se) > 0
    assert !missing(r(ci_lower))
    assert !missing(r(ci_upper))
    assert r(ci_lower) < r(ci_upper)
    assert r(N) == 500
    assert "`r(target_sample)'" == "treated"
    drop tau_ts3 _grf_yhat _grf_what
}
if _rc {
    display as error "FAIL: grf_ate targetsample(treated) valid results"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_ate targetsample(treated) valid results"
}

* ---- Test 37: grf_ate targetsample(control) returns valid results ----
capture noisily {
    grf_causal_forest y w x1-x5, gen(tau_ts4) ntrees(100) seed(42)
    grf_ate, targetsample(control)
    assert !missing(r(ate))
    assert !missing(r(se))
    assert r(se) > 0
    assert !missing(r(ci_lower))
    assert !missing(r(ci_upper))
    assert r(ci_lower) < r(ci_upper)
    assert r(N) == 500
    assert "`r(target_sample)'" == "control"
    drop tau_ts4 _grf_yhat _grf_what
}
if _rc {
    display as error "FAIL: grf_ate targetsample(control) valid results"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_ate targetsample(control) valid results"
}

* ---- Test 38: grf_ate targetsample(overlap) returns valid results ----
capture noisily {
    grf_causal_forest y w x1-x5, gen(tau_ts5) ntrees(100) seed(42)
    grf_ate, targetsample(overlap)
    assert !missing(r(ate))
    assert !missing(r(se))
    assert r(se) > 0
    assert !missing(r(ci_lower))
    assert !missing(r(ci_upper))
    assert r(ci_lower) < r(ci_upper)
    assert r(N) == 500
    assert "`r(target_sample)'" == "overlap"
    drop tau_ts5 _grf_yhat _grf_what
}
if _rc {
    display as error "FAIL: grf_ate targetsample(overlap) valid results"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_ate targetsample(overlap) valid results"
}

* ---- Test 39: grf_ate targetsample(invalid) errors ----
capture noisily {
    grf_causal_forest y w x1-x5, gen(tau_ts6) ntrees(100) seed(42)
    capture grf_ate, targetsample(invalid)
    assert _rc == 198
    drop tau_ts6 _grf_yhat _grf_what
}
if _rc {
    display as error "FAIL: grf_ate targetsample(invalid) errors"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_ate targetsample(invalid) errors"
}

* ---- Test 40: grf_ate targetsample(overlap) differs from all ----
capture noisily {
    grf_causal_forest y w x1-x5, gen(tau_ts7) ntrees(100) seed(42)
    grf_ate, targetsample(all)
    local ate_all = r(ate)
    drop tau_ts7 _grf_yhat _grf_what

    grf_causal_forest y w x1-x5, gen(tau_ts7b) ntrees(100) seed(42)
    grf_ate, targetsample(overlap)
    local ate_overlap = r(ate)
    drop tau_ts7b _grf_yhat _grf_what

    * Overlap weighting should produce a different ATE
    assert reldif(`ate_all', `ate_overlap') > 1e-6
}
if _rc {
    display as error "FAIL: grf_ate targetsample(overlap) differs from all"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_ate targetsample(overlap) differs from all"
}

* ---- Test 41: grf_ate treated vs control ATEs differ ----
capture noisily {
    grf_causal_forest y w x1-x5, gen(tau_ts8) ntrees(100) seed(42)
    grf_ate, targetsample(treated)
    local ate_treated = r(ate)
    drop tau_ts8 _grf_yhat _grf_what

    grf_causal_forest y w x1-x5, gen(tau_ts8b) ntrees(100) seed(42)
    grf_ate, targetsample(control)
    local ate_control = r(ate)
    drop tau_ts8b _grf_yhat _grf_what

    * With heterogeneous treatment effects, treated/control ATEs should differ
    assert reldif(`ate_treated', `ate_control') > 1e-6
}
if _rc {
    display as error "FAIL: grf_ate treated vs control ATEs differ"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_ate treated vs control ATEs differ"
}

* ---- Test 42: grf_ate targetsample(treated) with if restriction ----
capture noisily {
    grf_causal_forest y w x1-x5, gen(tau_ts9) ntrees(100) seed(42)
    grf_ate if x1 > 0, targetsample(treated)
    local saved_ate = r(ate)
    local saved_N = r(N)
    assert !missing(`saved_ate')
    assert `saved_N' > 0
    assert `saved_N' < 500
    drop tau_ts9 _grf_yhat _grf_what
}
if _rc {
    display as error "FAIL: grf_ate targetsample(treated) with if"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_ate targetsample(treated) with if"
}

* ---- Test 43: grf_ate targetsample(control) with in restriction ----
capture noisily {
    grf_causal_forest y w x1-x5, gen(tau_ts10) ntrees(100) seed(42)
    grf_ate in 1/300, targetsample(control)
    assert !missing(r(ate))
    assert r(N) == 300
    drop tau_ts10 _grf_yhat _grf_what
}
if _rc {
    display as error "FAIL: grf_ate targetsample(control) with in"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_ate targetsample(control) with in"
}

* ---- Test 44: grf_ate targetsample(overlap) SE is positive ----
capture noisily {
    grf_causal_forest y w x1-x5, gen(tau_ts11) ntrees(100) seed(42)
    grf_ate, targetsample(overlap)
    assert r(se) > 0
    assert r(pvalue) >= 0 & r(pvalue) <= 1
    drop tau_ts11 _grf_yhat _grf_what
}
if _rc {
    display as error "FAIL: grf_ate targetsample(overlap) SE positive"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_ate targetsample(overlap) SE positive"
}

* ---- Test 45: grf_ate all four target.sample SEs are positive ----
capture noisily {
    grf_causal_forest y w x1-x5, gen(tau_ts12a) ntrees(100) seed(42)
    grf_ate, targetsample(all)
    local se_all = r(se)
    drop tau_ts12a _grf_yhat _grf_what

    grf_causal_forest y w x1-x5, gen(tau_ts12b) ntrees(100) seed(42)
    grf_ate, targetsample(treated)
    local se_treated = r(se)
    drop tau_ts12b _grf_yhat _grf_what

    grf_causal_forest y w x1-x5, gen(tau_ts12c) ntrees(100) seed(42)
    grf_ate, targetsample(control)
    local se_control = r(se)
    drop tau_ts12c _grf_yhat _grf_what

    grf_causal_forest y w x1-x5, gen(tau_ts12d) ntrees(100) seed(42)
    grf_ate, targetsample(overlap)
    local se_overlap = r(se)
    drop tau_ts12d _grf_yhat _grf_what

    assert `se_all' > 0
    assert `se_treated' > 0
    assert `se_control' > 0
    assert `se_overlap' > 0
}
if _rc {
    display as error "FAIL: grf_ate all target.sample SEs positive"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_ate all target.sample SEs positive"
}

* ---- Test 46: grf_ate targetsample(all) CI contains ATE ----
capture noisily {
    grf_causal_forest y w x1-x5, gen(tau_ts13) ntrees(100) seed(42)
    grf_ate, targetsample(all)
    assert r(ci_lower) <= r(ate)
    assert r(ci_upper) >= r(ate)
    drop tau_ts13 _grf_yhat _grf_what
}
if _rc {
    display as error "FAIL: grf_ate targetsample(all) CI contains ATE"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_ate targetsample(all) CI contains ATE"
}

* ---- Test 47: grf_ate targetsample(treated) CI contains ATE ----
capture noisily {
    grf_causal_forest y w x1-x5, gen(tau_ts14) ntrees(100) seed(42)
    grf_ate, targetsample(treated)
    assert r(ci_lower) <= r(ate)
    assert r(ci_upper) >= r(ate)
    drop tau_ts14 _grf_yhat _grf_what
}
if _rc {
    display as error "FAIL: grf_ate targetsample(treated) CI contains ATE"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_ate targetsample(treated) CI contains ATE"
}

* ---- Test 48: grf_ate targetsample(control) CI contains ATE ----
capture noisily {
    grf_causal_forest y w x1-x5, gen(tau_ts15) ntrees(100) seed(42)
    grf_ate, targetsample(control)
    assert r(ci_lower) <= r(ate)
    assert r(ci_upper) >= r(ate)
    drop tau_ts15 _grf_yhat _grf_what
}
if _rc {
    display as error "FAIL: grf_ate targetsample(control) CI contains ATE"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_ate targetsample(control) CI contains ATE"
}

* ---- Test 49: grf_ate targetsample(overlap) CI contains ATE ----
capture noisily {
    grf_causal_forest y w x1-x5, gen(tau_ts16) ntrees(100) seed(42)
    grf_ate, targetsample(overlap)
    assert r(ci_lower) <= r(ate)
    assert r(ci_upper) >= r(ate)
    drop tau_ts16 _grf_yhat _grf_what
}
if _rc {
    display as error "FAIL: grf_ate targetsample(overlap) CI contains ATE"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_ate targetsample(overlap) CI contains ATE"
}

* ---- Test 50: grf_ate target_sample returned macro ----
capture noisily {
    grf_causal_forest y w x1-x5, gen(tau_ts17) ntrees(100) seed(42)
    grf_ate, targetsample(treated)
    assert "`r(target_sample)'" == "treated"
    drop tau_ts17 _grf_yhat _grf_what

    grf_causal_forest y w x1-x5, gen(tau_ts17b) ntrees(100) seed(42)
    grf_ate, targetsample(control)
    assert "`r(target_sample)'" == "control"
    drop tau_ts17b _grf_yhat _grf_what

    grf_causal_forest y w x1-x5, gen(tau_ts17c) ntrees(100) seed(42)
    grf_ate, targetsample(overlap)
    assert "`r(target_sample)'" == "overlap"
    drop tau_ts17c _grf_yhat _grf_what

    grf_causal_forest y w x1-x5, gen(tau_ts17d) ntrees(100) seed(42)
    grf_ate, targetsample(all)
    assert "`r(target_sample)'" == "all"
    drop tau_ts17d _grf_yhat _grf_what
}
if _rc {
    display as error "FAIL: grf_ate target_sample returned macro"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_ate target_sample returned macro"
}

* ---- Test 51: grf_ate targetsample(treated) pvalue in [0,1] ----
capture noisily {
    grf_causal_forest y w x1-x5, gen(tau_ts18) ntrees(100) seed(42)
    grf_ate, targetsample(treated)
    assert r(pvalue) >= 0 & r(pvalue) <= 1
    drop tau_ts18 _grf_yhat _grf_what
}
if _rc {
    display as error "FAIL: grf_ate targetsample(treated) pvalue in [0,1]"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_ate targetsample(treated) pvalue in [0,1]"
}

* ---- Test 52: grf_ate targetsample(control) pvalue in [0,1] ----
capture noisily {
    grf_causal_forest y w x1-x5, gen(tau_ts19) ntrees(100) seed(42)
    grf_ate, targetsample(control)
    assert r(pvalue) >= 0 & r(pvalue) <= 1
    drop tau_ts19 _grf_yhat _grf_what
}
if _rc {
    display as error "FAIL: grf_ate targetsample(control) pvalue in [0,1]"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_ate targetsample(control) pvalue in [0,1]"
}

* ---- Test 53: grf_ate default (no option) returns target_sample=="all" ----
capture noisily {
    grf_causal_forest y w x1-x5, gen(tau_ts20) ntrees(100) seed(42)
    grf_ate
    assert "`r(target_sample)'" == "all"
    drop tau_ts20 _grf_yhat _grf_what
}
if _rc {
    display as error "FAIL: grf_ate default returns target_sample all"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_ate default returns target_sample all"
}

* ============================================================
* grf_best_linear_projection vcov.type tests
* ============================================================

* ---- Test 54: grf_best_linear_projection vcovtype(HC0) ----
capture noisily {
    grf_causal_forest y w x1-x5, gen(tau_vc1) ntrees(100) seed(42)
    grf_best_linear_projection x1 x2, vcovtype(HC0)
    assert e(N) > 0
    assert "`e(vcov_type)'" == "HC0"
    drop tau_vc1 _grf_yhat _grf_what
}
if _rc {
    display as error "FAIL: grf_best_linear_projection vcovtype(HC0)"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_best_linear_projection vcovtype(HC0)"
}

* ---- Test 55: grf_best_linear_projection vcovtype(HC1) ----
capture noisily {
    grf_causal_forest y w x1-x5, gen(tau_vc2) ntrees(100) seed(42)
    grf_best_linear_projection x1 x2, vcovtype(HC1)
    assert e(N) > 0
    assert "`e(vcov_type)'" == "HC1"
    drop tau_vc2 _grf_yhat _grf_what
}
if _rc {
    display as error "FAIL: grf_best_linear_projection vcovtype(HC1)"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_best_linear_projection vcovtype(HC1)"
}

* ---- Test 56: grf_best_linear_projection vcovtype(HC2) ----
capture noisily {
    grf_causal_forest y w x1-x5, gen(tau_vc3) ntrees(100) seed(42)
    grf_best_linear_projection x1 x2, vcovtype(HC2)
    assert e(N) > 0
    assert "`e(vcov_type)'" == "HC2"
    drop tau_vc3 _grf_yhat _grf_what
}
if _rc {
    display as error "FAIL: grf_best_linear_projection vcovtype(HC2)"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_best_linear_projection vcovtype(HC2)"
}

* ---- Test 57: grf_best_linear_projection vcovtype(HC3) ----
capture noisily {
    grf_causal_forest y w x1-x5, gen(tau_vc4) ntrees(100) seed(42)
    grf_best_linear_projection x1 x2, vcovtype(HC3)
    assert e(N) > 0
    assert "`e(vcov_type)'" == "HC3"
    drop tau_vc4 _grf_yhat _grf_what
}
if _rc {
    display as error "FAIL: grf_best_linear_projection vcovtype(HC3)"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_best_linear_projection vcovtype(HC3)"
}

* ---- Test 58: grf_best_linear_projection default vcov is HC3 ----
capture noisily {
    grf_causal_forest y w x1-x5, gen(tau_vc5) ntrees(100) seed(42)
    grf_best_linear_projection x1 x2
    assert "`e(vcov_type)'" == "HC3"
    drop tau_vc5 _grf_yhat _grf_what
}
if _rc {
    display as error "FAIL: grf_best_linear_projection default vcov is HC3"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_best_linear_projection default vcov is HC3"
}

* ---- Test 59: grf_best_linear_projection HC0 vs HC3 SEs differ ----
capture noisily {
    grf_causal_forest y w x1-x5, gen(tau_vc6a) ntrees(100) seed(42)
    grf_best_linear_projection x1 x2, vcovtype(HC0)
    local se_hc0 = _se[x1]
    drop tau_vc6a _grf_yhat _grf_what

    grf_causal_forest y w x1-x5, gen(tau_vc6b) ntrees(100) seed(42)
    grf_best_linear_projection x1 x2, vcovtype(HC3)
    local se_hc3 = _se[x1]
    drop tau_vc6b _grf_yhat _grf_what

    * HC0 and HC3 should produce different SEs (HC3 is typically larger)
    assert reldif(`se_hc0', `se_hc3') > 1e-6
}
if _rc {
    display as error "FAIL: grf_best_linear_projection HC0 vs HC3 SEs differ"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_best_linear_projection HC0 vs HC3 SEs differ"
}

* ---- Test 60: grf_best_linear_projection HC1 vs HC2 SEs differ ----
capture noisily {
    grf_causal_forest y w x1-x5, gen(tau_vc7a) ntrees(100) seed(42)
    grf_best_linear_projection x1 x2, vcovtype(HC1)
    local se_hc1 = _se[x1]
    drop tau_vc7a _grf_yhat _grf_what

    grf_causal_forest y w x1-x5, gen(tau_vc7b) ntrees(100) seed(42)
    grf_best_linear_projection x1 x2, vcovtype(HC2)
    local se_hc2 = _se[x1]
    drop tau_vc7b _grf_yhat _grf_what

    * HC1 and HC2 should produce different SEs
    assert reldif(`se_hc1', `se_hc2') > 1e-6
}
if _rc {
    display as error "FAIL: grf_best_linear_projection HC1 vs HC2 SEs differ"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_best_linear_projection HC1 vs HC2 SEs differ"
}

* ---- Test 61: grf_best_linear_projection all HC SEs are positive ----
capture noisily {
    grf_causal_forest y w x1-x5, gen(tau_vc8a) ntrees(100) seed(42)
    grf_best_linear_projection x1, vcovtype(HC0)
    local se_hc0 = _se[x1]
    assert `se_hc0' > 0
    drop tau_vc8a _grf_yhat _grf_what

    grf_causal_forest y w x1-x5, gen(tau_vc8b) ntrees(100) seed(42)
    grf_best_linear_projection x1, vcovtype(HC1)
    local se_hc1 = _se[x1]
    assert `se_hc1' > 0
    drop tau_vc8b _grf_yhat _grf_what

    grf_causal_forest y w x1-x5, gen(tau_vc8c) ntrees(100) seed(42)
    grf_best_linear_projection x1, vcovtype(HC2)
    local se_hc2 = _se[x1]
    assert `se_hc2' > 0
    drop tau_vc8c _grf_yhat _grf_what

    grf_causal_forest y w x1-x5, gen(tau_vc8d) ntrees(100) seed(42)
    grf_best_linear_projection x1, vcovtype(HC3)
    local se_hc3 = _se[x1]
    assert `se_hc3' > 0
    drop tau_vc8d _grf_yhat _grf_what
}
if _rc {
    display as error "FAIL: grf_best_linear_projection all HC SEs positive"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_best_linear_projection all HC SEs positive"
}

* ---- Test 62: grf_best_linear_projection vcovtype(HC4) errors ----
capture noisily {
    grf_causal_forest y w x1-x5, gen(tau_vc9) ntrees(100) seed(42)
    capture grf_best_linear_projection x1, vcovtype(HC4)
    assert _rc == 198
    drop tau_vc9 _grf_yhat _grf_what
}
if _rc {
    display as error "FAIL: grf_best_linear_projection vcovtype(HC4) errors"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_best_linear_projection vcovtype(HC4) errors"
}

* ---- Test 63: grf_best_linear_projection targetsample(overlap) ----
capture noisily {
    grf_causal_forest y w x1-x5, gen(tau_vc10) ntrees(100) seed(42)
    grf_best_linear_projection x1 x2, targetsample(overlap)
    assert e(N) > 0
    assert "`e(target_sample)'" == "overlap"
    drop tau_vc10 _grf_yhat _grf_what
}
if _rc {
    display as error "FAIL: grf_best_linear_projection targetsample(overlap)"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_best_linear_projection targetsample(overlap)"
}

* ---- Test 64: grf_best_linear_projection targetsample(overlap) differs from all ----
capture noisily {
    grf_causal_forest y w x1-x5, gen(tau_vc11a) ntrees(100) seed(42)
    grf_best_linear_projection x1 x2, targetsample(all)
    local coef_all = _b[x1]
    drop tau_vc11a _grf_yhat _grf_what

    grf_causal_forest y w x1-x5, gen(tau_vc11b) ntrees(100) seed(42)
    grf_best_linear_projection x1 x2, targetsample(overlap)
    local coef_overlap = _b[x1]
    drop tau_vc11b _grf_yhat _grf_what

    * Overlap vs all should give different coefficients
    assert reldif(`coef_all', `coef_overlap') > 1e-6
}
if _rc {
    display as error "FAIL: grf_best_linear_projection overlap vs all coefficients differ"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_best_linear_projection overlap vs all coefficients differ"
}

* ---- Test 65: grf_best_linear_projection vcovtype with targetsample(overlap) ----
capture noisily {
    grf_causal_forest y w x1-x5, gen(tau_vc12) ntrees(100) seed(42)
    grf_best_linear_projection x1, vcovtype(HC2) targetsample(overlap)
    assert e(N) > 0
    assert "`e(vcov_type)'" == "HC2"
    assert "`e(target_sample)'" == "overlap"
    local se_val = _se[x1]
    assert `se_val' > 0
    drop tau_vc12 _grf_yhat _grf_what
}
if _rc {
    display as error "FAIL: grf_best_linear_projection vcovtype+targetsample combined"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_best_linear_projection vcovtype+targetsample combined"
}

* ---- Test 66: grf_best_linear_projection HC3 SEs >= HC0 SEs (typical) ----
capture noisily {
    grf_causal_forest y w x1-x5, gen(tau_vc13a) ntrees(100) seed(42)
    grf_best_linear_projection x1 x2, vcovtype(HC0)
    local se_hc0_x1 = _se[x1]
    drop tau_vc13a _grf_yhat _grf_what

    grf_causal_forest y w x1-x5, gen(tau_vc13b) ntrees(100) seed(42)
    grf_best_linear_projection x1 x2, vcovtype(HC3)
    local se_hc3_x1 = _se[x1]
    drop tau_vc13b _grf_yhat _grf_what

    * HC3 SEs are typically larger than HC0 SEs
    assert `se_hc3_x1' >= `se_hc0_x1'
}
if _rc {
    display as error "FAIL: grf_best_linear_projection HC3 >= HC0 SEs"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_best_linear_projection HC3 >= HC0 SEs"
}

* ---- Test 67: grf_best_linear_projection e(vcov_type) stored correctly for each ----
capture noisily {
    local all_ok 1
    foreach hc in HC0 HC1 HC2 HC3 {
        grf_causal_forest y w x1-x5, gen(tau_vc14_`hc') ntrees(100) seed(42)
        grf_best_linear_projection x1, vcovtype(`hc')
        if "`e(vcov_type)'" != "`hc'" {
            local all_ok 0
        }
        drop tau_vc14_`hc' _grf_yhat _grf_what
    }
    assert `all_ok' == 1
}
if _rc {
    display as error "FAIL: grf_best_linear_projection e(vcov_type) stored for each HC"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_best_linear_projection e(vcov_type) stored for each HC"
}

* ---- Test 68: grf_best_linear_projection with if and vcovtype ----
capture noisily {
    grf_causal_forest y w x1-x5, gen(tau_vc15) ntrees(100) seed(42)
    grf_best_linear_projection x1 if x1 > 0, vcovtype(HC1)
    assert e(N) > 0
    assert e(N) < 500
    assert "`e(vcov_type)'" == "HC1"
    drop tau_vc15 _grf_yhat _grf_what
}
if _rc {
    display as error "FAIL: grf_best_linear_projection with if and vcovtype"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_best_linear_projection with if and vcovtype"
}

* ============================================================
* grf_rate compliance.score tests
* ============================================================

* ---- Test 69: grf_rate without compliancescore (backward compat) ----
capture noisily {
    grf_causal_forest y w x1-x5, gen(tau_cs1) ntrees(100) seed(42)
    grf_rate tau_cs1, bootstrap(50) seed(42)
    assert !missing(r(estimate))
    assert !missing(r(std_err))
    assert r(std_err) > 0
    drop tau_cs1 _grf_yhat _grf_what
}
if _rc {
    display as error "FAIL: grf_rate without compliancescore"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_rate without compliancescore"
}

* ---- Test 70: grf_rate with compliancescore ----
capture noisily {
    * Generate a compliance-like score variable
    gen comp_score = abs(rnormal()) / 2
    replace comp_score = 0.01 if comp_score < 0.01

    grf_causal_forest y w x1-x5, gen(tau_cs2) ntrees(100) seed(42)
    grf_rate tau_cs2, compliancescore(comp_score) bootstrap(50) seed(42)
    assert !missing(r(estimate))
    assert !missing(r(std_err))
    assert r(std_err) > 0
    assert r(n) > 0
    drop tau_cs2 _grf_yhat _grf_what
}
if _rc {
    display as error "FAIL: grf_rate with compliancescore"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_rate with compliancescore"
}

* ---- Test 71: grf_rate compliancescore changes estimate ----
capture noisily {
    grf_causal_forest y w x1-x5, gen(tau_cs3) ntrees(100) seed(42)
    grf_rate tau_cs3, bootstrap(50) seed(42)
    local est_no_comp = r(estimate)
    drop tau_cs3 _grf_yhat _grf_what

    grf_causal_forest y w x1-x5, gen(tau_cs3b) ntrees(100) seed(42)
    grf_rate tau_cs3b, compliancescore(comp_score) bootstrap(50) seed(42)
    local est_with_comp = r(estimate)
    drop tau_cs3b _grf_yhat _grf_what

    * Compliance weighting should produce a different estimate
    assert reldif(`est_no_comp', `est_with_comp') > 1e-6
}
if _rc {
    display as error "FAIL: grf_rate compliancescore changes estimate"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_rate compliancescore changes estimate"
}

* ---- Test 72: grf_rate compliancescore returns compliance_score_var ----
capture noisily {
    grf_causal_forest y w x1-x5, gen(tau_cs4) ntrees(100) seed(42)
    grf_rate tau_cs4, compliancescore(comp_score) bootstrap(50) seed(42)
    assert "`r(compliance_score_var)'" == "comp_score"
    drop tau_cs4 _grf_yhat _grf_what
}
if _rc {
    display as error "FAIL: grf_rate compliancescore returns compliance_score_var"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_rate compliancescore returns compliance_score_var"
}

* ---- Test 73: grf_rate without compliancescore has no compliance_score_var ----
capture noisily {
    grf_causal_forest y w x1-x5, gen(tau_cs5) ntrees(100) seed(42)
    grf_rate tau_cs5, bootstrap(50) seed(42)
    assert "`r(compliance_score_var)'" == ""
    drop tau_cs5 _grf_yhat _grf_what
}
if _rc {
    display as error "FAIL: grf_rate without compliancescore has empty macro"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_rate without compliancescore has empty macro"
}

* ---- Test 74: grf_rate compliancescore with QINI target ----
capture noisily {
    grf_causal_forest y w x1-x5, gen(tau_cs6) ntrees(100) seed(42)
    grf_rate tau_cs6, compliancescore(comp_score) target(QINI) bootstrap(50) seed(42)
    assert !missing(r(estimate))
    assert !missing(r(std_err))
    assert "`r(target)'" == "QINI"
    drop tau_cs6 _grf_yhat _grf_what
}
if _rc {
    display as error "FAIL: grf_rate compliancescore with QINI"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_rate compliancescore with QINI"
}

* ---- Test 75: grf_rate compliancescore with if restriction ----
capture noisily {
    grf_causal_forest y w x1-x5, gen(tau_cs7) ntrees(100) seed(42)
    grf_rate tau_cs7 if x1 > 0, compliancescore(comp_score) bootstrap(50) seed(42)
    assert !missing(r(estimate))
    assert r(n) < 500
    drop tau_cs7 _grf_yhat _grf_what
}
if _rc {
    display as error "FAIL: grf_rate compliancescore with if"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_rate compliancescore with if"
}

* ---- Test 76: grf_rate compliancescore SE is positive ----
capture noisily {
    grf_causal_forest y w x1-x5, gen(tau_cs8) ntrees(100) seed(42)
    grf_rate tau_cs8, compliancescore(comp_score) bootstrap(50) seed(42)
    assert r(std_err) > 0
    assert r(p_value) >= 0 & r(p_value) <= 1
    drop tau_cs8 _grf_yhat _grf_what
}
if _rc {
    display as error "FAIL: grf_rate compliancescore SE positive"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_rate compliancescore SE positive"
}

* ---- Test 77: grf_rate compliancescore with uniform scores ----
capture noisily {
    * Uniform compliance scores should still work
    gen comp_uniform = 1
    grf_causal_forest y w x1-x5, gen(tau_cs9) ntrees(100) seed(42)
    grf_rate tau_cs9, compliancescore(comp_uniform) bootstrap(50) seed(42)
    local est_uniform = r(estimate)
    drop tau_cs9 _grf_yhat _grf_what

    * With uniform compliance = 1, estimate should match no-compliance result
    grf_causal_forest y w x1-x5, gen(tau_cs9b) ntrees(100) seed(42)
    grf_rate tau_cs9b, bootstrap(50) seed(42)
    local est_none = r(estimate)
    drop tau_cs9b _grf_yhat _grf_what

    * Should be very close (both use same DR scores, just scaled by 1)
    assert reldif(`est_uniform', `est_none') < 0.05
    drop comp_uniform
}
if _rc {
    display as error "FAIL: grf_rate compliancescore with uniform scores"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_rate compliancescore with uniform scores"
}

* ---- Test 78: grf_rate compliancescore n_bootstrap returned ----
capture noisily {
    grf_causal_forest y w x1-x5, gen(tau_cs10) ntrees(100) seed(42)
    grf_rate tau_cs10, compliancescore(comp_score) bootstrap(50) seed(42)
    assert r(n_bootstrap) > 0
    assert r(n_bootstrap) <= 50
    drop tau_cs10 _grf_yhat _grf_what
}
if _rc {
    display as error "FAIL: grf_rate compliancescore n_bootstrap returned"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_rate compliancescore n_bootstrap returned"
}

* Clean up compliance score
capture drop comp_score

* ============================================================
* grf_best_linear_projection targetsample(treated/control) tests
* ============================================================

* ---- Test 79: BLP targetsample(treated) runs ----
capture noisily {
    grf_causal_forest y w x1-x5, gen(tau_blpt) ntrees(100) seed(42)
    grf_best_linear_projection x1 x2, targetsample(treated)
    assert "`e(target_sample)'" == "treated"
    assert !missing(e(N))
    drop tau_blpt _grf_yhat _grf_what
}
if _rc {
    display as error "FAIL: BLP targetsample(treated) runs"
    local errors = `errors' + 1
}
else {
    display as result "PASS: BLP targetsample(treated) runs"
}

* ---- Test 80: BLP targetsample(control) runs ----
capture noisily {
    grf_causal_forest y w x1-x5, gen(tau_blpc) ntrees(100) seed(42)
    grf_best_linear_projection x1 x2, targetsample(control)
    assert "`e(target_sample)'" == "control"
    assert !missing(e(N))
    drop tau_blpc _grf_yhat _grf_what
}
if _rc {
    display as error "FAIL: BLP targetsample(control) runs"
    local errors = `errors' + 1
}
else {
    display as result "PASS: BLP targetsample(control) runs"
}

* ---- Test 81: BLP treated vs all coefficients differ ----
capture noisily {
    grf_causal_forest y w x1-x5, gen(tau_blp81) ntrees(100) seed(42)

    * All
    grf_best_linear_projection x1 x2, targetsample(all)
    local b1_all = _b[x1]

    * Re-fit to get back causal forest e() results
    grf_causal_forest y w x1-x5, gen(tau_blp81b) ntrees(100) seed(42) replace

    * Treated
    grf_best_linear_projection x1 x2, targetsample(treated)
    local b1_treated = _b[x1]

    * The coefficients should exist and differ
    assert !missing(`b1_all')
    assert !missing(`b1_treated')
    * With overlap weighting vs uniform, coefficients should generally differ
    * (they could be equal by chance, but with enough data this is unlikely)
    display "  b1_all = `b1_all', b1_treated = `b1_treated'"

    drop tau_blp81 tau_blp81b _grf_yhat _grf_what
}
if _rc {
    display as error "FAIL: BLP treated vs all coefficients differ"
    local errors = `errors' + 1
}
else {
    display as result "PASS: BLP treated vs all coefficients differ"
}

* ---- Test 82: BLP control vs all coefficients differ ----
capture noisily {
    grf_causal_forest y w x1-x5, gen(tau_blp82) ntrees(100) seed(42)

    * All
    grf_best_linear_projection x1 x2, targetsample(all)
    local b1_all = _b[x1]

    * Re-fit to get back causal forest e() results
    grf_causal_forest y w x1-x5, gen(tau_blp82b) ntrees(100) seed(42) replace

    * Control
    grf_best_linear_projection x1 x2, targetsample(control)
    local b1_control = _b[x1]

    * The coefficients should exist and differ
    assert !missing(`b1_all')
    assert !missing(`b1_control')
    display "  b1_all = `b1_all', b1_control = `b1_control'"

    drop tau_blp82 tau_blp82b _grf_yhat _grf_what
}
if _rc {
    display as error "FAIL: BLP control vs all coefficients differ"
    local errors = `errors' + 1
}
else {
    display as result "PASS: BLP control vs all coefficients differ"
}

* ============================================================
* grf_get_scores: multi-arm causal forest tests
* ============================================================

* ---- Test 83: grf_get_scores after multi_arm_causal_forest ----
capture noisily {
    * Create multi-arm treatment data: 2 binary treatment indicators
    capture drop w1 w2
    gen w1 = (x1 + rnormal() > 0)
    gen w2 = (x2 + rnormal() > 0)

    * Fit multi-arm causal forest: y w1 w2 x1..x5, ntreat(2)
    grf_multi_arm_causal_forest y w1 w2 x1-x5, gen(mac_tau) ntreat(2) ///
        ntrees(100) seed(42)

    * Get DR scores — should create 2 score columns
    grf_get_scores, gen(mac_dr) replace

    * Verify 2 score columns exist
    confirm numeric variable mac_dr_t1
    confirm numeric variable mac_dr_t2

    * Verify scores are non-missing for observations in sample
    quietly count if !missing(mac_dr_t1)
    assert r(N) > 0
    quietly count if !missing(mac_dr_t2)
    assert r(N) > 0

    * Verify return results
    assert r(n_treat) == 2
    assert r(N) > 0
    assert "`r(forest_type)'" == "multi_causal"

    * Cleanup
    capture drop mac_tau_t1 mac_tau_t2
    capture drop mac_tau_t1_var mac_tau_t2_var
    capture drop mac_dr_t1 mac_dr_t2
    capture drop w1 w2
    capture drop _grf_mac_yhat _grf_mac_what1 _grf_mac_what2
}
if _rc {
    display as error "FAIL: grf_get_scores after multi_arm_causal_forest"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_get_scores after multi_arm_causal_forest"
}

* ---- Test 84: grf_get_scores multi-arm replace option ----
capture noisily {
    capture drop w1 w2
    gen w1 = (x1 + rnormal() > 0.2)
    gen w2 = (x2 + rnormal() > 0.2)

    grf_multi_arm_causal_forest y w1 w2 x1-x5, gen(mac2_tau) ntreat(2) ///
        ntrees(100) seed(42)

    * First call
    grf_get_scores, gen(mac2_dr)

    * Second call with replace — should not error
    grf_get_scores, gen(mac2_dr) replace

    confirm numeric variable mac2_dr_t1
    confirm numeric variable mac2_dr_t2

    capture drop mac2_tau_t1 mac2_tau_t2
    capture drop mac2_tau_t1_var mac2_tau_t2_var
    capture drop mac2_dr_t1 mac2_dr_t2
    capture drop w1 w2
    capture drop _grf_mac_yhat _grf_mac_what1 _grf_mac_what2
}
if _rc {
    display as error "FAIL: grf_get_scores multi-arm replace option"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_get_scores multi-arm replace option"
}

* ============================================================
* grf_get_scores: causal survival forest tests
* ============================================================

* ---- Test 85: grf_get_scores after causal_survival_forest ----
capture noisily {
    * Create survival data
    capture drop surv_time surv_status surv_treat
    gen surv_treat = (x1 + rnormal() > 0)
    gen surv_time = abs(2 + x1 + surv_treat * 0.5 + rnormal())
    replace surv_time = 0.01 if surv_time <= 0
    gen surv_status = (runiform() > 0.3)

    * Fit causal survival forest: time status treat x1..x5
    grf_causal_survival_forest surv_time surv_status surv_treat x1-x5, ///
        gen(cs_tau) ntrees(100) seed(42)

    * Get DR scores
    grf_get_scores, gen(cs_dr)

    * Verify score variable exists
    confirm numeric variable cs_dr

    * Verify scores are non-missing
    quietly count if !missing(cs_dr)
    assert r(N) > 0

    * Verify return results
    assert r(N) > 0
    assert "`r(forest_type)'" == "causal_survival"

    * Cleanup
    capture drop cs_tau cs_dr cs_tau_var
    capture drop surv_time surv_status surv_treat
    capture drop _grf_cs_what
}
if _rc {
    display as error "FAIL: grf_get_scores after causal_survival_forest"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_get_scores after causal_survival_forest"
}

* ---- Test 86: grf_get_scores causal_survival replace option ----
capture noisily {
    capture drop surv_time surv_status surv_treat
    gen surv_treat = (x1 + rnormal() > 0)
    gen surv_time = abs(2 + x1 + surv_treat * 0.5 + rnormal())
    replace surv_time = 0.01 if surv_time <= 0
    gen surv_status = (runiform() > 0.3)

    grf_causal_survival_forest surv_time surv_status surv_treat x1-x5, ///
        gen(cs2_tau) ntrees(100) seed(42)

    * First call
    grf_get_scores, gen(cs2_dr)

    * Second call with replace
    grf_get_scores, gen(cs2_dr) replace

    confirm numeric variable cs2_dr
    quietly count if !missing(cs2_dr)
    assert r(N) > 0

    capture drop cs2_tau cs2_dr cs2_tau_var
    capture drop surv_time surv_status surv_treat
    capture drop _grf_cs_what
}
if _rc {
    display as error "FAIL: grf_get_scores causal_survival replace option"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_get_scores causal_survival replace option"
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
