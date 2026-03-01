* test_tune_extended.do -- Test grf_tune for all 13 forest types
* Tests that tuning runs and returns valid hyperparameters for each forest type

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
gen z = (rnormal() > 0)
gen y_binary = (y > 0)
gen y2 = x2 + rnormal()
* Multi-arm treatment (0, 1, 2)
gen treat_multi = floor(3 * runiform())
* Survival data
gen T_surv = -ln(runiform()) / exp(0.5 * x1)
gen C_surv = -ln(runiform()) / 0.5
gen time_surv = min(T_surv, C_surv)
gen status_surv = (T_surv <= C_surv)

* ============================================================
* Test 1: grf_tune regression forest
* ============================================================

* ---- Test 1: grf_tune foresttype("regression") ----
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

* ============================================================
* Test 2: grf_tune causal forest
* ============================================================

* ---- Test 2: grf_tune foresttype("causal") ----
capture noisily {
    grf_tune y w x1-x5, foresttype("causal") numreps(10) tunetrees(50) seed(42)
    assert !missing(r(best_mtry))
    assert r(best_mtry) >= 1
    assert !missing(r(best_min_node_size))
    assert r(best_min_node_size) >= 1
    assert !missing(r(best_sample_fraction))
    assert r(best_sample_fraction) > 0 & r(best_sample_fraction) < 1
    assert !missing(r(best_mse))
    assert r(best_mse) > 0
    assert r(N) == 500
    assert "`r(forest_type)'" == "causal"
}
if _rc {
    display as error "FAIL: grf_tune causal"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_tune causal"
}

* ============================================================
* Test 3: grf_tune quantile forest
* ============================================================

* ---- Test 3: grf_tune foresttype("quantile") ----
capture noisily {
    grf_tune y x1-x5, foresttype("quantile") numreps(10) tunetrees(50) seed(42)
    assert !missing(r(best_mtry))
    assert r(best_mtry) >= 1
    assert !missing(r(best_min_node_size))
    assert r(best_min_node_size) >= 1
    assert !missing(r(best_sample_fraction))
    assert r(best_sample_fraction) > 0 & r(best_sample_fraction) < 1
    assert !missing(r(best_mse))
    assert r(best_mse) > 0
    assert r(N) == 500
    assert "`r(forest_type)'" == "quantile"
}
if _rc {
    display as error "FAIL: grf_tune quantile"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_tune quantile"
}

* ============================================================
* Test 4: grf_tune probability forest
* ============================================================

* ---- Test 4: grf_tune foresttype("probability") ----
capture noisily {
    grf_tune y_binary x1-x5, foresttype("probability") numreps(10) tunetrees(50) seed(42)
    assert !missing(r(best_mtry))
    assert r(best_mtry) >= 1
    assert !missing(r(best_min_node_size))
    assert r(best_min_node_size) >= 1
    assert !missing(r(best_sample_fraction))
    assert r(best_sample_fraction) > 0 & r(best_sample_fraction) < 1
    assert !missing(r(best_mse))
    assert r(best_mse) > 0
    assert r(N) == 500
    assert "`r(forest_type)'" == "probability"
}
if _rc {
    display as error "FAIL: grf_tune probability"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_tune probability"
}

* ============================================================
* Test 5: grf_tune instrumental forest
* ============================================================

* ---- Test 5: grf_tune foresttype("instrumental") ----
capture noisily {
    grf_tune y w z x1-x5, foresttype("instrumental") numreps(10) tunetrees(50) seed(42)
    assert !missing(r(best_mtry))
    assert r(best_mtry) >= 1
    assert !missing(r(best_min_node_size))
    assert r(best_min_node_size) >= 1
    assert !missing(r(best_sample_fraction))
    assert r(best_sample_fraction) > 0 & r(best_sample_fraction) < 1
    assert !missing(r(best_mse))
    assert r(best_mse) > 0
    assert r(N) == 500
    assert "`r(forest_type)'" == "instrumental"
}
if _rc {
    display as error "FAIL: grf_tune instrumental"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_tune instrumental"
}

* ============================================================
* Test 6: grf_tune survival forest
* ============================================================

* ---- Test 6: grf_tune foresttype("survival") ----
capture noisily {
    grf_tune time_surv status_surv x1-x5, foresttype("survival") numreps(10) tunetrees(50) seed(42)
    assert !missing(r(best_mtry))
    assert r(best_mtry) >= 1
    assert !missing(r(best_min_node_size))
    assert r(best_min_node_size) >= 1
    assert !missing(r(best_sample_fraction))
    assert r(best_sample_fraction) > 0 & r(best_sample_fraction) < 1
    assert !missing(r(best_mse))
    assert r(best_mse) > 0
    assert r(N) == 500
    assert "`r(forest_type)'" == "survival"
}
if _rc {
    display as error "FAIL: grf_tune survival"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_tune survival"
}

* ============================================================
* Test 7: grf_tune ll_regression forest
* ============================================================

* ---- Test 7: grf_tune foresttype("ll_regression") ----
capture noisily {
    grf_tune y x1-x5, foresttype("ll_regression") numreps(10) tunetrees(50) seed(42)
    assert !missing(r(best_mtry))
    assert r(best_mtry) >= 1
    assert !missing(r(best_min_node_size))
    assert r(best_min_node_size) >= 1
    assert !missing(r(best_sample_fraction))
    assert r(best_sample_fraction) > 0 & r(best_sample_fraction) < 1
    assert !missing(r(best_mse))
    assert r(best_mse) > 0
    assert r(N) == 500
    assert "`r(forest_type)'" == "ll_regression"
}
if _rc {
    display as error "FAIL: grf_tune ll_regression"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_tune ll_regression"
}

* ============================================================
* Test 8: grf_tune lm_forest (requires xvars)
* ============================================================

* ---- Test 8: grf_tune foresttype("lm_forest") ----
capture noisily {
    grf_tune y w, foresttype("lm_forest") xvars(x1 x2 x3 x4 x5) numreps(10) tunetrees(50) seed(42)
    assert !missing(r(best_mtry))
    assert r(best_mtry) >= 1
    assert !missing(r(best_min_node_size))
    assert r(best_min_node_size) >= 1
    assert !missing(r(best_sample_fraction))
    assert r(best_sample_fraction) > 0 & r(best_sample_fraction) < 1
    assert !missing(r(best_mse))
    assert r(best_mse) > 0
    assert r(N) == 500
    assert "`r(forest_type)'" == "lm_forest"
}
if _rc {
    display as error "FAIL: grf_tune lm_forest"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_tune lm_forest"
}

* ============================================================
* Test 9: grf_tune multi_arm_causal forest (requires ntreat)
* ============================================================

* ---- Test 9: grf_tune foresttype("multi_arm_causal") ----
capture noisily {
    grf_tune y treat_multi x1-x5, foresttype("multi_arm_causal") ntreat(1) numreps(10) tunetrees(50) seed(42)
    assert !missing(r(best_mtry))
    assert r(best_mtry) >= 1
    assert !missing(r(best_min_node_size))
    assert r(best_min_node_size) >= 1
    assert !missing(r(best_sample_fraction))
    assert r(best_sample_fraction) > 0 & r(best_sample_fraction) < 1
    assert !missing(r(best_mse))
    assert r(best_mse) > 0
    assert r(N) == 500
    assert "`r(forest_type)'" == "multi_arm_causal"
}
if _rc {
    display as error "FAIL: grf_tune multi_arm_causal"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_tune multi_arm_causal"
}

* ============================================================
* Test 10: grf_tune multi_regression forest (requires ndep)
* ============================================================

* ---- Test 10: grf_tune foresttype("multi_regression") ----
capture noisily {
    grf_tune y y2 x1-x5, foresttype("multi_regression") ndep(2) numreps(10) tunetrees(50) seed(42)
    assert !missing(r(best_mtry))
    assert r(best_mtry) >= 1
    assert !missing(r(best_min_node_size))
    assert r(best_min_node_size) >= 1
    assert !missing(r(best_sample_fraction))
    assert r(best_sample_fraction) > 0 & r(best_sample_fraction) < 1
    assert !missing(r(best_mse))
    assert r(best_mse) > 0
    assert r(N) == 500
    assert "`r(forest_type)'" == "multi_regression"
}
if _rc {
    display as error "FAIL: grf_tune multi_regression"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_tune multi_regression"
}

* ============================================================
* Test 11: grf_tune causal_survival forest
* ============================================================

* ---- Test 11: grf_tune foresttype("causal_survival") ----
capture noisily {
    grf_tune time_surv status_surv w x1-x5, foresttype("causal_survival") numreps(10) tunetrees(50) seed(42)
    assert !missing(r(best_mtry))
    assert r(best_mtry) >= 1
    assert !missing(r(best_min_node_size))
    assert r(best_min_node_size) >= 1
    assert !missing(r(best_sample_fraction))
    assert r(best_sample_fraction) > 0 & r(best_sample_fraction) < 1
    assert !missing(r(best_mse))
    assert r(best_mse) > 0
    assert r(N) == 500
    assert "`r(forest_type)'" == "causal_survival"
}
if _rc {
    display as error "FAIL: grf_tune causal_survival"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_tune causal_survival"
}

* ============================================================
* Test 12: grf_tune boosted_regression forest
* ============================================================

* ---- Test 12: grf_tune foresttype("boosted_regression") ----
capture noisily {
    grf_tune y x1-x5, foresttype("boosted_regression") numreps(10) tunetrees(50) seed(42)
    assert !missing(r(best_mtry))
    assert r(best_mtry) >= 1
    assert !missing(r(best_min_node_size))
    assert r(best_min_node_size) >= 1
    assert !missing(r(best_sample_fraction))
    assert r(best_sample_fraction) > 0 & r(best_sample_fraction) < 1
    assert !missing(r(best_mse))
    assert r(best_mse) > 0
    assert r(N) == 500
    assert "`r(forest_type)'" == "boosted_regression"
}
if _rc {
    display as error "FAIL: grf_tune boosted_regression"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_tune boosted_regression"
}

* ============================================================
* Test 13: grf_tune probability with nclasses(3)
* ============================================================

* ---- Test 13: grf_tune probability with nclasses(3) ----
capture noisily {
    * Create 3-class outcome
    capture drop y3class
    gen y3class = floor(3 * runiform())
    grf_tune y3class x1-x5, foresttype("probability") nclasses(3) ///
        numreps(10) tunetrees(50) seed(42)
    assert !missing(r(best_mtry))
    assert r(best_mtry) >= 1
    assert !missing(r(best_mse))
    assert r(best_mse) > 0
    assert "`r(forest_type)'" == "probability"
    drop y3class
}
if _rc {
    display as error "FAIL: grf_tune probability nclasses(3)"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_tune probability nclasses(3)"
}

* ============================================================
* Test 14: grf_tune survival with numfailures(10)
* ============================================================

* ---- Test 14: grf_tune survival numfailures(10) ----
capture noisily {
    grf_tune time_surv status_surv x1-x5, foresttype("survival") ///
        numfailures(10) numreps(10) tunetrees(50) seed(42)
    assert !missing(r(best_mtry))
    assert r(best_mtry) >= 1
    assert !missing(r(best_mse))
    assert r(best_mse) > 0
    assert "`r(forest_type)'" == "survival"
}
if _rc {
    display as error "FAIL: grf_tune survival numfailures(10)"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_tune survival numfailures(10)"
}

* ============================================================
* Test 15: grf_tune survival with predtype(0) (Nelson-Aalen)
* ============================================================

* ---- Test 15: grf_tune survival predtype(0) ----
capture noisily {
    grf_tune time_surv status_surv x1-x5, foresttype("survival") ///
        predtype(0) numreps(10) tunetrees(50) seed(42)
    assert !missing(r(best_mtry))
    assert r(best_mtry) >= 1
    assert !missing(r(best_mse))
    assert r(best_mse) > 0
    assert "`r(forest_type)'" == "survival"
}
if _rc {
    display as error "FAIL: grf_tune survival predtype(0)"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_tune survival predtype(0)"
}

* ============================================================
* Test 16: grf_tune instrumental with reducedformweight(0.5)
* ============================================================

* ---- Test 16: grf_tune instrumental reducedformweight(0.5) ----
capture noisily {
    grf_tune y w z x1-x5, foresttype("instrumental") ///
        reducedformweight(0.5) numreps(10) tunetrees(50) seed(42)
    assert !missing(r(best_mtry))
    assert r(best_mtry) >= 1
    assert !missing(r(best_mse))
    assert r(best_mse) > 0
    assert "`r(forest_type)'" == "instrumental"
}
if _rc {
    display as error "FAIL: grf_tune instrumental reducedformweight(0.5)"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_tune instrumental reducedformweight(0.5)"
}

* ============================================================
* Test 17: grf_tune causal_survival with horizon(1.5)
* ============================================================

* ---- Test 17: grf_tune causal_survival horizon(1.5) ----
capture noisily {
    grf_tune time_surv status_surv w x1-x5, foresttype("causal_survival") ///
        horizon(1.5) numreps(10) tunetrees(50) seed(42)
    assert !missing(r(best_mtry))
    assert r(best_mtry) >= 1
    assert !missing(r(best_mse))
    assert r(best_mse) > 0
    assert "`r(forest_type)'" == "causal_survival"
}
if _rc {
    display as error "FAIL: grf_tune causal_survival horizon(1.5)"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_tune causal_survival horizon(1.5)"
}

* ============================================================
* Test 18: grf_tune causal_survival with target(2)
* ============================================================

* ---- Test 18: grf_tune causal_survival target(2) ----
capture noisily {
    grf_tune time_surv status_surv w x1-x5, foresttype("causal_survival") ///
        target(2) numreps(10) tunetrees(50) seed(42)
    assert !missing(r(best_mtry))
    assert r(best_mtry) >= 1
    assert !missing(r(best_mse))
    assert r(best_mse) > 0
    assert "`r(forest_type)'" == "causal_survival"
}
if _rc {
    display as error "FAIL: grf_tune causal_survival target(2)"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_tune causal_survival target(2)"
}

* ============================================================
* Test 19: grf_tune returns consistent results with same seed
* ============================================================

* ---- Test 19: grf_tune reproducibility with same seed ----
capture noisily {
    grf_tune y x1-x5, foresttype("regression") numreps(10) tunetrees(50) seed(99)
    local mtry_1 = r(best_mtry)
    local mns_1 = r(best_min_node_size)
    local sf_1 = r(best_sample_fraction)
    local mse_1 = r(best_mse)

    grf_tune y x1-x5, foresttype("regression") numreps(10) tunetrees(50) seed(99)
    local mtry_2 = r(best_mtry)
    local mns_2 = r(best_min_node_size)
    local sf_2 = r(best_sample_fraction)
    local mse_2 = r(best_mse)

    assert `mtry_1' == `mtry_2'
    assert `mns_1' == `mns_2'
    assert reldif(`sf_1', `sf_2') < 1e-10
    assert reldif(`mse_1', `mse_2') < 1e-10
}
if _rc {
    display as error "FAIL: grf_tune reproducibility"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_tune reproducibility"
}

* ============================================================
* Summary
* ============================================================
display ""
if `errors' > 0 {
    display as error "FAILED: `errors' errors in tune extended tests"
    exit 1
}
else {
    display as result "ALL TUNE EXTENDED TESTS PASSED (19 tests: 12 forest types + 6 option tests + reproducibility)"
}
