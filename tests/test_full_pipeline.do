* test_full_pipeline.do -- End-to-end integration test
* Exercises the complete grf Stata workflow: data generation, forest fitting,
* post-estimation (ATE, BLP, RATE, APE, scores, calibration, prediction,
* variable importance), and internal consistency checks.
*
* Run: do tests/test_full_pipeline.do

clear all
set more off

local errors = 0

* ============================================================
* Step 1: Generate causal data
* ============================================================
capture noisily {
    grf_generate_causal_data, n(800) p(5) dgp(simple) seed(42)
    assert _N == 800
    confirm variable Y
    confirm variable W
    confirm variable X1
    confirm variable X5
    confirm variable tau_true
}
if _rc {
    display as error "FAIL: Step 1 - generate causal data"
    local errors = `errors' + 1
}
else {
    display as result "PASS: Step 1 - generate causal data"
}

* ============================================================
* Step 2: Fit causal forest on training sample (first 600 obs)
* ============================================================
capture noisily {
    grf_causal_forest Y W X1-X5 in 1/600, gen(cate) ntrees(500) seed(42)
    assert e(N) == 600
    assert "`e(forest_type)'" == "causal"
    quietly count if !missing(cate) & _n <= 600
    assert r(N) == 600
}
if _rc {
    display as error "FAIL: Step 2 - fit causal forest"
    local errors = `errors' + 1
}
else {
    display as result "PASS: Step 2 - fit causal forest"
}

* ============================================================
* Step 3: ATE (default = all)
* ============================================================
capture noisily {
    grf_ate
    assert !missing(r(ate))
    assert !missing(r(se))
    assert r(se) > 0
    local ate_all = r(ate)
    local ate_se  = r(se)
}
if _rc {
    display as error "FAIL: Step 3 - ATE (all)"
    local errors = `errors' + 1
}
else {
    display as result "PASS: Step 3 - ATE (all), estimate = `ate_all'"
}

* ============================================================
* Step 4: ATE with overlap weighting
* ============================================================
capture noisily {
    grf_ate, targetsample(overlap)
    assert !missing(r(ate))
    assert !missing(r(se))
    local ate_overlap = r(ate)
}
if _rc {
    display as error "FAIL: Step 4 - ATE (overlap)"
    local errors = `errors' + 1
}
else {
    display as result "PASS: Step 4 - ATE (overlap), estimate = `ate_overlap'"
}

* ============================================================
* Step 5: Get DR scores
* ============================================================
capture noisily {
    grf_get_scores, gen(scores)
    quietly count if !missing(scores)
    assert r(N) > 0
}
if _rc {
    display as error "FAIL: Step 5 - get DR scores"
    local errors = `errors' + 1
}
else {
    display as result "PASS: Step 5 - get DR scores"
}

* ============================================================
* Step 6: Best linear projection
* NOTE: BLP is an eclass command -- it replaces e() results.
*       After this step, e(cmd) is no longer grf_causal_forest.
* ============================================================
capture noisily {
    grf_best_linear_projection X1 X2, vcovtype(HC3)
    assert e(N) > 0
}
if _rc {
    display as error "FAIL: Step 6 - best linear projection"
    local errors = `errors' + 1
}
else {
    display as result "PASS: Step 6 - best linear projection"
}

* ============================================================
* Step 7: RATE (re-fit causal forest since BLP replaced e())
* ============================================================
capture noisily {
    grf_causal_forest Y W X1-X5 in 1/600, gen(cate2) ntrees(500) seed(42) replace
    grf_rate cate2, target(AUTOC) bootstrap(50) seed(42)
    assert !missing(r(estimate))
    assert !missing(r(std_err))
    assert r(std_err) > 0
}
if _rc {
    display as error "FAIL: Step 7 - RATE (AUTOC)"
    local errors = `errors' + 1
}
else {
    display as result "PASS: Step 7 - RATE (AUTOC)"
}

* ============================================================
* Step 8: Average partial effect (skip if command doesn't exist)
* ============================================================
capture which grf_average_partial_effect
if _rc {
    display as text "SKIP: Step 8 - average partial effect (command not found)"
}
else {
    capture noisily {
        * Re-fit causal forest (RATE may have altered state)
        grf_causal_forest Y W X1-X5 in 1/600, gen(cate3) ntrees(500) seed(42) replace
        grf_average_partial_effect
        assert !missing(r(estimate))
        assert !missing(r(std_err))
    }
    if _rc {
        display as error "FAIL: Step 8 - average partial effect"
        local errors = `errors' + 1
    }
    else {
        display as result "PASS: Step 8 - average partial effect"
    }
}

* ============================================================
* Step 9: Predict on test observations (obs 601-800)
* ============================================================
capture noisily {
    * Re-fit causal forest for clean e() state
    grf_causal_forest Y W X1-X5 in 1/600, gen(cate4) ntrees(500) seed(42) replace
    grf_predict, gen(cate_new)

    * Check that test obs (601-800) have predictions
    quietly count if !missing(cate_new) & _n > 600
    assert r(N) > 0
}
if _rc {
    display as error "FAIL: Step 9 - predict on test data"
    local errors = `errors' + 1
}
else {
    display as result "PASS: Step 9 - predict on test data"
}

* ============================================================
* Step 10: Variable importance
* ============================================================
capture noisily {
    grf_variable_importance Y X1-X5 in 1/600, ntrees(500) seed(42)
    assert r(N) == 600
    matrix vi = r(importance)
    assert colsof(vi) == 5
    * All importance values should be non-negative
    forvalues j = 1/5 {
        assert vi[1,`j'] >= 0
    }
    * Sum should be close to 1
    local vi_sum = 0
    forvalues j = 1/5 {
        local vi_sum = `vi_sum' + vi[1,`j']
    }
    assert reldif(`vi_sum', 1.0) < 0.01
}
if _rc {
    display as error "FAIL: Step 10 - variable importance"
    local errors = `errors' + 1
}
else {
    display as result "PASS: Step 10 - variable importance"
}

* ============================================================
* Step 11: Test calibration (re-fit causal forest first)
* ============================================================
capture noisily {
    grf_causal_forest Y W X1-X5 in 1/600, gen(cate5) ntrees(500) seed(42) replace
    grf_test_calibration
    assert !missing(r(b_mean))
    assert !missing(r(se_mean))
    assert !missing(r(b_diff))
    assert !missing(r(se_diff))
    assert r(se_mean) > 0
    assert r(se_diff) > 0
}
if _rc {
    display as error "FAIL: Step 11 - test calibration"
    local errors = `errors' + 1
}
else {
    display as result "PASS: Step 11 - test calibration"
}

* ============================================================
* Step 12: Consistency check -- mean(scores) should approximate ATE
* The DR scores from Step 5 should have mean ≈ ATE from Step 3
* ============================================================
capture noisily {
    quietly summarize scores in 1/600
    local score_mean = r(mean)
    local diff = abs(`score_mean' - `ate_all')
    local tolerance = 3 * `ate_se'
    display as text "  Score mean: `score_mean', ATE: `ate_all', tolerance: `tolerance'"
    assert `diff' < `tolerance'
}
if _rc {
    display as error "FAIL: Step 12 - consistency check (mean scores ≈ ATE)"
    local errors = `errors' + 1
}
else {
    display as result "PASS: Step 12 - consistency check (mean scores ≈ ATE)"
}

* ============================================================
* Summary
* ============================================================
display ""
if `errors' > 0 {
    display as error "FAILED: `errors' pipeline step(s) did not pass"
    exit 1
}
else {
    display as result "ALL PIPELINE TESTS PASSED (12 steps)"
}
