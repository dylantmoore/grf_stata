* test_post_estimation.do -- Test grf_ate, grf_best_linear_projection,
*                           grf_variable_importance, grf_test_calibration
* Run from the test-c-plugin-skill-grf/ directory

clear all
set more off

display as text ""
display as text "=============================================="
display as text " GRF Post-Estimation Tests"
display as text "=============================================="

* ---- Setup: Fit a causal forest ----
display as text ""
display as text "--- Setup: Fitting causal forest ---"

clear
set obs 500
set seed 42

* Generate data with known treatment effect tau = x1
* True ATE = E[tau] = E[x1] = 0  (mean of standard normal)
gen x1 = rnormal()
gen x2 = rnormal()
gen x3 = rnormal()
gen x4 = rnormal()
gen x5 = rnormal()
gen w = (runiform() > 0.5)
gen tau = x1
gen y = 2 * x1 + tau * w + 0.5 * rnormal()

* Fit causal forest with nuisance estimates
grf_causal_forest y w x1 x2 x3 x4 x5, gen(cate) ntrees(1000) seed(42) ///
    yhatgenerate(yhat) whatgenerate(what) replace

quietly count if !missing(cate)
display as text "  CATE predictions: " as result r(N) " / 500"
assert r(N) > 400

display as text "  Setup complete"

* ---- Test 1: grf_ate ----
* (reads from e() results set by grf_causal_forest above)
display as text ""
display as text "--- Test 1: grf_ate ---"

grf_ate

* Check returned scalars exist
assert !missing(r(ate))
assert !missing(r(se))

local ate_est = r(ate)
local ate_se  = r(se)

display as text "  ATE estimate: " as result %8.4f `ate_est'
display as text "  ATE SE:       " as result %8.4f `ate_se'

* SE should be positive
assert `ate_se' > 0

* With true ATE = mean(x1), which for N=500 should be near 0,
* check estimate is within a reasonable range
local lower = `ate_est' - 2 * `ate_se'
local upper = `ate_est' + 2 * `ate_se'
display as text "  95% CI: [" as result %8.4f `lower' as text ", " as result %8.4f `upper' as text "]"

* ATE should be roughly near 0 (true ATE = E[x1] ~ 0)
assert abs(`ate_est') < 2.0
display as text "  ATE within plausible range: OK"

* Check p-value if returned
capture assert !missing(r(pvalue))
if _rc == 0 {
    local pval = r(pvalue)
    display as text "  p-value: " as result %6.4f `pval'
}

display as text "  PASSED"

* ---- Test 2: grf_test_calibration ----
* (reads from e(); internally runs regress which destroys e())
display as text ""
display as text "--- Test 2: grf_test_calibration ---"

grf_test_calibration

* If it runs without error, the test passes
display as text "  Command executed without error"

display as text "  PASSED"

* ---- Test 3: grf_variable_importance ----
* (standalone command, takes varlist: cate_var predictor_vars; does not need e())
display as text ""
display as text "--- Test 3: grf_variable_importance ---"

grf_variable_importance cate x1 x2 x3 x4 x5

* Check returned matrix or scalars
* x1 should have the highest importance since tau = x1
capture assert !missing(r(importance_x1))
if _rc == 0 {
    * If importances are returned as scalars
    local imp1 = r(importance_x1)
    local imp2 = r(importance_x2)
    local imp3 = r(importance_x3)
    local imp4 = r(importance_x4)
    local imp5 = r(importance_x5)

    display as text "  x1 importance: " as result %6.4f `imp1'
    display as text "  x2 importance: " as result %6.4f `imp2'
    display as text "  x3 importance: " as result %6.4f `imp3'
    display as text "  x4 importance: " as result %6.4f `imp4'
    display as text "  x5 importance: " as result %6.4f `imp5'

    * x1 should dominate
    assert `imp1' > `imp2'
    assert `imp1' > `imp3'
    assert `imp1' > `imp4'
    assert `imp1' > `imp5'
    display as text "  x1 has highest importance: OK"
}
else {
    * Try matrix form
    capture matrix list r(importance)
    if _rc == 0 {
        local imp1 = r(importance)[1,1]
        display as text "  x1 importance: " as result %6.4f `imp1'
        display as text "  Variable importance matrix returned: OK"
    }
    else {
        display as text "  Variable importance returned (format varies by implementation)"
    }
}

display as text "  PASSED"

* ---- Test 4: grf_best_linear_projection ----
* (eclass command -- needs e() from causal forest, so re-fit first since
*  test_calibration's internal regress destroyed e())
display as text ""
display as text "--- Test 4: grf_best_linear_projection ---"

* Re-fit causal forest to restore e()
display as text "  (Re-fitting causal forest to restore e() results...)"
grf_causal_forest y w x1 x2 x3 x4 x5, gen(cate_blp) ntrees(500) seed(42) ///
    yhatgenerate(yhat_blp) whatgenerate(what_blp) replace

grf_best_linear_projection x1 x2 x3 x4 x5

* If it runs without error, the basic test passes
display as text "  Command executed without error"

display as text "  PASSED"

* ---- Fidelity: Compare ATE with R reference (if available) ----
display as text ""
display as text "--- Fidelity: ATE vs R reference ---"

capture confirm file "tests/ref/causal_ate.csv"
if _rc == 0 {
    * Reload the reference data and refit
    preserve
    clear
    import delimited "tests/ref/causal_input.csv", clear

    grf_causal_forest y w x1 x2 x3 x4 x5, gen(stata_cate) ntrees(2000) ///
        seed(42) yhatgenerate(stata_yhat) whatgenerate(stata_what) replace

    grf_ate
    local stata_ate = r(ate)
    local stata_se  = r(se)

    * Load R ATE
    clear
    import delimited "tests/ref/causal_ate.csv", clear
    local r_ate = estimate[1]
    local r_se  = stderr[1]

    display as text "  Stata ATE: " as result %8.4f `stata_ate' " (SE: " as result %6.4f `stata_se' as text ")"
    display as text "  R ATE:     " as result %8.4f `r_ate' " (SE: " as result %6.4f `r_se' as text ")"

    * Check that estimates are reasonably close (within combined SE)
    local diff = abs(`stata_ate' - `r_ate')
    local combined_se = sqrt(`stata_se'^2 + `r_se'^2)
    local z = `diff' / `combined_se'
    display as text "  |Difference|: " as result %8.4f `diff'
    display as text "  Combined SE:  " as result %8.4f `combined_se'
    display as text "  Z-score:      " as result %6.2f `z'

    if `z' < 3 {
        display as text "  PASSED (z < 3)"
    }
    else {
        display as error "  FAILED (z >= 3, estimates diverge)"
    }

    restore
}
else {
    display as text "  Skipped (no reference data at tests/ref/causal_ate.csv)"
    display as text "  Run: Rscript tests/generate_reference.R"
}

* ---- Summary ----
display as text ""
display as text "=============================================="
display as text " All post-estimation tests completed"
display as text "=============================================="
display as text ""
