* test_regression.do -- Test grf_regression_forest against R reference data
* Run from the test-c-plugin-skill-grf/ directory

clear all
set more off

display as text ""
display as text "=============================================="
display as text " GRF Regression Forest Tests"
display as text "=============================================="

* ---- Test 1: Basic functionality (synthetic data) ----
display as text ""
display as text "--- Test 1: Basic functionality ---"

clear
set obs 500
set seed 42

* Generate test data
gen x1 = rnormal()
gen x2 = rnormal()
gen x3 = rnormal()
gen x4 = rnormal()
gen x5 = rnormal()
gen y = 2 * x1 + x2^2 + 0.5 * rnormal()

* Run regression forest
grf_regression_forest y x1 x2 x3 x4 x5, gen(pred) ntrees(500) seed(42) replace

* Check predictions exist
quietly count if !missing(pred)
display as text "  Predictions: " as result r(N) " / 500"
assert r(N) > 400

* Check predictions are reasonable (correlated with Y)
corr y pred
assert r(rho) > 0.5
display as text "  Correlation(Y, pred): " as result %6.4f r(rho) " (> 0.50 required)"

* Check prediction range is reasonable
summarize pred
assert r(sd) > 0.1
display as text "  Prediction SD: " as result %6.4f r(sd)

display as text "  PASSED"

* ---- Test 2: Variance estimation ----
display as text ""
display as text "--- Test 2: Variance estimation ---"

grf_regression_forest y x1 x2 x3 x4 x5, gen(pred2) ntrees(500) seed(42) ///
    estimatevariance vargen(pred2_var) replace

quietly count if !missing(pred2_var)
display as text "  Variance estimates: " as result r(N) " / 500"
assert r(N) > 400

* Variance should be positive
summarize pred2_var
assert r(min) >= 0
display as text "  Min variance: " as result %9.6f r(min) " (>= 0 required)"
display as text "  Mean variance: " as result %9.6f r(mean)

display as text "  PASSED"

* ---- Test 3: if/in restrictions ----
display as text ""
display as text "--- Test 3: if/in restrictions ---"

grf_regression_forest y x1 x2 x3 x4 x5 if x1 > 0, gen(pred3) ntrees(200) seed(42) replace

quietly count if !missing(pred3) & x1 > 0
local n_pred = r(N)
quietly count if x1 > 0
local n_eligible = r(N)
display as text "  Predictions: " as result `n_pred' " / " as result `n_eligible' " (if x1 > 0)"
assert `n_pred' > 0

* Predictions should only exist for x1 > 0
quietly count if !missing(pred3) & x1 <= 0
assert r(N) == 0
display as text "  No predictions for x1 <= 0: OK"

display as text "  PASSED"

* ---- Test 4: Options ----
display as text ""
display as text "--- Test 4: Various options ---"

* No honesty
grf_regression_forest y x1 x2 x3 x4 x5, gen(pred4a) ntrees(200) seed(42) ///
    nohonesty replace
quietly count if !missing(pred4a)
assert r(N) > 400
display as text "  No honesty: " as result r(N) " predictions - PASSED"

* Different sample fraction
grf_regression_forest y x1 x2 x3 x4 x5, gen(pred4b) ntrees(200) seed(42) ///
    samplefrac(0.7) replace
quietly count if !missing(pred4b)
assert r(N) > 400
display as text "  Sample frac 0.7: " as result r(N) " predictions - PASSED"

* Different min node size
grf_regression_forest y x1 x2 x3 x4 x5, gen(pred4c) ntrees(200) seed(42) ///
    minnodesize(20) replace
quietly count if !missing(pred4c)
assert r(N) > 400
display as text "  Min node 20: " as result r(N) " predictions - PASSED"

display as text "  ALL OPTIONS TESTS PASSED"

* ---- Test 5: Fidelity vs R reference (if available) ----
display as text ""
display as text "--- Test 5: Fidelity vs R reference ---"

capture confirm file "tests/ref/regression_input.csv"
if _rc == 0 {
    clear
    import delimited "tests/ref/regression_input.csv", clear

    * Run forest on same data
    grf_regression_forest y x1 x2 x3 x4 x5, gen(stata_pred) ntrees(2000) ///
        seed(42) estimatevariance vargen(stata_var) replace

    * Load R predictions
    preserve
    import delimited "tests/ref/regression_output.csv", clear
    rename prediction r_pred
    rename variance r_var
    gen n = _n
    tempfile rref
    save `rref'
    restore

    gen n = _n
    merge 1:1 n using `rref', nogenerate

    * Compare predictions
    corr stata_pred r_pred
    local r_corr = r(rho)
    display as text "  Stata-R prediction correlation: " as result %6.4f `r_corr'

    if `r_corr' >= 0.95 {
        display as text "  PASSED (r >= 0.95)"
    }
    else if `r_corr' >= 0.90 {
        display as text "  MARGINAL (0.90 <= r < 0.95)"
    }
    else {
        display as error "  FAILED (r < 0.90)"
    }
}
else {
    display as text "  Skipped (no reference data at tests/ref/regression_input.csv)"
    display as text "  Run: Rscript tests/generate_reference.R"
}

* ---- Summary ----
display as text ""
display as text "=============================================="
display as text " All regression forest tests completed"
display as text "=============================================="
display as text ""
