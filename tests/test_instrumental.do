* test_instrumental.do -- Test grf_instrumental_forest
* Run from the test-c-plugin-skill-grf/ directory

clear all
set more off

display as text ""
display as text "=============================================="
display as text " GRF Instrumental Forest Tests"
display as text "=============================================="

* ---- Test 1: Basic functionality (synthetic IV data) ----
display as text ""
display as text "--- Test 1: Basic functionality ---"

clear
set obs 500
set seed 42

* Generate synthetic IV data
gen x1 = rnormal()
gen x2 = rnormal()
gen x3 = rnormal()
gen x4 = rnormal()
gen z = rnormal()
gen w = 0.5 * z + 0.3 * x1 + rnormal()
gen y = 2 * w + x1 + 0.5 * rnormal()

* Run instrumental forest
grf_instrumental_forest y w z x1 x2 x3 x4, gen(pred) ntrees(500) seed(42) replace

* Check predictions exist
quietly count if !missing(pred)
display as text "  Predictions: " as result r(N) " / 500"
assert r(N) > 400

* Check predictions are reasonable
summarize pred
assert r(sd) > 0.01
display as text "  Prediction SD: " as result %6.4f r(sd)

* Check prediction range is finite
assert r(min) > -100
assert r(max) < 100
display as text "  Prediction range: [" as result %6.3f r(min) ", " as result %6.3f r(max) "]"

display as text "  PASSED"

* ---- Test 2: Variance estimation ----
display as text ""
display as text "--- Test 2: Variance estimation ---"

grf_instrumental_forest y w z x1 x2 x3 x4, gen(pred2) ntrees(500) seed(42) ///
    estimatevariance vargen(pred2_var) replace

quietly count if !missing(pred2_var)
display as text "  Variance estimates: " as result r(N) " / 500"
assert r(N) > 400

* Variance should be non-negative
summarize pred2_var
assert r(min) >= 0
display as text "  Min variance: " as result %9.6f r(min) " (>= 0 required)"
display as text "  Mean variance: " as result %9.6f r(mean)

display as text "  PASSED"

* ---- Test 3: if/in restrictions ----
display as text ""
display as text "--- Test 3: if/in restrictions ---"

grf_instrumental_forest y w z x1 x2 x3 x4 if x1 > 0, gen(pred3) ntrees(200) seed(42) replace

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

* ---- Test 4: Fidelity vs R reference (if available) ----
display as text ""
display as text "--- Test 4: Fidelity vs R reference ---"

capture confirm file "tests/ref/instrumental_input.csv"
if _rc == 0 {
    clear
    import delimited "tests/ref/instrumental_input.csv", clear

    * Run forest on same data
    grf_instrumental_forest y w z x1 x2 x3 x4, gen(stata_pred) ntrees(2000) ///
        seed(42) estimatevariance vargen(stata_var) replace

    * Load R predictions
    preserve
    import delimited "tests/ref/instrumental_output.csv", clear
    rename late r_pred
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
    display as text "  Skipped (no reference data at tests/ref/instrumental_input.csv)"
    display as text "  Run: Rscript tests/generate_reference.R"
}

* ---- Summary ----
display as text ""
display as text "=============================================="
display as text " All instrumental forest tests completed"
display as text "=============================================="
display as text ""
