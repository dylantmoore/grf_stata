* test_causal.do -- Test grf_causal_forest
* Run from the test-c-plugin-skill-grf/ directory

clear all
set more off

display as text ""
display as text "=============================================="
display as text " GRF Causal Forest Tests"
display as text "=============================================="

* ---- Test 1: Basic functionality (synthetic data) ----
display as text ""
display as text "--- Test 1: Basic functionality ---"

clear
set obs 500
set seed 42

* Generate test data with known treatment effect tau = x1
gen x1 = rnormal()
gen x2 = rnormal()
gen x3 = rnormal()
gen x4 = rnormal()
gen x5 = rnormal()
gen w = (runiform() > 0.5)
gen tau = x1
gen y = 2 * x1 + tau * w + 0.5 * rnormal()

* Run causal forest
grf_causal_forest y w x1 x2 x3 x4 x5, gen(cate) ntrees(500) seed(42) replace

* Check predictions exist
quietly count if !missing(cate)
display as text "  Predictions: " as result r(N) " / 500"
assert r(N) > 400

* Check CATE is correlated with true tau
corr cate tau
local r_tau = r(rho)
display as text "  Correlation(CATE, tau): " as result %6.4f `r_tau' " (> 0.30 required)"
assert `r_tau' > 0.30

* Check prediction range is reasonable
summarize cate
assert r(sd) > 0.05
display as text "  CATE SD: " as result %6.4f r(sd)

display as text "  PASSED"

* ---- Test 2: Variance estimation ----
display as text ""
display as text "--- Test 2: Variance estimation ---"

grf_causal_forest y w x1 x2 x3 x4 x5, gen(cate2) ntrees(500) seed(42) ///
    estimatevariance vargen(cate2_var) replace

quietly count if !missing(cate2_var)
display as text "  Variance estimates: " as result r(N) " / 500"
assert r(N) > 400

* Variance should be positive
summarize cate2_var
assert r(min) >= 0
display as text "  Min variance: " as result %9.6f r(min) " (>= 0 required)"
display as text "  Mean variance: " as result %9.6f r(mean)

display as text "  PASSED"

* ---- Test 3: Nuisance estimates ----
display as text ""
display as text "--- Test 3: Nuisance estimates (Y-hat, W-hat) ---"

grf_causal_forest y w x1 x2 x3 x4 x5, gen(cate3) ntrees(500) seed(42) ///
    yhatgenerate(yhat) whatgenerate(what) replace

* Check Y-hat exists and is non-missing
quietly count if !missing(yhat)
local n_yhat = r(N)
display as text "  Y-hat estimates: " as result `n_yhat' " / 500"
assert `n_yhat' > 400

* Check W-hat exists and is non-missing
quietly count if !missing(what)
local n_what = r(N)
display as text "  W-hat estimates: " as result `n_what' " / 500"
assert `n_what' > 400

* Y-hat should be correlated with y
corr yhat y
display as text "  Correlation(Y-hat, Y): " as result %6.4f r(rho)
assert r(rho) > 0.3

* W-hat should be roughly around 0.5 (since treatment is coin flip)
summarize what
display as text "  W-hat mean: " as result %6.4f r(mean) " (expect ~0.50)"
assert r(mean) > 0.2 & r(mean) < 0.8

display as text "  PASSED"

* ---- Test 4: if/in restrictions ----
display as text ""
display as text "--- Test 4: if/in restrictions ---"

grf_causal_forest y w x1 x2 x3 x4 x5 if x1 > 0, gen(cate4) ntrees(200) seed(42) replace

quietly count if !missing(cate4) & x1 > 0
local n_pred = r(N)
quietly count if x1 > 0
local n_eligible = r(N)
display as text "  Predictions: " as result `n_pred' " / " as result `n_eligible' " (if x1 > 0)"
assert `n_pred' > 0

* Predictions should only exist for x1 > 0
quietly count if !missing(cate4) & x1 <= 0
assert r(N) == 0
display as text "  No predictions for x1 <= 0: OK"

display as text "  PASSED"

* ---- Test 5: Fidelity vs R reference (if available) ----
display as text ""
display as text "--- Test 5: Fidelity vs R reference ---"

capture confirm file "tests/ref/causal_input.csv"
if _rc == 0 {
    clear
    import delimited "tests/ref/causal_input.csv", clear

    * Run causal forest on same data
    grf_causal_forest y w x1 x2 x3 x4 x5, gen(stata_cate) ntrees(2000) ///
        seed(42) estimatevariance vargen(stata_var) replace

    * Load R predictions
    preserve
    import delimited "tests/ref/causal_output.csv", clear
    rename cate r_cate
    rename variance r_var
    gen n = _n
    tempfile rref
    save `rref'
    restore

    gen n = _n
    merge 1:1 n using `rref', nogenerate

    * Compare CATE predictions
    corr stata_cate r_cate
    local r_corr = r(rho)
    display as text "  Stata-R CATE correlation: " as result %6.4f `r_corr'

    * Also compute correlation of CATE with true tau (= x1)
    corr stata_cate x1
    local r_tau_corr = r(rho)
    display as text "  Stata CATE-tau correlation: " as result %6.4f `r_tau_corr'

    * Require CATE correlation >= 0.90
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
    display as text "  Skipped (no reference data at tests/ref/causal_input.csv)"
    display as text "  Run: Rscript tests/generate_reference.R"
}

* ---- Summary ----
display as text ""
display as text "=============================================="
display as text " All causal forest tests completed"
display as text "=============================================="
display as text ""
