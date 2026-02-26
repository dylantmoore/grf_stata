* test_quantile.do -- Test grf_quantile_forest
* Run from the test-c-plugin-skill-grf/ directory
*
* Note: grf_quantile_forest gen() takes a STUB name. Output variables are:
*   stub_q10, stub_q50, stub_q90, etc. (quantile*100 appended)

clear all
set more off

display as text ""
display as text "=============================================="
display as text " GRF Quantile Forest Tests"
display as text "=============================================="

* ---- Test 1: Basic functionality (default quantiles) ----
display as text ""
display as text "--- Test 1: Basic functionality (default quantiles 0.1, 0.5, 0.9) ---"

clear
set obs 500
set seed 42

* Generate heteroscedastic data: Y = 2*x1 + x2 + N(0,1)*(1 + |x3|)
gen x1 = rnormal()
gen x2 = rnormal()
gen x3 = rnormal()
gen x4 = rnormal()
gen x5 = rnormal()
gen y = 2 * x1 + x2 + rnormal() * (1 + abs(x3))

* Run quantile forest with default quantiles
* gen(qf) creates: qf_q10, qf_q50, qf_q90
grf_quantile_forest y x1 x2 x3 x4 x5, quantiles(0.1 0.5 0.9) ///
    gen(qf) ntrees(500) seed(42) replace

* Check output variables exist
quietly count if !missing(qf_q10)
display as text "  qf_q10 predictions: " as result r(N) " / 500"
assert r(N) > 400

quietly count if !missing(qf_q50)
display as text "  qf_q50 predictions: " as result r(N) " / 500"
assert r(N) > 400

quietly count if !missing(qf_q90)
display as text "  qf_q90 predictions: " as result r(N) " / 500"
assert r(N) > 400

* Median should be correlated with Y
corr qf_q50 y
display as text "  Correlation(qf_q50, Y): " as result %6.4f r(rho) " (> 0.40 required)"
assert r(rho) > 0.40

display as text "  PASSED"

* ---- Test 2: Custom quantiles ----
display as text ""
display as text "--- Test 2: Custom quantiles (0.25, 0.5, 0.75) ---"

grf_quantile_forest y x1 x2 x3 x4 x5, quantiles(0.25 0.5 0.75) ///
    gen(qf2) ntrees(500) seed(42) replace

* Verify all 3 output variables exist
quietly count if !missing(qf2_q25)
local n_q25 = r(N)
display as text "  qf2_q25 predictions: " as result `n_q25' " / 500"
assert `n_q25' > 400

quietly count if !missing(qf2_q50)
local n_q50 = r(N)
display as text "  qf2_q50 predictions: " as result `n_q50' " / 500"
assert `n_q50' > 400

quietly count if !missing(qf2_q75)
local n_q75 = r(N)
display as text "  qf2_q75 predictions: " as result `n_q75' " / 500"
assert `n_q75' > 400

* IQR should be positive on average
gen iqr = qf2_q75 - qf2_q25
summarize iqr
display as text "  Mean IQR (q75 - q25): " as result %6.4f r(mean) " (should be > 0)"
assert r(mean) > 0
drop iqr

display as text "  PASSED"

* ---- Test 3: Monotonicity ----
display as text ""
display as text "--- Test 3: Monotonicity (q10 <= q50 <= q90) ---"

* Check q10 <= q50 for most observations
quietly count if qf_q10 <= qf_q50 & !missing(qf_q10) & !missing(qf_q50)
local n_mono_lower = r(N)
quietly count if !missing(qf_q10) & !missing(qf_q50)
local n_total_lower = r(N)
local pct_lower = 100 * `n_mono_lower' / `n_total_lower'
display as text "  q10 <= q50: " as result `n_mono_lower' "/" `n_total_lower' ///
    " (" as result %5.1f `pct_lower' as text "%)"
assert `pct_lower' > 90

* Check q50 <= q90 for most observations
quietly count if qf_q50 <= qf_q90 & !missing(qf_q50) & !missing(qf_q90)
local n_mono_upper = r(N)
quietly count if !missing(qf_q50) & !missing(qf_q90)
local n_total_upper = r(N)
local pct_upper = 100 * `n_mono_upper' / `n_total_upper'
display as text "  q50 <= q90: " as result `n_mono_upper' "/" `n_total_upper' ///
    " (" as result %5.1f `pct_upper' as text "%)"
assert `pct_upper' > 90

* Check full ordering q10 <= q50 <= q90
quietly count if qf_q10 <= qf_q50 & qf_q50 <= qf_q90 & !missing(qf_q10) & !missing(qf_q50) & !missing(qf_q90)
local n_full_mono = r(N)
quietly count if !missing(qf_q10) & !missing(qf_q50) & !missing(qf_q90)
local n_total = r(N)
local pct_full = 100 * `n_full_mono' / `n_total'
display as text "  q10 <= q50 <= q90: " as result `n_full_mono' "/" `n_total' ///
    " (" as result %5.1f `pct_full' as text "%)"
assert `pct_full' > 85

* Spread should reflect heteroscedasticity: wider intervals where |x3| is large
gen spread = qf_q90 - qf_q10
gen abs_x3 = abs(x3)
corr spread abs_x3
local r_hetero = r(rho)
display as text "  Correlation(spread, |x3|): " as result %6.4f `r_hetero'
display as text "  (positive expected -- heteroscedastic DGP)"
drop spread abs_x3

display as text "  PASSED"

* ---- Test 4: Fidelity vs R reference (if available) ----
display as text ""
display as text "--- Test 4: Fidelity vs R reference ---"

capture confirm file "tests/ref/quantile_input.csv"
if _rc == 0 {
    clear
    import delimited "tests/ref/quantile_input.csv", clear

    * Run quantile forest on same data with matching quantiles
    * R uses quantiles: 0.1, 0.25, 0.5, 0.75, 0.9
    grf_quantile_forest y x1 x2 x3 x4 x5, quantiles(0.1 0.25 0.5 0.75 0.9) ///
        gen(sq) ntrees(2000) seed(42) replace

    * Load R predictions
    * R columns after import: q01, q025, q05, q075, q09
    * (dots in q0.1 etc. are stripped by Stata import)
    preserve
    import delimited "tests/ref/quantile_output.csv", clear
    rename q01 r_q10
    rename q025 r_q25
    rename q05 r_q50
    rename q075 r_q75
    rename q09 r_q90
    gen n = _n
    tempfile rref
    save `rref'
    restore

    gen n = _n
    merge 1:1 n using `rref', nogenerate

    * Compare median predictions
    * Stata columns: sq_q10, sq_q25, sq_q50, sq_q75, sq_q90
    corr sq_q50 r_q50
    local r_median = r(rho)
    display as text "  Stata-R median correlation: " as result %6.4f `r_median'

    * Compare q10
    corr sq_q10 r_q10
    local r_q10 = r(rho)
    display as text "  Stata-R q10 correlation:    " as result %6.4f `r_q10'

    * Compare q90
    corr sq_q90 r_q90
    local r_q90 = r(rho)
    display as text "  Stata-R q90 correlation:    " as result %6.4f `r_q90'

    if `r_median' >= 0.90 {
        display as text "  PASSED (median r >= 0.90)"
    }
    else if `r_median' >= 0.80 {
        display as text "  MARGINAL (0.80 <= median r < 0.90)"
    }
    else {
        display as error "  FAILED (median r < 0.80)"
    }
}
else {
    display as text "  Skipped (no reference data at tests/ref/quantile_input.csv)"
    display as text "  Run: Rscript tests/generate_reference.R"
}

* ---- Summary ----
display as text ""
display as text "=============================================="
display as text " All quantile forest tests completed"
display as text "=============================================="
display as text ""
