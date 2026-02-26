* test_survival.do -- Test grf_survival_forest against R reference data
* Run from the test-c-plugin-skill-grf/ directory

clear all
set more off

display as text ""
display as text "=============================================="
display as text " GRF Survival Forest Tests"
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
gen double time = exp(0.5 * x1 + rnormal())
gen status = (runiform() > 0.3)

* Run survival forest
grf_survival_forest time status x1 x2 x3, gen(surv) ntrees(500) seed(42) noutput(10) replace

* Check that surv_s1 through surv_s10 exist and have non-missing values
forvalues i = 1/10 {
    quietly count if !missing(surv_s`i')
    display as text "  surv_s`i': " as result r(N) " / 500 non-missing"
    assert r(N) > 400

    * Survival probabilities must be in [0, 1]
    summarize surv_s`i', meanonly
    assert r(min) >= 0
    assert r(max) <= 1
}

display as text "  All 10 survival columns exist with values in [0, 1]"
display as text "  PASSED"

* ---- Test 2: Kaplan-Meier prediction type ----
display as text ""
display as text "--- Test 2: Kaplan-Meier prediction type (predtype=1) ---"

grf_survival_forest time status x1 x2 x3, gen(km) ntrees(500) seed(42) ///
    noutput(10) predtype(1) replace

* Check predictions exist
forvalues i = 1/10 {
    quietly count if !missing(km_s`i')
    assert r(N) > 400
}

* Kaplan-Meier survival probabilities should also be in [0, 1]
forvalues i = 1/10 {
    summarize km_s`i', meanonly
    assert r(min) >= 0
    assert r(max) <= 1
}

display as text "  All 10 KM columns exist with values in [0, 1]"
display as text "  PASSED"

* ---- Test 3: if/in restrictions ----
display as text ""
display as text "--- Test 3: if/in restrictions ---"

grf_survival_forest time status x1 x2 x3 if x1 > 0, gen(surv3) ntrees(200) ///
    seed(42) noutput(5) replace

quietly count if !missing(surv3_s1) & x1 > 0
local n_pred = r(N)
quietly count if x1 > 0
local n_eligible = r(N)
display as text "  Predictions: " as result `n_pred' " / " as result `n_eligible' " (if x1 > 0)"
assert `n_pred' > 0

* Predictions should only exist for x1 > 0
quietly count if !missing(surv3_s1) & x1 <= 0
assert r(N) == 0
display as text "  No predictions for x1 <= 0: OK"

* Check all 5 output columns for the restricted sample
forvalues i = 1/5 {
    summarize surv3_s`i' if x1 > 0, meanonly
    assert r(min) >= 0
    assert r(max) <= 1
}

display as text "  All 5 columns valid in [0, 1] for restricted sample"
display as text "  PASSED"

* ---- Test 4: Fidelity vs R reference (if available) ----
display as text ""
display as text "--- Test 4: Fidelity vs R reference ---"

capture confirm file "tests/ref/survival_input.csv"
if _rc == 0 {
    clear
    import delimited "tests/ref/survival_input.csv", clear

    * Run forest on same data (R uses x1-x5)
    grf_survival_forest time status x1 x2 x3 x4 x5, gen(stata_surv) ntrees(2000) ///
        seed(42) noutput(10) replace

    * Load R predictions
    preserve
    import delimited "tests/ref/survival_output.csv", clear
    * R generates columns named t1, t2, ..., t10 (failure time columns)
    gen n = _n
    tempfile rref
    save `rref'
    restore

    gen n = _n
    merge 1:1 n using `rref', nogenerate

    * Compare survival curves column by column
    * R columns are t1..t10, Stata columns are stata_surv_s1..stata_surv_s10
    local all_pass = 1
    forvalues i = 1/10 {
        capture confirm variable t`i'
        if _rc == 0 {
            corr stata_surv_s`i' t`i'
            local r_corr = r(rho)
            display as text "  Column s`i' Stata-R correlation: " as result %6.4f `r_corr'

            if `r_corr' < 0.90 {
                local all_pass = 0
            }
        }
    }

    if `all_pass' == 1 {
        display as text "  PASSED (all columns r >= 0.90)"
    }
    else {
        display as error "  FAILED (one or more columns r < 0.90)"
    }
}
else {
    display as text "  Skipped (no reference data at tests/ref/survival_input.csv)"
    display as text "  Run: Rscript tests/generate_reference.R"
}

* ---- Summary ----
display as text ""
display as text "=============================================="
display as text " All survival forest tests completed"
display as text "=============================================="
display as text ""
