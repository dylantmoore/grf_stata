* test_probability.do -- Test grf_probability_forest against R reference data
* Run from the test-c-plugin-skill-grf/ directory

clear all
set more off

display as text ""
display as text "=============================================="
display as text " GRF Probability Forest Tests"
display as text "=============================================="

* ---- Test 1: Basic functionality (2 classes) ----
display as text ""
display as text "--- Test 1: Basic functionality (2 classes) ---"

clear
set obs 500
set seed 42

* Generate test data
gen x1 = rnormal()
gen x2 = rnormal()
gen x3 = rnormal()
gen x4 = rnormal()
gen y = (x1 + x2 > 0)

* Run probability forest
grf_probability_forest y x1 x2 x3 x4, gen(prob) ntrees(500) seed(42) replace

* Check class probability variables exist
quietly count if !missing(prob_c0)
display as text "  prob_c0 predictions: " as result r(N) " / 500"
assert r(N) > 400

quietly count if !missing(prob_c1)
display as text "  prob_c1 predictions: " as result r(N) " / 500"
assert r(N) > 400

* Check probabilities are between 0 and 1
summarize prob_c0
assert r(min) >= 0 & r(max) <= 1
display as text "  prob_c0 range: [" as result %6.4f r(min) ", " as result %6.4f r(max) "]"

summarize prob_c1
assert r(min) >= 0 & r(max) <= 1
display as text "  prob_c1 range: [" as result %6.4f r(min) ", " as result %6.4f r(max) "]"

* Check probabilities sum to approximately 1
gen prob_sum = prob_c0 + prob_c1
summarize prob_sum
display as text "  Prob sum mean: " as result %9.6f r(mean) " (should be ~1.0)"
assert abs(r(mean) - 1) < 0.01
assert abs(r(min) - 1) < 0.05
assert abs(r(max) - 1) < 0.05
drop prob_sum

display as text "  PASSED"

* ---- Test 2: Three classes ----
display as text ""
display as text "--- Test 2: Three classes ---"

clear
set obs 600
set seed 42

gen x1 = rnormal()
gen x2 = rnormal()
gen x3 = rnormal()
gen x4 = rnormal()

* Create 3-class outcome based on x1 thresholds
gen y = 0
replace y = 1 if x1 >= -0.5 & x1 < 0.5
replace y = 2 if x1 >= 0.5

* Run probability forest
grf_probability_forest y x1 x2 x3 x4, gen(prob3) ntrees(500) seed(42) replace

* Check all three class variables exist
quietly count if !missing(prob3_c0)
display as text "  prob3_c0 predictions: " as result r(N) " / 600"
assert r(N) > 400

quietly count if !missing(prob3_c1)
display as text "  prob3_c1 predictions: " as result r(N) " / 600"
assert r(N) > 400

quietly count if !missing(prob3_c2)
display as text "  prob3_c2 predictions: " as result r(N) " / 600"
assert r(N) > 400

* Check all probabilities are non-negative
summarize prob3_c0
assert r(min) >= 0 & r(max) <= 1
display as text "  prob3_c0 range: [" as result %6.4f r(min) ", " as result %6.4f r(max) "]"

summarize prob3_c1
assert r(min) >= 0 & r(max) <= 1
display as text "  prob3_c1 range: [" as result %6.4f r(min) ", " as result %6.4f r(max) "]"

summarize prob3_c2
assert r(min) >= 0 & r(max) <= 1
display as text "  prob3_c2 range: [" as result %6.4f r(min) ", " as result %6.4f r(max) "]"

* Check probabilities sum to approximately 1
gen prob3_sum = prob3_c0 + prob3_c1 + prob3_c2
summarize prob3_sum
display as text "  Prob sum mean: " as result %9.6f r(mean) " (should be ~1.0)"
assert abs(r(mean) - 1) < 0.01
assert abs(r(min) - 1) < 0.05
assert abs(r(max) - 1) < 0.05
drop prob3_sum

display as text "  PASSED"

* ---- Test 3: if/in restrictions ----
display as text ""
display as text "--- Test 3: if/in restrictions ---"

grf_probability_forest y x1 x2 x3 x4 if x1 > 0, gen(prob_if) ntrees(200) seed(42) replace

quietly count if !missing(prob_if_c0) & x1 > 0
local n_pred = r(N)
quietly count if x1 > 0
local n_eligible = r(N)
display as text "  Predictions: " as result `n_pred' " / " as result `n_eligible' " (if x1 > 0)"
assert `n_pred' > 0

* Predictions should only exist for x1 > 0
quietly count if !missing(prob_if_c0) & x1 <= 0
assert r(N) == 0
display as text "  No predictions for x1 <= 0: OK"

display as text "  PASSED"

* ---- Test 4: Fidelity vs R reference (if available) ----
display as text ""
display as text "--- Test 4: Fidelity vs R reference ---"

capture confirm file "tests/ref/probability_input.csv"
if _rc == 0 {
    clear
    import delimited "tests/ref/probability_input.csv", clear

    * Run forest on same data
    grf_probability_forest y x1 x2 x3 x4, gen(stata_prob) ntrees(2000) ///
        seed(42) replace

    * Load R predictions
    preserve
    import delimited "tests/ref/probability_output.csv", clear

    * R generates columns: class0, class1, class2
    rename class0 r_class0
    rename class1 r_class1
    gen n = _n
    tempfile rref
    save `rref'
    restore

    gen n = _n
    merge 1:1 n using `rref', nogenerate

    * Compare class 0 probabilities
    corr stata_prob_c0 r_class0
    local r_corr0 = r(rho)
    display as text "  Stata-R class 0 correlation: " as result %6.4f `r_corr0'

    * Compare class 1 probabilities
    corr stata_prob_c1 r_class1
    local r_corr1 = r(rho)
    display as text "  Stata-R class 1 correlation: " as result %6.4f `r_corr1'

    local r_corr_min = min(`r_corr0', `r_corr1')
    if `r_corr_min' >= 0.95 {
        display as text "  PASSED (min r >= 0.95)"
    }
    else if `r_corr_min' >= 0.90 {
        display as text "  MARGINAL (0.90 <= min r < 0.95)"
    }
    else {
        display as error "  FAILED (min r < 0.90)"
    }
}
else {
    display as text "  Skipped (no reference data at tests/ref/probability_input.csv)"
    display as text "  Run: Rscript tests/generate_reference.R"
}

* ---- Summary ----
display as text ""
display as text "=============================================="
display as text " All probability forest tests completed"
display as text "=============================================="
display as text ""
