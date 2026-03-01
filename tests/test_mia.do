* test_mia.do -- Tests for Missing Indicator Action (MIA) passthrough
*
* Tests that the grf plugin correctly handles missing covariates via
* grf's native MIA splitting, and that the nomia option reverts to
* casewise deletion.
*
* Run: do tests/test_mia.do

clear all
set more off

local errors = 0

* ============================================================
* 1. Regression forest with missing covariates (MIA enabled)
* ============================================================
display as text ""
display as text "=== MIA Test 1: Regression forest with missing X1 ==="

capture noisily {
    clear
    set obs 500
    set seed 42
    gen x1 = rnormal()
    gen x2 = rnormal()
    gen x3 = rnormal()
    gen y = 2*x1 + x2 + 0.5*rnormal()

    * Introduce 10% missing in x1
    replace x1 = . if runiform() < 0.10

    * Count non-missing and total
    count
    local n_total = r(N)
    count if !missing(x1)
    local n_complete = r(N)

    * MIA mode (default) -- should use all N observations
    grf_regression_forest y x1 x2 x3, gen(pred_mia) ntrees(500) seed(42) replace
    local n_used = e(N)

    * With MIA, all observations should be used
    assert `n_used' == `n_total'

    * Predictions should exist for all observations
    count if !missing(pred_mia)
    assert r(N) == `n_total'
}
if _rc {
    display as error "FAIL: regression forest with MIA"
    local errors = `errors' + 1
}
else {
    display as result "PASS: regression forest with MIA (N=`n_used' == `n_total')"
}

* ============================================================
* 2. Regression forest with nomia (casewise deletion)
* ============================================================
display as text ""
display as text "=== MIA Test 2: Regression forest with nomia ==="

capture noisily {
    * nomia option -- should drop observations with missing X1
    grf_regression_forest y x1 x2 x3, gen(pred_nomia) ntrees(500) seed(42) nomia replace
    local n_nomia = e(N)

    * With nomia, N should equal the number of complete cases
    assert `n_nomia' == `n_complete'
    assert `n_nomia' < `n_total'
}
if _rc {
    display as error "FAIL: regression forest with nomia"
    local errors = `errors' + 1
}
else {
    display as result "PASS: regression forest with nomia (N=`n_nomia' == `n_complete')"
}

* ============================================================
* 3. Causal forest with missing covariates (MIA enabled)
* ============================================================
display as text ""
display as text "=== MIA Test 3: Causal forest with missing X1 ==="

capture noisily {
    clear
    set obs 500
    set seed 42
    gen x1 = rnormal()
    gen x2 = rnormal()
    gen w = rbinomial(1, 0.5)
    gen y = 2*x1 + x1*w + 0.5*rnormal()

    * Introduce 10% missing in x1
    replace x1 = . if runiform() < 0.10

    count
    local n_total = r(N)

    * MIA mode (default)
    grf_causal_forest y w x1 x2, gen(cate_mia) ntrees(500) seed(42) replace
    local n_used = e(N)

    assert `n_used' == `n_total'
}
if _rc {
    display as error "FAIL: causal forest with MIA"
    local errors = `errors' + 1
}
else {
    display as result "PASS: causal forest with MIA (N=`n_used' == `n_total')"
}

* ============================================================
* 4. Causal forest with nomia
* ============================================================
display as text ""
display as text "=== MIA Test 4: Causal forest with nomia ==="

capture noisily {
    count if !missing(x1)
    local n_complete = r(N)

    grf_causal_forest y w x1 x2, gen(cate_nomia) ntrees(500) seed(42) nomia replace
    local n_nomia = e(N)

    assert `n_nomia' == `n_complete'
    assert `n_nomia' < `n_total'
}
if _rc {
    display as error "FAIL: causal forest with nomia"
    local errors = `errors' + 1
}
else {
    display as result "PASS: causal forest with nomia (N=`n_nomia' == `n_complete')"
}

* ============================================================
* 5. Missing outcome (Y) should always drop observation
* ============================================================
display as text ""
display as text "=== MIA Test 5: Missing Y always drops ==="

capture noisily {
    clear
    set obs 500
    set seed 42
    gen x1 = rnormal()
    gen x2 = rnormal()
    gen y = 2*x1 + x2 + 0.5*rnormal()

    * Make 5% of Y missing
    replace y = . if runiform() < 0.05

    count if !missing(y)
    local n_complete_y = r(N)
    count
    local n_total = r(N)

    * Even with MIA (default), missing Y should be dropped
    grf_regression_forest y x1 x2, gen(pred_missy) ntrees(500) seed(42) replace
    local n_used = e(N)

    assert `n_used' == `n_complete_y'
    assert `n_used' < `n_total'
}
if _rc {
    display as error "FAIL: missing Y always drops observation"
    local errors = `errors' + 1
}
else {
    display as result "PASS: missing Y always drops (N=`n_used' == `n_complete_y')"
}

* ============================================================
* 6. Predictions from MIA vs nomia should differ
* ============================================================
display as text ""
display as text "=== MIA Test 6: MIA vs nomia predictions differ ==="

capture noisily {
    clear
    set obs 500
    set seed 42
    gen x1 = rnormal()
    gen x2 = rnormal()
    gen y = 2*x1 + x2 + 0.5*rnormal()
    replace x1 = . if runiform() < 0.15

    grf_regression_forest y x1 x2, gen(p_mia) ntrees(500) seed(42) replace
    grf_regression_forest y x1 x2, gen(p_nomia) ntrees(500) seed(42) nomia replace

    * Compare predictions on complete cases only
    gen diff = p_mia - p_nomia if !missing(p_mia) & !missing(p_nomia)
    summarize diff
    * With different training sets, predictions should differ somewhat
    assert r(sd) > 0
}
if _rc {
    display as error "FAIL: MIA vs nomia predictions differ"
    local errors = `errors' + 1
}
else {
    display as result "PASS: MIA vs nomia predictions differ"
}

* ============================================================
* Summary
* ============================================================
display ""
if `errors' > 0 {
    display as error "FAILED: `errors' MIA test(s) did not pass"
    exit 1
}
else {
    display as result "ALL MIA TESTS PASSED"
}
