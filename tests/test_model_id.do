* test_model_id.do -- e(model_id) coverage for fit commands

clear all
set more off

local errors = 0

clear
set obs 300
set seed 42
forvalues j = 1/5 {
    gen x`j' = rnormal()
}
gen y = 2*x1 + x2^2 + rnormal()
gen w = (x1 + rnormal() > 0)
gen z = (rnormal() > 0)
gen status = (runiform() > 0.3)
gen time = exp(x1 + 0.5*x2 + rnormal())

* ---- Test 1: model_id exists and increases across fits ----
capture noisily {
    grf_regression_forest y x1-x5, gen(m1) ntrees(100) seed(42)
    local id1 = e(model_id)
    assert !missing(`id1')

    grf_causal_forest y w x1-x5, gen(m2) ntrees(100) seed(42)
    local id2 = e(model_id)
    assert !missing(`id2')
    assert `id2' > `id1'

    grf_survival_forest time status x1-x5, gen(m3) ntrees(100) seed(42)
    local id3 = e(model_id)
    assert !missing(`id3')
    assert `id3' > `id2'

    drop m1 m2 m3_s*
    capture drop _grf_yhat _grf_what
}
if _rc {
    display as error "FAIL: e(model_id) monotone checks"
    local errors = `errors' + 1
}
else {
    display as result "PASS: e(model_id) monotone checks"
}

* ---- Test 2: utility summaries expose model_id ----
capture noisily {
    grf_regression_forest y x1-x5, gen(m4) ntrees(100) seed(99)
    local id4 = e(model_id)

    grf_forest_summary
    assert r(model_id) == `id4'

    grf_tree_summary, tree(1)
    assert r(model_id) == `id4'

    drop m4
}
if _rc {
    display as error "FAIL: summary model_id exposure"
    local errors = `errors' + 1
}
else {
    display as result "PASS: summary model_id exposure"
}

display ""
if `errors' > 0 {
    display as error "FAILED: `errors' model_id test(s) did not pass"
    exit 1
}
else {
    display as result "ALL MODEL_ID TESTS PASSED"
}
