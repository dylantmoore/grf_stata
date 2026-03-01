* TMLE clustered SE + equalizeclusterweights verification
clear all
set more off
local errors = 0

* E5.6: TMLE with clustered SEs
capture noisily {
    clear
    set obs 200
    set seed 42
    gen cluster_id = ceil(_n / 20)
    gen x1 = rnormal()
    gen x2 = rnormal()
    gen w = rbinomial(1, normal(0.3*x1))
    gen y = 2*x1 + x1*w + rnormal()

    grf_causal_forest y w x1 x2, gen(tau_cl) ntrees(100) seed(42) ///
        cluster(cluster_id)
    grf_ate, method(TMLE)
    assert !missing(r(ate))
    assert !missing(r(se))
    assert r(se) > 0
    display as text "  ATE = " r(ate) "  SE = " r(se)
    drop tau_cl _grf_yhat _grf_what
}
if _rc {
    display as error "FAIL: E5.6 TMLE clustered SEs"
    local errors = `errors' + 1
}
else {
    display as result "PASS: E5.6 TMLE clustered SEs"
}

* E4.1: equalize.cluster.weights regression forest
capture noisily {
    clear
    set obs 300
    set seed 42
    gen cluster_id = cond(_n <= 100, 1, 2)
    gen x1 = rnormal()
    gen x2 = rnormal()
    gen y = 2*x1 + x2 + rnormal()

    grf_regression_forest y x1 x2, gen(pred_eq) ntrees(100) seed(42) ///
        cluster(cluster_id) equalizeclusterweights
    assert !missing(pred_eq) in 1
    drop pred_eq
}
if _rc {
    display as error "FAIL: E4.1 equalizeclusterweights regression"
    local errors = `errors' + 1
}
else {
    display as result "PASS: E4.1 equalizeclusterweights regression"
}

* Summary
display as text ""
if `errors' > 0 {
    display as error "FAILED: `errors' error(s)"
    exit 1
}
else {
    display as result "ALL CLUSTER TESTS PASSED"
}
