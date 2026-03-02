* test_options_causal_survival.do -- Comprehensive option tests for grf_causal_survival_forest

clear all
set more off

local errors = 0

* ============================================================
* Setup causal survival data
* Syntax: grf_causal_survival_forest time status treat X1..Xp, gen(...)
* ============================================================
clear
set obs 500
set seed 42
forvalues j = 1/5 {
    gen x`j' = rnormal()
}
gen time = exp(x1 + 0.5*x2 + rnormal())
gen status = (uniform() > 0.3)
gen treat = (x1 + rnormal() > 0)

* ---- Test 1: mtry(3) ----
capture noisily {
    grf_causal_survival_forest time status treat x1-x5, gen(cs1) ntrees(100) seed(42) mtry(3)
    assert !missing(cs1) in 1
    assert e(mtry) == 3
    drop cs1 _grf_cs_what _grf_cs_yhat _grf_cs_shat _grf_cs_chat _grf_cs_numer _grf_cs_denom
}
if _rc {
    display as error "FAIL: mtry(3)"
    local errors = `errors' + 1
}
else {
    display as result "PASS: mtry(3)"
}

* ---- Test 2: minnodesize(20) ----
capture noisily {
    grf_causal_survival_forest time status treat x1-x5, gen(cs2) ntrees(100) seed(42) minnodesize(20)
    assert !missing(cs2) in 1
    assert e(min_node) == 20
    drop cs2 _grf_cs_what _grf_cs_yhat _grf_cs_shat _grf_cs_chat _grf_cs_numer _grf_cs_denom
}
if _rc {
    display as error "FAIL: minnodesize(20)"
    local errors = `errors' + 1
}
else {
    display as result "PASS: minnodesize(20)"
}

* ---- Test 3: samplefrac(0.7) ----
capture noisily {
    grf_causal_survival_forest time status treat x1-x5, gen(cs3) ntrees(100) seed(42) samplefrac(0.7)
    assert !missing(cs3) in 1
    assert reldif(e(sample_fraction), 0.7) < 1e-6
    drop cs3 _grf_cs_what _grf_cs_yhat _grf_cs_shat _grf_cs_chat _grf_cs_numer _grf_cs_denom
}
if _rc {
    display as error "FAIL: samplefrac(0.7)"
    local errors = `errors' + 1
}
else {
    display as result "PASS: samplefrac(0.7)"
}

* ---- Test 4: nohonesty ----
capture noisily {
    grf_causal_survival_forest time status treat x1-x5, gen(cs4) ntrees(100) seed(42) nohonesty
    assert !missing(cs4) in 1
    assert e(honesty) == 0
    drop cs4 _grf_cs_what _grf_cs_yhat _grf_cs_shat _grf_cs_chat _grf_cs_numer _grf_cs_denom
}
if _rc {
    display as error "FAIL: nohonesty"
    local errors = `errors' + 1
}
else {
    display as result "PASS: nohonesty"
}

* ---- Test 5: honestyfrac(0.7) ----
capture noisily {
    grf_causal_survival_forest time status treat x1-x5, gen(cs5) ntrees(100) seed(42) honestyfrac(0.7)
    assert !missing(cs5) in 1
    assert reldif(e(honesty_fraction), 0.7) < 1e-6
    drop cs5 _grf_cs_what _grf_cs_yhat _grf_cs_shat _grf_cs_chat _grf_cs_numer _grf_cs_denom
}
if _rc {
    display as error "FAIL: honestyfrac(0.7)"
    local errors = `errors' + 1
}
else {
    display as result "PASS: honestyfrac(0.7)"
}

* ---- Test 6: nohonestyprune ----
capture noisily {
    grf_causal_survival_forest time status treat x1-x5, gen(cs6) ntrees(100) seed(42) nohonestyprune
    assert !missing(cs6) in 1
    assert e(honesty_prune) == 0
    drop cs6 _grf_cs_what _grf_cs_yhat _grf_cs_shat _grf_cs_chat _grf_cs_numer _grf_cs_denom
}
if _rc {
    display as error "FAIL: nohonestyprune"
    local errors = `errors' + 1
}
else {
    display as result "PASS: nohonestyprune"
}

* ---- Test 7: alpha(0.1) ----
capture noisily {
    grf_causal_survival_forest time status treat x1-x5, gen(cs7) ntrees(100) seed(42) alpha(0.1)
    assert !missing(cs7) in 1
    assert reldif(e(alpha), 0.1) < 1e-6
    drop cs7 _grf_cs_what _grf_cs_yhat _grf_cs_shat _grf_cs_chat _grf_cs_numer _grf_cs_denom
}
if _rc {
    display as error "FAIL: alpha(0.1)"
    local errors = `errors' + 1
}
else {
    display as result "PASS: alpha(0.1)"
}

* ---- Test 8: imbalancepenalty(0.5) ----
capture noisily {
    grf_causal_survival_forest time status treat x1-x5, gen(cs8) ntrees(100) seed(42) imbalancepenalty(0.5)
    assert !missing(cs8) in 1
    assert reldif(e(imbalance_penalty), 0.5) < 1e-6
    drop cs8 _grf_cs_what _grf_cs_yhat _grf_cs_shat _grf_cs_chat _grf_cs_numer _grf_cs_denom
}
if _rc {
    display as error "FAIL: imbalancepenalty(0.5)"
    local errors = `errors' + 1
}
else {
    display as result "PASS: imbalancepenalty(0.5)"
}

* ---- Test 9: cigroupsize(2) ----
capture noisily {
    grf_causal_survival_forest time status treat x1-x5, gen(cs9) ntrees(100) seed(42) cigroupsize(2) samplefrac(0.4)
    assert !missing(cs9) in 1
    assert e(ci_group_size) == 2
    drop cs9 _grf_cs_what _grf_cs_yhat _grf_cs_shat _grf_cs_chat _grf_cs_numer _grf_cs_denom
}
if _rc {
    display as error "FAIL: cigroupsize(2)"
    local errors = `errors' + 1
}
else {
    display as result "PASS: cigroupsize(2)"
}

* ---- Test 10: numthreads(2) ----
capture noisily {
    grf_causal_survival_forest time status treat x1-x5, gen(cs10) ntrees(100) seed(42) numthreads(2)
    assert !missing(cs10) in 1
    drop cs10 _grf_cs_what _grf_cs_yhat _grf_cs_shat _grf_cs_chat _grf_cs_numer _grf_cs_denom
}
if _rc {
    display as error "FAIL: numthreads(2)"
    local errors = `errors' + 1
}
else {
    display as result "PASS: numthreads(2)"
}

* ---- Test 11: target(1) RMST (default) ----
capture noisily {
    grf_causal_survival_forest time status treat x1-x5, gen(cs11) ntrees(100) seed(42) target(1)
    assert !missing(cs11) in 1
    assert e(target) == 1
    drop cs11 _grf_cs_what _grf_cs_yhat _grf_cs_shat _grf_cs_chat _grf_cs_numer _grf_cs_denom
}
if _rc {
    display as error "FAIL: target(1) RMST"
    local errors = `errors' + 1
}
else {
    display as result "PASS: target(1) RMST"
}

* ---- Test 12: target(2) survival probability ----
capture noisily {
    grf_causal_survival_forest time status treat x1-x5, gen(cs12) ntrees(100) seed(42) target(2)
    assert !missing(cs12) in 1
    assert e(target) == 2
    drop cs12 _grf_cs_what _grf_cs_yhat _grf_cs_shat _grf_cs_chat _grf_cs_numer _grf_cs_denom
}
if _rc {
    display as error "FAIL: target(2) survival probability"
    local errors = `errors' + 1
}
else {
    display as result "PASS: target(2) survival probability"
}

* ---- Test 13: estimatevariance ----
capture noisily {
    grf_causal_survival_forest time status treat x1-x5, gen(cs13) ntrees(100) seed(42) estimatevariance
    assert !missing(cs13) in 1
    assert !missing(cs13_var) in 1
    assert cs13_var[1] > 0
    assert "`e(variance_var)'" == "cs13_var"
    drop cs13 cs13_var _grf_cs_what _grf_cs_yhat _grf_cs_shat _grf_cs_chat _grf_cs_numer _grf_cs_denom
}
if _rc {
    display as error "FAIL: estimatevariance"
    local errors = `errors' + 1
}
else {
    display as result "PASS: estimatevariance"
}

* ---- Test 14: estimatevariance + vargenerate ----
capture noisily {
    grf_causal_survival_forest time status treat x1-x5, gen(cs14) ntrees(100) seed(42) ///
        estimatevariance vargenerate(csvar14)
    assert !missing(cs14) in 1
    assert !missing(csvar14) in 1
    assert "`e(variance_var)'" == "csvar14"
    drop cs14 csvar14 _grf_cs_what _grf_cs_yhat _grf_cs_shat _grf_cs_chat _grf_cs_numer _grf_cs_denom
}
if _rc {
    display as error "FAIL: estimatevariance + vargenerate"
    local errors = `errors' + 1
}
else {
    display as result "PASS: estimatevariance + vargenerate"
}

* ---- Test 15: nostabilizesplits ----
capture noisily {
    grf_causal_survival_forest time status treat x1-x5, gen(cs15) ntrees(100) seed(42) nostabilizesplits
    assert !missing(cs15) in 1
    assert e(stabilize) == 0
    drop cs15 _grf_cs_what _grf_cs_yhat _grf_cs_shat _grf_cs_chat _grf_cs_numer _grf_cs_denom
}
if _rc {
    display as error "FAIL: nostabilizesplits"
    local errors = `errors' + 1
}
else {
    display as result "PASS: nostabilizesplits"
}

* ---- Test 16: horizon(5) explicit ----
capture noisily {
    grf_causal_survival_forest time status treat x1-x5, gen(cs16) ntrees(100) seed(42) horizon(5)
    assert !missing(cs16) in 1
    assert reldif(e(horizon), 5) < 1e-6
    drop cs16 _grf_cs_what _grf_cs_yhat _grf_cs_shat _grf_cs_chat _grf_cs_numer _grf_cs_denom
}
if _rc {
    display as error "FAIL: horizon(5)"
    local errors = `errors' + 1
}
else {
    display as result "PASS: horizon(5)"
}

* ---- Test 17: if restriction ----
capture noisily {
    grf_causal_survival_forest time status treat x1-x5 if x1 > 0, gen(cs17) ntrees(100) seed(42)
    quietly count if x1 > 0
    local n_sub = r(N)
    assert e(N) == `n_sub'
    drop cs17 _grf_cs_what _grf_cs_yhat _grf_cs_shat _grf_cs_chat _grf_cs_numer _grf_cs_denom
}
if _rc {
    display as error "FAIL: if restriction"
    local errors = `errors' + 1
}
else {
    display as result "PASS: if restriction"
}

* ---- Test 18: replace ----
capture noisily {
    grf_causal_survival_forest time status treat x1-x5, gen(cs18) ntrees(100) seed(42)
    grf_causal_survival_forest time status treat x1-x5, gen(cs18) ntrees(100) seed(42) replace
    assert !missing(cs18) in 1
    drop cs18 _grf_cs_what _grf_cs_yhat _grf_cs_shat _grf_cs_chat _grf_cs_numer _grf_cs_denom
}
if _rc {
    display as error "FAIL: replace"
    local errors = `errors' + 1
}
else {
    display as result "PASS: replace"
}

* ---- Test 19: all major non-default options combined ----
capture noisily {
    grf_causal_survival_forest time status treat x1-x5, gen(cs19) ntrees(100) seed(123) ///
        mtry(3) minnodesize(20) samplefrac(0.4) honestyfrac(0.7)                        ///
        nohonestyprune alpha(0.1) imbalancepenalty(0.5) cigroupsize(2)                   ///
        numthreads(2) nostabilizesplits target(2) horizon(5)                             ///
        estimatevariance vargenerate(csvar19)
    assert !missing(cs19) in 1
    assert !missing(csvar19) in 1
    assert e(mtry) == 3
    assert e(min_node) == 20
    assert reldif(e(sample_fraction), 0.4) < 1e-6
    assert reldif(e(honesty_fraction), 0.7) < 1e-6
    assert e(honesty_prune) == 0
    assert reldif(e(alpha), 0.1) < 1e-6
    assert reldif(e(imbalance_penalty), 0.5) < 1e-6
    assert e(ci_group_size) == 2
    assert e(stabilize) == 0
    assert e(target) == 2
    assert reldif(e(horizon), 5) < 1e-6
    drop cs19 csvar19 _grf_cs_what _grf_cs_yhat _grf_cs_shat _grf_cs_chat _grf_cs_numer _grf_cs_denom
}
if _rc {
    display as error "FAIL: all non-default options combined"
    local errors = `errors' + 1
}
else {
    display as result "PASS: all non-default options combined"
}

* ---- Test 20: default e() values and nuisance pointers ----
capture noisily {
    grf_causal_survival_forest time status treat x1-x5, gen(cs20) ntrees(100) seed(42)
    assert e(honesty) == 1
    assert e(honesty_prune) == 1
    assert e(stabilize) == 1
    assert e(target) == 1
    assert e(mtry) == 0
    assert e(min_node) == 15
    assert reldif(e(sample_fraction), 0.5) < 1e-6
    assert reldif(e(alpha), 0.05) < 1e-6
    assert !missing(e(horizon))
    assert e(horizon) > 0
    assert !missing(e(ate))
    assert !missing(e(ate_se))
    assert "`e(forest_type)'" == "causal_survival"
    assert "`e(timevar)'" == "time"
    assert "`e(statusvar)'" == "status"
    assert "`e(treatvar)'" == "treat"
    assert "`e(what_var)'" == "_grf_cs_what"
    assert "`e(yhat_var)'" == "_grf_cs_yhat"
    assert "`e(shat_var)'" == "_grf_cs_shat"
    assert "`e(chat_var)'" == "_grf_cs_chat"
    assert "`e(numer_var)'" == "_grf_cs_numer"
    assert "`e(denom_var)'" == "_grf_cs_denom"
    assert "`e(nuisance_mode)'" == "auto"
    drop cs20 _grf_cs_what _grf_cs_yhat _grf_cs_shat _grf_cs_chat _grf_cs_numer _grf_cs_denom
}
if _rc {
    display as error "FAIL: default e() values"
    local errors = `errors' + 1
}
else {
    display as result "PASS: default e() values"
}

* ---- Test 21: full-input nuisance mode requires full set and runs ----
capture noisily {
    gen double what_in = min(max(0.5 + 0.2*x1, 0.01), 0.99)
    gen double yhat_in = min(time, 4)
    gen double shat_in = min(max((time > 4) * 0.9 + 0.05, 0.01), 0.99)
    gen double chat_in = min(max(0.8 + 0.1*x2, 0.05), 0.99)

    grf_causal_survival_forest time status treat x1-x5, gen(cs21) ntrees(100) seed(42) ///
        horizon(4) whatinput(what_in) yhatinput(yhat_in) shatinput(shat_in) chatinput(chat_in)

    assert !missing(cs21) in 1
    assert "`e(nuisance_mode)'" == "full_input"
    assert !missing(_grf_cs_what) in 1
    assert !missing(_grf_cs_yhat) in 1
    assert !missing(_grf_cs_shat) in 1
    assert !missing(_grf_cs_chat) in 1

    drop cs21 what_in yhat_in shat_in chat_in ///
        _grf_cs_what _grf_cs_yhat _grf_cs_shat _grf_cs_chat _grf_cs_numer _grf_cs_denom
}
if _rc {
    display as error "FAIL: full-input nuisance mode"
    local errors = `errors' + 1
}
else {
    display as result "PASS: full-input nuisance mode"
}

* ---- Test 22: moment-input numer()/denom() mode ----
capture noisily {
    gen double numer_pre = (treat - 0.5) * status * min(time, 3)
    gen double denom_pre = (treat - 0.5)^2 + 1e-4

    grf_causal_survival_forest time status treat x1-x5, gen(cs22) ntrees(100) seed(42) ///
        horizon(3) numer(numer_pre) denom(denom_pre)

    assert !missing(cs22) in 1
    assert "`e(nuisance_mode)'" == "moment_input"
    assert !missing(_grf_cs_numer) in 1
    assert !missing(_grf_cs_denom) in 1

    drop cs22 numer_pre denom_pre ///
        _grf_cs_what _grf_cs_yhat _grf_cs_shat _grf_cs_chat _grf_cs_numer _grf_cs_denom
}
if _rc {
    display as error "FAIL: moment-input numer()/denom() mode"
    local errors = `errors' + 1
}
else {
    display as result "PASS: moment-input numer()/denom() mode"
}

* ---- Test 23: partial nuisance input should fail ----
capture noisily {
    gen double what_only = min(max(0.5 + 0.2*x1, 0.01), 0.99)
    capture grf_causal_survival_forest time status treat x1-x5, gen(cs23) ntrees(100) seed(42) whatinput(what_only)
    assert _rc != 0
    drop what_only
}
if _rc {
    display as error "FAIL: partial nuisance input rejection"
    local errors = `errors' + 1
}
else {
    display as result "PASS: partial nuisance input rejection"
}

* ---- Test 24: numer()/denom() cannot combine with nuisance input ----
capture noisily {
    gen double numer_conflict = (treat - 0.5) * status
    gen double denom_conflict = (treat - 0.5)^2 + 1e-4
    gen double what_conflict = min(max(0.5 + 0.1*x1, 0.01), 0.99)
    capture grf_causal_survival_forest time status treat x1-x5, gen(cs24) ntrees(100) seed(42) ///
        numer(numer_conflict) denom(denom_conflict) whatinput(what_conflict)
    assert _rc != 0
    drop numer_conflict denom_conflict what_conflict
}
if _rc {
    display as error "FAIL: numer()/denom() conflict rejection"
    local errors = `errors' + 1
}
else {
    display as result "PASS: numer()/denom() conflict rejection"
}

* ---- Test 25: nuisance generate outputs ----
capture noisily {
    grf_causal_survival_forest time status treat x1-x5, gen(cs25) ntrees(100) seed(42) ///
        whatgenerate(what_out) yhatgenerate(yhat_out) shatgenerate(shat_out) chatgenerate(chat_out)

    assert !missing(cs25) in 1
    assert !missing(what_out) in 1
    assert !missing(yhat_out) in 1
    assert !missing(shat_out) in 1
    assert !missing(chat_out) in 1
    assert "`e(what_generate)'" == "what_out"
    assert "`e(yhat_generate)'" == "yhat_out"
    assert "`e(shat_generate)'" == "shat_out"
    assert "`e(chat_generate)'" == "chat_out"

    drop cs25 what_out yhat_out shat_out chat_out ///
        _grf_cs_what _grf_cs_yhat _grf_cs_shat _grf_cs_chat _grf_cs_numer _grf_cs_denom
}
if _rc {
    display as error "FAIL: nuisance generate outputs"
    local errors = `errors' + 1
}
else {
    display as result "PASS: nuisance generate outputs"
}

* ============================================================
* Summary
* ============================================================
display ""
if `errors' > 0 {
    display as error "FAILED: `errors' errors in causal survival forest option tests"
    exit 1
}
else {
    display as result "ALL CAUSAL SURVIVAL FOREST OPTION TESTS PASSED"
}
