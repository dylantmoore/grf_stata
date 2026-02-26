* test_options_boosted.do -- Comprehensive option tests for grf_boosted_regression_forest
* Tests every option individually and in combination

clear all
set more off

local errors = 0

* ============================================================
* Setup data
* ============================================================
clear
set obs 500
set seed 42
forvalues j = 1/5 {
    gen x`j' = rnormal()
}
gen y = 2*x1 + x2^2 + rnormal()

* ---- Test 1: mtry(3) ----
capture noisily {
    grf_boosted_regression_forest y x1-x5, gen(bp1) ntrees(100) seed(42) mtry(3)
    assert !missing(bp1) in 1
    assert e(mtry) == 3
    drop bp1
}
if _rc {
    display as error "FAIL: mtry(3)"
    local errors = `errors' + 1
}
else {
    display as result "PASS: mtry(3)"
}

* ---- Test 2: minnodesize(10) ----
capture noisily {
    grf_boosted_regression_forest y x1-x5, gen(bp2) ntrees(100) seed(42) minnodesize(10)
    assert !missing(bp2) in 1
    assert e(min_node) == 10
    drop bp2
}
if _rc {
    display as error "FAIL: minnodesize(10)"
    local errors = `errors' + 1
}
else {
    display as result "PASS: minnodesize(10)"
}

* ---- Test 3: samplefrac(0.7) ----
* NOTE: boosted default is auto-tune (booststeps=0) which uses ci_group_size>=2,
*       requiring samplefrac < 0.5. Use booststeps(3) to bypass auto-tune.
capture noisily {
    grf_boosted_regression_forest y x1-x5, gen(bp3) ntrees(100) seed(42) samplefrac(0.7) booststeps(3)
    assert !missing(bp3) in 1
    assert reldif(e(sample_fraction), 0.7) < 1e-6
    drop bp3
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
    grf_boosted_regression_forest y x1-x5, gen(bp4) ntrees(100) seed(42) nohonesty
    assert !missing(bp4) in 1
    assert e(honesty) == 0
    drop bp4
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
    grf_boosted_regression_forest y x1-x5, gen(bp5) ntrees(100) seed(42) honestyfrac(0.7)
    assert !missing(bp5) in 1
    assert reldif(e(honesty_fraction), 0.7) < 1e-6
    drop bp5
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
    grf_boosted_regression_forest y x1-x5, gen(bp6) ntrees(100) seed(42) nohonestyprune
    assert !missing(bp6) in 1
    assert e(honesty_prune) == 0
    drop bp6
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
    grf_boosted_regression_forest y x1-x5, gen(bp7) ntrees(100) seed(42) alpha(0.1)
    assert !missing(bp7) in 1
    assert reldif(e(alpha), 0.1) < 1e-6
    drop bp7
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
    grf_boosted_regression_forest y x1-x5, gen(bp8) ntrees(100) seed(42) imbalancepenalty(0.5)
    assert !missing(bp8) in 1
    assert reldif(e(imbalance_penalty), 0.5) < 1e-6
    drop bp8
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
    grf_boosted_regression_forest y x1-x5, gen(bp9) ntrees(100) seed(42) cigroupsize(2) samplefrac(0.4)
    assert !missing(bp9) in 1
    assert e(ci_group_size) == 2
    drop bp9
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
    grf_boosted_regression_forest y x1-x5, gen(bp10) ntrees(100) seed(42) numthreads(2)
    assert !missing(bp10) in 1
    drop bp10
}
if _rc {
    display as error "FAIL: numthreads(2)"
    local errors = `errors' + 1
}
else {
    display as result "PASS: numthreads(2)"
}

* ---- Test 11: nostabilizesplits ----
capture noisily {
    grf_boosted_regression_forest y x1-x5, gen(bp11) ntrees(100) seed(42) nostabilizesplits
    assert !missing(bp11) in 1
    assert e(stabilize) == 0
    drop bp11
}
if _rc {
    display as error "FAIL: nostabilizesplits"
    local errors = `errors' + 1
}
else {
    display as result "PASS: nostabilizesplits"
}

* ---- Test 12: estimatevariance ----
* NOTE: Boosted regression is iterative -- variance estimation across boosting
*       steps is not meaningful. The option is accepted but _var will be missing.
capture noisily {
    grf_boosted_regression_forest y x1-x5, gen(bp12) ntrees(100) seed(42) estimatevariance booststeps(2)
    assert !missing(bp12) in 1
    * Variance variable should be created (even if all missing for boosted)
    confirm variable bp12_var
    drop bp12 bp12_var
}
if _rc {
    display as error "FAIL: estimatevariance"
    local errors = `errors' + 1
}
else {
    display as result "PASS: estimatevariance"
}

* ---- Test 13: estimatevariance + vargenerate ----
* NOTE: Same as Test 12 — variance is all missing for boosted, but vargenerate works.
capture noisily {
    grf_boosted_regression_forest y x1-x5, gen(bp13) ntrees(100) seed(42) ///
        estimatevariance vargenerate(bvar13) booststeps(2)
    assert !missing(bp13) in 1
    confirm variable bvar13
    assert "`e(variance_var)'" == "bvar13"
    drop bp13 bvar13
}
if _rc {
    display as error "FAIL: estimatevariance + vargenerate"
    local errors = `errors' + 1
}
else {
    display as result "PASS: estimatevariance + vargenerate"
}

* ---- Test 14: booststeps(3) explicit ----
capture noisily {
    grf_boosted_regression_forest y x1-x5, gen(bp14) ntrees(100) seed(42) booststeps(3)
    assert !missing(bp14) in 1
    assert e(boost_steps) == 3
    drop bp14
}
if _rc {
    display as error "FAIL: booststeps(3)"
    local errors = `errors' + 1
}
else {
    display as result "PASS: booststeps(3)"
}

* ---- Test 15: boostmaxsteps(3) ----
capture noisily {
    grf_boosted_regression_forest y x1-x5, gen(bp15) ntrees(100) seed(42) boostmaxsteps(3)
    assert !missing(bp15) in 1
    assert e(boost_max_steps) == 3
    assert e(boost_steps) <= 3
    drop bp15
}
if _rc {
    display as error "FAIL: boostmaxsteps(3)"
    local errors = `errors' + 1
}
else {
    display as result "PASS: boostmaxsteps(3)"
}

* ---- Test 16: boosterrorreduction(0.90) ----
capture noisily {
    grf_boosted_regression_forest y x1-x5, gen(bp16) ntrees(100) seed(42) boosterrorreduction(0.90)
    assert !missing(bp16) in 1
    assert reldif(e(boost_error_reduction), 0.90) < 1e-6
    drop bp16
}
if _rc {
    display as error "FAIL: boosterrorreduction(0.90)"
    local errors = `errors' + 1
}
else {
    display as result "PASS: boosterrorreduction(0.90)"
}

* ---- Test 17: boosttreestune(5) ----
capture noisily {
    grf_boosted_regression_forest y x1-x5, gen(bp17) ntrees(100) seed(42) boosttreestune(5)
    assert !missing(bp17) in 1
    drop bp17
}
if _rc {
    display as error "FAIL: boosttreestune(5)"
    local errors = `errors' + 1
}
else {
    display as result "PASS: boosttreestune(5)"
}

* ---- Test 18: if restriction ----
capture noisily {
    grf_boosted_regression_forest y x1-x5 if x1 > 0, gen(bp18) ntrees(100) seed(42)
    quietly count if x1 > 0
    local n_sub = r(N)
    assert e(N) == `n_sub'
    drop bp18
}
if _rc {
    display as error "FAIL: if restriction"
    local errors = `errors' + 1
}
else {
    display as result "PASS: if restriction"
}

* ---- Test 19: replace ----
capture noisily {
    grf_boosted_regression_forest y x1-x5, gen(bp19) ntrees(100) seed(42)
    grf_boosted_regression_forest y x1-x5, gen(bp19) ntrees(100) seed(42) replace
    assert !missing(bp19) in 1
    drop bp19
}
if _rc {
    display as error "FAIL: replace"
    local errors = `errors' + 1
}
else {
    display as result "PASS: replace"
}

* ---- Test 20: All non-default options combined ----
capture noisily {
    grf_boosted_regression_forest y x1-x5, gen(bp20) ntrees(100) seed(123) ///
        mtry(3) minnodesize(10) samplefrac(0.4) honestyfrac(0.7)          ///
        nohonestyprune alpha(0.1) imbalancepenalty(0.5) cigroupsize(2)     ///
        numthreads(2) nostabilizesplits estimatevariance                   ///
        vargenerate(bvar20) booststeps(2)
    assert !missing(bp20) in 1
    confirm variable bvar20
    assert e(mtry) == 3
    assert e(min_node) == 10
    assert reldif(e(sample_fraction), 0.4) < 1e-6
    assert reldif(e(honesty_fraction), 0.7) < 1e-6
    assert e(honesty_prune) == 0
    assert reldif(e(alpha), 0.1) < 1e-6
    assert reldif(e(imbalance_penalty), 0.5) < 1e-6
    assert e(ci_group_size) == 2
    assert e(stabilize) == 0
    assert e(boost_steps) == 2
    drop bp20 bvar20
}
if _rc {
    display as error "FAIL: all non-default options combined"
    local errors = `errors' + 1
}
else {
    display as result "PASS: all non-default options combined"
}

* ---- Test 21: Verify default e() values ----
capture noisily {
    grf_boosted_regression_forest y x1-x5, gen(bp21) ntrees(100) seed(42)
    assert e(honesty) == 1
    assert e(honesty_prune) == 1
    assert e(stabilize) == 1
    assert e(mtry) == 0
    assert e(min_node) == 5
    assert reldif(e(sample_fraction), 0.5) < 1e-6
    assert reldif(e(alpha), 0.05) < 1e-6
    assert e(boost_max_steps) == 5
    assert reldif(e(boost_error_reduction), 0.97) < 1e-6
    assert !missing(e(boost_steps))
    assert e(boost_steps) >= 0
    assert "`e(forest_type)'" == "boosted_regression"
    drop bp21
}
if _rc {
    display as error "FAIL: default e() values"
    local errors = `errors' + 1
}
else {
    display as result "PASS: default e() values"
}

* ---- Test 22: booststeps(0) auto-tune mode ----
capture noisily {
    grf_boosted_regression_forest y x1-x5, gen(bp22) ntrees(100) seed(42) booststeps(0)
    assert !missing(bp22) in 1
    * booststeps(0) = auto-tune, e(boost_steps) is the actual steps used
    assert e(boost_steps) >= 0
    drop bp22
}
if _rc {
    display as error "FAIL: booststeps(0) auto-tune"
    local errors = `errors' + 1
}
else {
    display as result "PASS: booststeps(0) auto-tune"
}

* ============================================================
* Summary
* ============================================================
display ""
if `errors' > 0 {
    display as error "FAILED: `errors' errors in boosted regression forest option tests"
    exit 1
}
else {
    display as result "ALL BOOSTED REGRESSION FOREST OPTION TESTS PASSED"
}
