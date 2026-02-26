* test_options_multi_arm.do -- Comprehensive option tests for grf_multi_arm_causal_forest
* Tests every option individually and in combination

clear all
set more off

local errors = 0

* ============================================================
* Setup data
* Syntax: grf_multi_arm_causal_forest depvar treat1 treat2 ... X1..Xp, gen(...) ntreat(K)
* ============================================================
clear
set obs 500
set seed 42
forvalues j = 1/5 {
    gen x`j' = rnormal()
}
gen y = 2*x1 + x2^2 + rnormal()
gen t1 = (x1 + rnormal() > 0)
gen t2 = (x2 + rnormal() > 0.5)

* ---- Test 1: mtry(3) ----
capture noisily {
    grf_multi_arm_causal_forest y t1 t2 x1-x5, gen(mac1) ntreat(2) ntrees(100) seed(42) mtry(3)
    assert !missing(mac1_t1) in 1
    assert !missing(mac1_t2) in 1
    assert e(mtry) == 3
    drop mac1_t1 mac1_t2 _grf_mac_yhat _grf_mac_what1 _grf_mac_what2
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
    grf_multi_arm_causal_forest y t1 t2 x1-x5, gen(mac2) ntreat(2) ntrees(100) seed(42) minnodesize(10)
    assert !missing(mac2_t1) in 1
    assert e(min_node) == 10
    drop mac2_t1 mac2_t2 _grf_mac_yhat _grf_mac_what1 _grf_mac_what2
}
if _rc {
    display as error "FAIL: minnodesize(10)"
    local errors = `errors' + 1
}
else {
    display as result "PASS: minnodesize(10)"
}

* ---- Test 3: samplefrac(0.7) ----
capture noisily {
    grf_multi_arm_causal_forest y t1 t2 x1-x5, gen(mac3) ntreat(2) ntrees(100) seed(42) samplefrac(0.7)
    assert !missing(mac3_t1) in 1
    assert reldif(e(sample_fraction), 0.7) < 1e-6
    drop mac3_t1 mac3_t2 _grf_mac_yhat _grf_mac_what1 _grf_mac_what2
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
    grf_multi_arm_causal_forest y t1 t2 x1-x5, gen(mac4) ntreat(2) ntrees(100) seed(42) nohonesty
    assert !missing(mac4_t1) in 1
    assert e(honesty) == 0
    drop mac4_t1 mac4_t2 _grf_mac_yhat _grf_mac_what1 _grf_mac_what2
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
    grf_multi_arm_causal_forest y t1 t2 x1-x5, gen(mac5) ntreat(2) ntrees(100) seed(42) honestyfrac(0.7)
    assert !missing(mac5_t1) in 1
    assert reldif(e(honesty_fraction), 0.7) < 1e-6
    drop mac5_t1 mac5_t2 _grf_mac_yhat _grf_mac_what1 _grf_mac_what2
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
    grf_multi_arm_causal_forest y t1 t2 x1-x5, gen(mac6) ntreat(2) ntrees(100) seed(42) nohonestyprune
    assert !missing(mac6_t1) in 1
    assert e(honesty_prune) == 0
    drop mac6_t1 mac6_t2 _grf_mac_yhat _grf_mac_what1 _grf_mac_what2
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
    grf_multi_arm_causal_forest y t1 t2 x1-x5, gen(mac7) ntreat(2) ntrees(100) seed(42) alpha(0.1)
    assert !missing(mac7_t1) in 1
    assert reldif(e(alpha), 0.1) < 1e-6
    drop mac7_t1 mac7_t2 _grf_mac_yhat _grf_mac_what1 _grf_mac_what2
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
    grf_multi_arm_causal_forest y t1 t2 x1-x5, gen(mac8) ntreat(2) ntrees(100) seed(42) imbalancepenalty(0.5)
    assert !missing(mac8_t1) in 1
    assert reldif(e(imbalance_penalty), 0.5) < 1e-6
    drop mac8_t1 mac8_t2 _grf_mac_yhat _grf_mac_what1 _grf_mac_what2
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
    grf_multi_arm_causal_forest y t1 t2 x1-x5, gen(mac9) ntreat(2) ntrees(100) seed(42) cigroupsize(2) samplefrac(0.4)
    assert !missing(mac9_t1) in 1
    assert e(ci_group_size) == 2
    drop mac9_t1 mac9_t2 _grf_mac_yhat _grf_mac_what1 _grf_mac_what2
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
    grf_multi_arm_causal_forest y t1 t2 x1-x5, gen(mac10) ntreat(2) ntrees(100) seed(42) numthreads(2)
    assert !missing(mac10_t1) in 1
    drop mac10_t1 mac10_t2 _grf_mac_yhat _grf_mac_what1 _grf_mac_what2
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
    grf_multi_arm_causal_forest y t1 t2 x1-x5, gen(mac11) ntreat(2) ntrees(100) seed(42) nostabilizesplits
    assert !missing(mac11_t1) in 1
    assert e(stabilize) == 0
    drop mac11_t1 mac11_t2 _grf_mac_yhat _grf_mac_what1 _grf_mac_what2
}
if _rc {
    display as error "FAIL: nostabilizesplits"
    local errors = `errors' + 1
}
else {
    display as result "PASS: nostabilizesplits"
}

* ---- Test 12: estimatevariance ----
capture noisily {
    grf_multi_arm_causal_forest y t1 t2 x1-x5, gen(mac12) ntreat(2) ntrees(100) seed(42) estimatevariance
    assert !missing(mac12_t1) in 1
    assert !missing(mac12_t2) in 1
    assert !missing(mac12_t1_var) in 1
    assert !missing(mac12_t2_var) in 1
    assert mac12_t1_var[1] > 0
    assert mac12_t2_var[1] > 0
    drop mac12_t1 mac12_t2 mac12_t1_var mac12_t2_var _grf_mac_yhat _grf_mac_what1 _grf_mac_what2
}
if _rc {
    display as error "FAIL: estimatevariance"
    local errors = `errors' + 1
}
else {
    display as result "PASS: estimatevariance"
}

* ---- Test 13: if restriction ----
capture noisily {
    grf_multi_arm_causal_forest y t1 t2 x1-x5 if x1 > 0, gen(mac13) ntreat(2) ntrees(100) seed(42)
    quietly count if x1 > 0
    local n_sub = r(N)
    assert e(N) == `n_sub'
    drop mac13_t1 mac13_t2 _grf_mac_yhat _grf_mac_what1 _grf_mac_what2
}
if _rc {
    display as error "FAIL: if restriction"
    local errors = `errors' + 1
}
else {
    display as result "PASS: if restriction"
}

* ---- Test 14: replace ----
capture noisily {
    grf_multi_arm_causal_forest y t1 t2 x1-x5, gen(mac14) ntreat(2) ntrees(100) seed(42)
    grf_multi_arm_causal_forest y t1 t2 x1-x5, gen(mac14) ntreat(2) ntrees(100) seed(42) replace
    assert !missing(mac14_t1) in 1
    drop mac14_t1 mac14_t2 _grf_mac_yhat _grf_mac_what1 _grf_mac_what2
}
if _rc {
    display as error "FAIL: replace"
    local errors = `errors' + 1
}
else {
    display as result "PASS: replace"
}

* ---- Test 15: All non-default options combined ----
capture noisily {
    grf_multi_arm_causal_forest y t1 t2 x1-x5, gen(mac15) ntreat(2) ntrees(100) seed(123) ///
        mtry(3) minnodesize(10) samplefrac(0.4) honestyfrac(0.7)                            ///
        nohonestyprune alpha(0.1) imbalancepenalty(0.5) cigroupsize(2)                       ///
        numthreads(2) nostabilizesplits estimatevariance
    assert !missing(mac15_t1) in 1
    assert !missing(mac15_t2) in 1
    assert !missing(mac15_t1_var) in 1
    assert !missing(mac15_t2_var) in 1
    assert e(mtry) == 3
    assert e(min_node) == 10
    assert reldif(e(sample_fraction), 0.4) < 1e-6
    assert reldif(e(honesty_fraction), 0.7) < 1e-6
    assert e(honesty_prune) == 0
    assert reldif(e(alpha), 0.1) < 1e-6
    assert reldif(e(imbalance_penalty), 0.5) < 1e-6
    assert e(ci_group_size) == 2
    assert e(stabilize) == 0
    assert e(n_treat) == 2
    drop mac15_t1 mac15_t2 mac15_t1_var mac15_t2_var _grf_mac_yhat _grf_mac_what1 _grf_mac_what2
}
if _rc {
    display as error "FAIL: all non-default options combined"
    local errors = `errors' + 1
}
else {
    display as result "PASS: all non-default options combined"
}

* ---- Test 16: Verify default e() values ----
capture noisily {
    grf_multi_arm_causal_forest y t1 t2 x1-x5, gen(mac16) ntreat(2) ntrees(100) seed(42)
    assert e(honesty) == 1
    assert e(honesty_prune) == 1
    assert e(stabilize) == 1
    assert e(mtry) == 0
    assert e(min_node) == 5
    assert reldif(e(sample_fraction), 0.5) < 1e-6
    assert reldif(e(alpha), 0.05) < 1e-6
    assert e(n_treat) == 2
    assert "`e(forest_type)'" == "multi_causal"
    assert "`e(depvar)'" == "y"
    assert !missing(e(ate_1))
    assert !missing(e(ate_se_1))
    assert !missing(e(ate_2))
    assert !missing(e(ate_se_2))
    drop mac16_t1 mac16_t2 _grf_mac_yhat _grf_mac_what1 _grf_mac_what2
}
if _rc {
    display as error "FAIL: default e() values"
    local errors = `errors' + 1
}
else {
    display as result "PASS: default e() values"
}

* ---- Test 17: Single treatment arm (ntreat(1)) ----
capture noisily {
    grf_multi_arm_causal_forest y t1 x1-x5, gen(mac17) ntreat(1) ntrees(100) seed(42)
    assert !missing(mac17_t1) in 1
    assert e(n_treat) == 1
    drop mac17_t1 _grf_mac_yhat _grf_mac_what1
}
if _rc {
    display as error "FAIL: ntreat(1)"
    local errors = `errors' + 1
}
else {
    display as result "PASS: ntreat(1)"
}

* ============================================================
* Summary
* ============================================================
display ""
if `errors' > 0 {
    display as error "FAILED: `errors' errors in multi-arm causal forest option tests"
    exit 1
}
else {
    display as result "ALL MULTI-ARM CAUSAL FOREST OPTION TESTS PASSED"
}
