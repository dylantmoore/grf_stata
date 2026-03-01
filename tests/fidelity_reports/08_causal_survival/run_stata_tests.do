* Stata fidelity tests for grf_causal_survival_forest (tests 1-16)
* Work dir: /tmp/grf_stata/tests/fidelity_reports/08_causal_survival/

adopath ++ /tmp/grf_stata
set more off
set linesize 120

local outdir "/tmp/grf_stata/tests/fidelity_reports/08_causal_survival"

* DGP constants from R
local horizon       = 0.532911631814285
local horizon_q75   = 1.14894156643613
local horizon_q25   = 0.21734178687488
local horizon_q90   = 2.00056696149032
local horizon_hc    = 0.350721944153199
local horizon_unbal = 0.585108221029108

* ---- Load data ----
import delimited "`outdir'/test_data.csv", clear varnames(1) case(lower)
describe

* ---- Test 1: Default RMST target, horizon=median(Y) ----
di _n "=== Test 1: Default RMST, horizon=median(Y) ==="
capture drop tau
grf_causal_survival_forest time status w x1 x2 x3 x4 x5, ///
    gen(tau) ntrees(500) seed(42) horizon(`horizon') target(1) replace
if _rc == 0 {
    di "Test 1: Stata OK"
    export delimited tau using "`outdir'/s_pred_01.csv", replace
}
else {
    di "Test 1: FAILED rc=" _rc
}

* ---- Test 2: Explicit horizon = Q75 ----
di _n "=== Test 2: Explicit horizon=Q75(Y) ==="
capture drop tau
grf_causal_survival_forest time status w x1 x2 x3 x4 x5, ///
    gen(tau) ntrees(500) seed(42) horizon(`horizon_q75') target(1) replace
if _rc == 0 {
    di "Test 2: Stata OK"
    export delimited tau using "`outdir'/s_pred_02.csv", replace
}
else {
    di "Test 2: FAILED rc=" _rc
}

* ---- Test 3: Survival probability target ----
di _n "=== Test 3: Survival probability target ==="
capture drop tau
grf_causal_survival_forest time status w x1 x2 x3 x4 x5, ///
    gen(tau) ntrees(500) seed(42) horizon(`horizon') target(2) replace
if _rc == 0 {
    di "Test 3: Stata OK"
    export delimited tau using "`outdir'/s_pred_03.csv", replace
}
else {
    di "Test 3: FAILED rc=" _rc
}

* ---- Test 4: nostabilizesplits ----
di _n "=== Test 4: nostabilizesplits ==="
capture drop tau
grf_causal_survival_forest time status w x1 x2 x3 x4 x5, ///
    gen(tau) ntrees(500) seed(42) horizon(`horizon') target(1) nostabilizesplits replace
if _rc == 0 {
    di "Test 4: Stata OK"
    export delimited tau using "`outdir'/s_pred_04.csv", replace
}
else {
    di "Test 4: FAILED rc=" _rc
}

* ---- Test 5: User-supplied W.hat = 0.5 ----
di _n "=== Test 5: User-supplied W.hat=0.5 ==="
capture drop tau what_fixed
gen double what_fixed = 0.5
grf_causal_survival_forest time status w x1 x2 x3 x4 x5, ///
    gen(tau) ntrees(500) seed(42) horizon(`horizon') target(1) whatinput(what_fixed) replace
if _rc == 0 {
    di "Test 5: Stata OK"
    export delimited tau using "`outdir'/s_pred_05.csv", replace
}
else {
    di "Test 5: FAILED rc=" _rc
}

* ---- Test 6: With cluster() ----
di _n "=== Test 6: With cluster() ==="
capture drop tau
grf_causal_survival_forest time status w x1 x2 x3 x4 x5, ///
    gen(tau) ntrees(500) seed(42) horizon(`horizon') target(1) cluster(cluster_var) replace
if _rc == 0 {
    di "Test 6: Stata OK"
    export delimited tau using "`outdir'/s_pred_06.csv", replace
}
else {
    di "Test 6: FAILED rc=" _rc
}

* ---- Test 7: With weights ----
di _n "=== Test 7: With observation weights ==="
capture drop tau
grf_causal_survival_forest time status w x1 x2 x3 x4 x5, ///
    gen(tau) ntrees(500) seed(42) horizon(`horizon') target(1) weights(weights_var) replace
if _rc == 0 {
    di "Test 7: Stata OK"
    export delimited tau using "`outdir'/s_pred_07.csv", replace
}
else {
    di "Test 7: FAILED rc=" _rc
}

* ---- Test 8: nohonesty ----
di _n "=== Test 8: nohonesty ==="
capture drop tau
grf_causal_survival_forest time status w x1 x2 x3 x4 x5, ///
    gen(tau) ntrees(500) seed(42) horizon(`horizon') target(1) nohonesty replace
if _rc == 0 {
    di "Test 8: Stata OK"
    export delimited tau using "`outdir'/s_pred_08.csv", replace
}
else {
    di "Test 8: FAILED rc=" _rc
}

* ---- Test 9: mtry=2 ----
di _n "=== Test 9: mtry=2 ==="
capture drop tau
grf_causal_survival_forest time status w x1 x2 x3 x4 x5, ///
    gen(tau) ntrees(500) seed(42) horizon(`horizon') target(1) mtry(2) replace
if _rc == 0 {
    di "Test 9: Stata OK"
    export delimited tau using "`outdir'/s_pred_09.csv", replace
}
else {
    di "Test 9: FAILED rc=" _rc
}

* ---- Test 10: minnodesize=20 ----
di _n "=== Test 10: minnodesize=20 ==="
capture drop tau
grf_causal_survival_forest time status w x1 x2 x3 x4 x5, ///
    gen(tau) ntrees(500) seed(42) horizon(`horizon') target(1) minnodesize(20) replace
if _rc == 0 {
    di "Test 10: Stata OK"
    export delimited tau using "`outdir'/s_pred_10.csv", replace
}
else {
    di "Test 10: FAILED rc=" _rc
}

* ---- Test 11: Heavy censoring (~40% censored) ----
di _n "=== Test 11: Heavy censoring ==="
capture drop tau
grf_causal_survival_forest time_hc status_hc w x1 x2 x3 x4 x5, ///
    gen(tau) ntrees(500) seed(42) horizon(`horizon_hc') target(1) replace
if _rc == 0 {
    di "Test 11: Stata OK"
    export delimited tau using "`outdir'/s_pred_11.csv", replace
}
else {
    di "Test 11: FAILED rc=" _rc
}

* ---- Test 12: Balanced treatment 50/50 (same as Test 1 DGP) ----
di _n "=== Test 12: Balanced treatment 50/50 ==="
capture drop tau
grf_causal_survival_forest time status w x1 x2 x3 x4 x5, ///
    gen(tau) ntrees(500) seed(42) horizon(`horizon') target(1) replace
if _rc == 0 {
    di "Test 12: Stata OK"
    export delimited tau using "`outdir'/s_pred_12.csv", replace
}
else {
    di "Test 12: FAILED rc=" _rc
}

* ---- Test 13: Unbalanced treatment 70/30 ----
di _n "=== Test 13: Unbalanced treatment 70/30 ==="
capture drop tau
grf_causal_survival_forest time_unbal status_unbal w_unbal x1 x2 x3 x4 x5, ///
    gen(tau) ntrees(500) seed(42) horizon(`horizon_unbal') target(1) replace
if _rc == 0 {
    di "Test 13: Stata OK"
    export delimited tau using "`outdir'/s_pred_13.csv", replace
}
else {
    di "Test 13: FAILED rc=" _rc
}

* ---- Test 14: Fewer trees (nuisancetrees=100 via ntrees=100) ----
di _n "=== Test 14: Fewer nuisance trees (ntrees=100) ==="
capture drop tau
grf_causal_survival_forest time status w x1 x2 x3 x4 x5, ///
    gen(tau) ntrees(100) seed(42) horizon(`horizon') target(1) replace
if _rc == 0 {
    di "Test 14: Stata OK"
    export delimited tau using "`outdir'/s_pred_14.csv", replace
}
else {
    di "Test 14: FAILED rc=" _rc
}

* ---- Test 15: Small horizon (Q25 of failure times) ----
di _n "=== Test 15: Small horizon=Q25(failure) ==="
capture drop tau
grf_causal_survival_forest time status w x1 x2 x3 x4 x5, ///
    gen(tau) ntrees(500) seed(42) horizon(`horizon_q25') target(1) replace
if _rc == 0 {
    di "Test 15: Stata OK"
    export delimited tau using "`outdir'/s_pred_15.csv", replace
}
else {
    di "Test 15: FAILED rc=" _rc
}

* ---- Test 16: Large horizon (Q90 of failure times) ----
di _n "=== Test 16: Large horizon=Q90(failure) ==="
capture drop tau
grf_causal_survival_forest time status w x1 x2 x3 x4 x5, ///
    gen(tau) ntrees(500) seed(42) horizon(`horizon_q90') target(1) replace
if _rc == 0 {
    di "Test 16: Stata OK"
    export delimited tau using "`outdir'/s_pred_16.csv", replace
}
else {
    di "Test 16: FAILED rc=" _rc
}

di _n "===== ALL STATA TESTS COMPLETE ====="
