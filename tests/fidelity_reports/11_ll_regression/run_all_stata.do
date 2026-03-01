* Master Stata do-file: Run all 18 ll_regression_forest fidelity tests
* Each test loads R-generated CSV, runs Stata, exports results for comparison
set more off
adopath ++ /tmp/grf_stata

local OUTDIR "/tmp/grf_stata/tests/fidelity_reports/11_ll_regression"

display ""
display "============================================================"
display " GRF Local Linear Regression Forest: Stata Fidelity Tests"
display "============================================================"

* ============================================================
* Test 01: Default LL (all vars, lambda=0.1)
* Stata: grf_ll_regression_forest y x1-x5, gen(pred) ntrees(500) seed(42)
* (no llenable by default — LL correction at predict time via llvars)
* ============================================================
display ""
display "--- Test 01: Default LL (all vars, lllambda=0.1) ---"
import delimited "`OUTDIR'/test01_data.csv", clear
grf_ll_regression_forest y x1 x2 x3 x4 x5, gen(stata_pred) ntrees(500) seed(42) lllambda(0.1) llvars(x1 x2 x3 x4 x5)
keep stata_pred r_pred
export delimited "`OUTDIR'/test01_stata.csv", replace
display "  Test 01 done."

* ============================================================
* Test 02: llenable (enable_ll_split=TRUE)
* ============================================================
display ""
display "--- Test 02: llenable ---"
import delimited "`OUTDIR'/test02_data.csv", clear
grf_ll_regression_forest y x1 x2 x3 x4 x5, gen(stata_pred) ntrees(500) seed(42) ///
    llenable lllambda(0.1) llvars(x1 x2 x3 x4 x5)
keep stata_pred r_pred
export delimited "`OUTDIR'/test02_stata.csv", replace
display "  Test 02 done."

* ============================================================
* Test 03: llvars(x1 x2) — subset LL correction
* ============================================================
display ""
display "--- Test 03: llvars(x1 x2) ---"
import delimited "`OUTDIR'/test03_data.csv", clear
grf_ll_regression_forest y x1 x2 x3 x4 x5, gen(stata_pred) ntrees(500) seed(42) ///
    lllambda(0.1) llvars(x1 x2)
keep stata_pred r_pred
export delimited "`OUTDIR'/test03_stata.csv", replace
display "  Test 03 done."

* ============================================================
* Test 04: llsplitvars(x1) — restricted split vars for LL
* In Stata: llenable llvars(x1) controls split vars
* ============================================================
display ""
display "--- Test 04: llvars(x1) with llenable (split vars restricted) ---"
import delimited "`OUTDIR'/test04_data.csv", clear
grf_ll_regression_forest y x1 x2 x3 x4 x5, gen(stata_pred) ntrees(500) seed(42) ///
    llenable llvars(x1) lllambda(0.1)
keep stata_pred r_pred
export delimited "`OUTDIR'/test04_stata.csv", replace
display "  Test 04 done."

* ============================================================
* Test 05: lllambda=0.01
* ============================================================
display ""
display "--- Test 05: lllambda=0.01 ---"
import delimited "`OUTDIR'/test05_data.csv", clear
grf_ll_regression_forest y x1 x2 x3 x4 x5, gen(stata_pred) ntrees(500) seed(42) ///
    lllambda(0.01) llvars(x1 x2 x3 x4 x5)
keep stata_pred r_pred
export delimited "`OUTDIR'/test05_stata.csv", replace
display "  Test 05 done."

* ============================================================
* Test 06: lllambda=1.0
* ============================================================
display ""
display "--- Test 06: lllambda=1.0 ---"
import delimited "`OUTDIR'/test06_data.csv", clear
grf_ll_regression_forest y x1 x2 x3 x4 x5, gen(stata_pred) ntrees(500) seed(42) ///
    lllambda(1.0) llvars(x1 x2 x3 x4 x5)
keep stata_pred r_pred
export delimited "`OUTDIR'/test06_stata.csv", replace
display "  Test 06 done."

* ============================================================
* Test 07: lllambda=10.0
* ============================================================
display ""
display "--- Test 07: lllambda=10.0 ---"
import delimited "`OUTDIR'/test07_data.csv", clear
grf_ll_regression_forest y x1 x2 x3 x4 x5, gen(stata_pred) ntrees(500) seed(42) ///
    lllambda(10.0) llvars(x1 x2 x3 x4 x5)
keep stata_pred r_pred r_pred_noll
export delimited "`OUTDIR'/test07_stata.csv", replace
display "  Test 07 done."

* ============================================================
* Test 08: llweightpenalty
* ============================================================
display ""
display "--- Test 08: llweightpenalty ---"
import delimited "`OUTDIR'/test08_data.csv", clear
grf_ll_regression_forest y x1 x2 x3 x4 x5, gen(stata_pred) ntrees(500) seed(42) ///
    lllambda(0.1) llvars(x1 x2 x3 x4 x5) llweightpenalty
keep stata_pred r_pred
export delimited "`OUTDIR'/test08_stata.csv", replace
display "  Test 08 done."

* ============================================================
* Test 09: llcutoff=3
* ============================================================
display ""
display "--- Test 09: llcutoff=3 ---"
import delimited "`OUTDIR'/test09_data.csv", clear
grf_ll_regression_forest y x1 x2 x3 x4 x5, gen(stata_pred) ntrees(500) seed(42) ///
    llenable lllambda(0.1) llvars(x1 x2 x3 x4 x5) llcutoff(3)
keep stata_pred r_pred
export delimited "`OUTDIR'/test09_stata.csv", replace
display "  Test 09 done."

* ============================================================
* Test 10: cluster()
* ============================================================
display ""
display "--- Test 10: cluster() ---"
import delimited "`OUTDIR'/test10_data.csv", clear
grf_ll_regression_forest y x1 x2 x3 x4 x5, gen(stata_pred) ntrees(500) seed(42) ///
    lllambda(0.1) llvars(x1 x2 x3 x4 x5) cluster(clust)
keep stata_pred r_pred
export delimited "`OUTDIR'/test10_stata.csv", replace
display "  Test 10 done."

* ============================================================
* Test 11: weights()
* ============================================================
display ""
display "--- Test 11: weights() ---"
import delimited "`OUTDIR'/test11_data.csv", clear
grf_ll_regression_forest y x1 x2 x3 x4 x5, gen(stata_pred) ntrees(500) seed(42) ///
    lllambda(0.1) llvars(x1 x2 x3 x4 x5) weights(wt)
keep stata_pred r_pred
export delimited "`OUTDIR'/test11_stata.csv", replace
display "  Test 11 done."

* ============================================================
* Test 12: nohonesty
* ============================================================
display ""
display "--- Test 12: nohonesty ---"
import delimited "`OUTDIR'/test12_data.csv", clear
grf_ll_regression_forest y x1 x2 x3 x4 x5, gen(stata_pred) ntrees(500) seed(42) ///
    nohonesty lllambda(0.1) llvars(x1 x2 x3 x4 x5)
keep stata_pred r_pred
export delimited "`OUTDIR'/test12_stata.csv", replace
display "  Test 12 done."

* ============================================================
* Test 13: mtry=2
* ============================================================
display ""
display "--- Test 13: mtry=2 ---"
import delimited "`OUTDIR'/test13_data.csv", clear
grf_ll_regression_forest y x1 x2 x3 x4 x5, gen(stata_pred) ntrees(500) seed(42) ///
    mtry(2) lllambda(0.1) llvars(x1 x2 x3 x4 x5)
keep stata_pred r_pred
export delimited "`OUTDIR'/test13_stata.csv", replace
display "  Test 13 done."

* ============================================================
* Test 14: minnodesize=20
* ============================================================
display ""
display "--- Test 14: minnodesize=20 ---"
import delimited "`OUTDIR'/test14_data.csv", clear
grf_ll_regression_forest y x1 x2 x3 x4 x5, gen(stata_pred) ntrees(500) seed(42) ///
    minnodesize(20) lllambda(0.1) llvars(x1 x2 x3 x4 x5)
keep stata_pred r_pred
export delimited "`OUTDIR'/test14_stata.csv", replace
display "  Test 14 done."

* ============================================================
* Test 15: Combined llvars(x1 x2) + lllambda=0.5 + llweightpenalty
* ============================================================
display ""
display "--- Test 15: Combined (llvars + lllambda=0.5 + llweightpenalty) ---"
import delimited "`OUTDIR'/test15_data.csv", clear
grf_ll_regression_forest y x1 x2 x3 x4 x5, gen(stata_pred) ntrees(500) seed(42) ///
    llvars(x1 x2) lllambda(0.5) llweightpenalty
keep stata_pred r_pred
export delimited "`OUTDIR'/test15_stata.csv", replace
display "  Test 15 done."

* ============================================================
* Test 16: LL vs non-LL (linear DGP)
* ============================================================
display ""
display "--- Test 16: LL vs non-LL (linear DGP) ---"
import delimited "`OUTDIR'/test16_data.csv", clear
* LL predictions
grf_ll_regression_forest y x1 x2 x3 x4 x5, gen(stata_pred_ll) ntrees(500) seed(42) ///
    lllambda(0.1) llvars(x1 x2 x3 x4 x5)
* Standard regression forest (no LL)
grf_regression_forest y x1 x2 x3 x4 x5, gen(stata_pred_noll) ntrees(500) seed(42)
quietly gen mse_ll   = (stata_pred_ll   - y)^2
quietly gen mse_noll = (stata_pred_noll - y)^2
quietly summarize mse_ll
local stata_mse_ll = r(mean)
quietly summarize mse_noll
local stata_mse_noll = r(mean)
display "  Stata MSE LL:    " %9.4f `stata_mse_ll'
display "  Stata MSE no-LL: " %9.4f `stata_mse_noll'
keep stata_pred_ll stata_pred_noll r_pred r_pred_noll r_pred_std true_mu y
export delimited "`OUTDIR'/test16_stata.csv", replace
display "  Test 16 done."

* ============================================================
* Test 17: Nonlinear DGP (sin + quadratic)
* ============================================================
display ""
display "--- Test 17: Nonlinear DGP ---"
import delimited "`OUTDIR'/test17_data.csv", clear
grf_ll_regression_forest y x1 x2 x3 x4 x5, gen(stata_pred) ntrees(500) seed(42) ///
    lllambda(0.1) llvars(x1 x2 x3 x4 x5)
keep stata_pred r_pred
export delimited "`OUTDIR'/test17_stata.csv", replace
display "  Test 17 done."

* ============================================================
* Test 18: Pure Linear DGP (Y = sum(X))
* ============================================================
display ""
display "--- Test 18: Pure Linear DGP ---"
import delimited "`OUTDIR'/test18_data.csv", clear
grf_ll_regression_forest y x1 x2 x3 x4 x5, gen(stata_pred) ntrees(500) seed(42) ///
    lllambda(0.1) llvars(x1 x2 x3 x4 x5)
grf_regression_forest y x1 x2 x3 x4 x5, gen(stata_pred_std) ntrees(500) seed(42)
quietly gen mse_ll  = (stata_pred - y)^2
quietly gen mse_std = (stata_pred_std - y)^2
quietly summarize mse_ll
display "  Stata MSE LL:    " %9.4f r(mean)
quietly summarize mse_std
display "  Stata MSE std:   " %9.4f r(mean)
keep stata_pred stata_pred_std r_pred r_pred_std y
export delimited "`OUTDIR'/test18_stata.csv", replace
display "  Test 18 done."

display ""
display "============================================================"
display " All 18 Stata tests complete"
display "============================================================"
