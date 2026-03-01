* run_all_stata.do
* Run all 14 multi_regression_forest fidelity tests and export Stata predictions
adopath ++ /tmp/grf_stata

local OUTDIR "/tmp/grf_stata/tests/fidelity_reports/13_multi_regression"

* ==============================================================
* Test 01: 2 outcomes (Y1, Y2) -- default options
* ==============================================================
import delimited "`OUTDIR'/test01_data.csv", clear
grf_multi_regression_forest y1 y2 x1 x2 x3 x4 x5, ///
    gen(stata_pred) ndep(2) ntrees(500) seed(42)
keep stata_pred_y1 stata_pred_y2 r_pred_y1 r_pred_y2
export delimited "`OUTDIR'/test01_stata.csv", replace

* ==============================================================
* Test 02: 3 outcomes
* ==============================================================
import delimited "`OUTDIR'/test02_data.csv", clear
grf_multi_regression_forest y1 y2 y3 x1 x2 x3 x4 x5, ///
    gen(stata_pred) ndep(3) ntrees(500) seed(42)
keep stata_pred_y1 stata_pred_y2 stata_pred_y3 r_pred_y1 r_pred_y2 r_pred_y3
export delimited "`OUTDIR'/test02_stata.csv", replace

* ==============================================================
* Test 03: 5 outcomes
* ==============================================================
import delimited "`OUTDIR'/test03_data.csv", clear
grf_multi_regression_forest y1 y2 y3 y4 y5 x1 x2 x3 x4 x5, ///
    gen(stata_pred) ndep(5) ntrees(500) seed(42)
keep stata_pred_y1 stata_pred_y2 stata_pred_y3 stata_pred_y4 stata_pred_y5 ///
     r_pred_y1 r_pred_y2 r_pred_y3 r_pred_y4 r_pred_y5
export delimited "`OUTDIR'/test03_stata.csv", replace

* ==============================================================
* Test 04: Correlated outcomes
* ==============================================================
import delimited "`OUTDIR'/test04_data.csv", clear
grf_multi_regression_forest y1 y2 x1 x2 x3 x4 x5, ///
    gen(stata_pred) ndep(2) ntrees(500) seed(42)
keep stata_pred_y1 stata_pred_y2 r_pred_y1 r_pred_y2
export delimited "`OUTDIR'/test04_stata.csv", replace

* ==============================================================
* Test 05: Independent outcomes
* ==============================================================
import delimited "`OUTDIR'/test05_data.csv", clear
grf_multi_regression_forest y1 y2 x1 x2 x3 x4 x5, ///
    gen(stata_pred) ndep(2) ntrees(500) seed(42)
keep stata_pred_y1 stata_pred_y2 r_pred_y1 r_pred_y2
export delimited "`OUTDIR'/test05_stata.csv", replace

* ==============================================================
* Test 06: cluster()
* ==============================================================
import delimited "`OUTDIR'/test06_data.csv", clear
grf_multi_regression_forest y1 y2 x1 x2 x3 x4 x5, ///
    gen(stata_pred) ndep(2) ntrees(500) seed(42) cluster(cluster_id)
keep stata_pred_y1 stata_pred_y2 r_pred_y1 r_pred_y2
export delimited "`OUTDIR'/test06_stata.csv", replace

* ==============================================================
* Test 07: weights()
* ==============================================================
import delimited "`OUTDIR'/test07_data.csv", clear
grf_multi_regression_forest y1 y2 x1 x2 x3 x4 x5, ///
    gen(stata_pred) ndep(2) ntrees(500) seed(42) weights(wt)
keep stata_pred_y1 stata_pred_y2 r_pred_y1 r_pred_y2
export delimited "`OUTDIR'/test07_stata.csv", replace

* ==============================================================
* Test 08: nohonesty
* ==============================================================
import delimited "`OUTDIR'/test08_data.csv", clear
grf_multi_regression_forest y1 y2 x1 x2 x3 x4 x5, ///
    gen(stata_pred) ndep(2) ntrees(500) seed(42) nohonesty
keep stata_pred_y1 stata_pred_y2 r_pred_y1 r_pred_y2
export delimited "`OUTDIR'/test08_stata.csv", replace

* ==============================================================
* Test 09: mtry=2
* ==============================================================
import delimited "`OUTDIR'/test09_data.csv", clear
grf_multi_regression_forest y1 y2 x1 x2 x3 x4 x5, ///
    gen(stata_pred) ndep(2) ntrees(500) seed(42) mtry(2)
keep stata_pred_y1 stata_pred_y2 r_pred_y1 r_pred_y2
export delimited "`OUTDIR'/test09_stata.csv", replace

* ==============================================================
* Test 10: minnodesize=20
* ==============================================================
import delimited "`OUTDIR'/test10_data.csv", clear
grf_multi_regression_forest y1 y2 x1 x2 x3 x4 x5, ///
    gen(stata_pred) ndep(2) ntrees(500) seed(42) minnodesize(20)
keep stata_pred_y1 stata_pred_y2 r_pred_y1 r_pred_y2
export delimited "`OUTDIR'/test10_stata.csv", replace

* ==============================================================
* Test 11: samplefrac=0.3
* ==============================================================
import delimited "`OUTDIR'/test11_data.csv", clear
grf_multi_regression_forest y1 y2 x1 x2 x3 x4 x5, ///
    gen(stata_pred) ndep(2) ntrees(500) seed(42) samplefrac(0.3)
keep stata_pred_y1 stata_pred_y2 r_pred_y1 r_pred_y2
export delimited "`OUTDIR'/test11_stata.csv", replace

* ==============================================================
* Test 12: Combined cluster + weights + nohonesty
* ==============================================================
import delimited "`OUTDIR'/test12_data.csv", clear
grf_multi_regression_forest y1 y2 x1 x2 x3 x4 x5, ///
    gen(stata_pred) ndep(2) ntrees(500) seed(42) ///
    cluster(cluster_id) weights(wt) nohonesty
keep stata_pred_y1 stata_pred_y2 r_pred_y1 r_pred_y2
export delimited "`OUTDIR'/test12_stata.csv", replace

* ==============================================================
* Test 13: Linear outcomes
* ==============================================================
import delimited "`OUTDIR'/test13_data.csv", clear
grf_multi_regression_forest y1 y2 x1 x2 x3 x4 x5, ///
    gen(stata_pred) ndep(2) ntrees(500) seed(42)
keep stata_pred_y1 stata_pred_y2 r_pred_y1 r_pred_y2
export delimited "`OUTDIR'/test13_stata.csv", replace

* ==============================================================
* Test 14: Nonlinear outcomes
* ==============================================================
import delimited "`OUTDIR'/test14_data.csv", clear
grf_multi_regression_forest y1 y2 x1 x2 x3 x4 x5, ///
    gen(stata_pred) ndep(2) ntrees(500) seed(42)
keep stata_pred_y1 stata_pred_y2 r_pred_y1 r_pred_y2
export delimited "`OUTDIR'/test14_stata.csv", replace

exit
