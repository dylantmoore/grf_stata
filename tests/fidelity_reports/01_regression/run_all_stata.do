* Master do-file: run all regression forest fidelity tests
adopath ++ /tmp/grf_stata

* ---- Test 02: nohonesty ----
import delimited "/tmp/grf_stata/tests/fidelity_reports/01_regression/test02_data.csv", clear
grf_regression_forest y x1 x2 x3 x4 x5, gen(stata_pred) ntrees(500) seed(42) nohonesty
keep stata_pred r_pred
export delimited "/tmp/grf_stata/tests/fidelity_reports/01_regression/test02_stata.csv", replace

* ---- Test 03: mtry=2 ----
import delimited "/tmp/grf_stata/tests/fidelity_reports/01_regression/test03_data.csv", clear
grf_regression_forest y x1 x2 x3 x4 x5, gen(stata_pred) ntrees(500) seed(42) mtry(2)
keep stata_pred r_pred
export delimited "/tmp/grf_stata/tests/fidelity_reports/01_regression/test03_stata.csv", replace

* ---- Test 04: minnodesize=20 ----
import delimited "/tmp/grf_stata/tests/fidelity_reports/01_regression/test04_data.csv", clear
grf_regression_forest y x1 x2 x3 x4 x5, gen(stata_pred) ntrees(500) seed(42) minnodesize(20)
keep stata_pred r_pred
export delimited "/tmp/grf_stata/tests/fidelity_reports/01_regression/test04_stata.csv", replace

* ---- Test 05: samplefrac=0.3 ----
import delimited "/tmp/grf_stata/tests/fidelity_reports/01_regression/test05_data.csv", clear
grf_regression_forest y x1 x2 x3 x4 x5, gen(stata_pred) ntrees(500) seed(42) samplefrac(0.3)
keep stata_pred r_pred
export delimited "/tmp/grf_stata/tests/fidelity_reports/01_regression/test05_stata.csv", replace

* ---- Test 06: honestyfrac=0.7 ----
import delimited "/tmp/grf_stata/tests/fidelity_reports/01_regression/test06_data.csv", clear
grf_regression_forest y x1 x2 x3 x4 x5, gen(stata_pred) ntrees(500) seed(42) honestyfrac(0.7)
keep stata_pred r_pred
export delimited "/tmp/grf_stata/tests/fidelity_reports/01_regression/test06_stata.csv", replace

* ---- Test 07: alpha=0.15 ----
import delimited "/tmp/grf_stata/tests/fidelity_reports/01_regression/test07_data.csv", clear
grf_regression_forest y x1 x2 x3 x4 x5, gen(stata_pred) ntrees(500) seed(42) alpha(0.15)
keep stata_pred r_pred
export delimited "/tmp/grf_stata/tests/fidelity_reports/01_regression/test07_stata.csv", replace

* ---- Test 08: imbalancepenalty=1.0 ----
import delimited "/tmp/grf_stata/tests/fidelity_reports/01_regression/test08_data.csv", clear
grf_regression_forest y x1 x2 x3 x4 x5, gen(stata_pred) ntrees(500) seed(42) imbalancepenalty(1.0)
keep stata_pred r_pred
export delimited "/tmp/grf_stata/tests/fidelity_reports/01_regression/test08_stata.csv", replace

* ---- Test 09: estimatevariance ----
import delimited "/tmp/grf_stata/tests/fidelity_reports/01_regression/test09_data.csv", clear
grf_regression_forest y x1 x2 x3 x4 x5, gen(stata_pred) ntrees(500) seed(42) estimatevariance vargenerate(stata_var) cigroupsize(2)
keep stata_pred stata_var r_pred r_var
export delimited "/tmp/grf_stata/tests/fidelity_reports/01_regression/test09_stata.csv", replace

* ---- Test 10: cluster ----
import delimited "/tmp/grf_stata/tests/fidelity_reports/01_regression/test10_data.csv", clear
grf_regression_forest y x1 x2 x3 x4 x5, gen(stata_pred) ntrees(500) seed(42) cluster(cluster_id)
keep stata_pred r_pred
export delimited "/tmp/grf_stata/tests/fidelity_reports/01_regression/test10_stata.csv", replace

* ---- Test 11: weights ----
import delimited "/tmp/grf_stata/tests/fidelity_reports/01_regression/test11_data.csv", clear
grf_regression_forest y x1 x2 x3 x4 x5, gen(stata_pred) ntrees(500) seed(42) weights(wts)
keep stata_pred r_pred
export delimited "/tmp/grf_stata/tests/fidelity_reports/01_regression/test11_stata.csv", replace

* ---- Test 12: equalizeclusterweights ----
import delimited "/tmp/grf_stata/tests/fidelity_reports/01_regression/test12_data.csv", clear
grf_regression_forest y x1 x2 x3 x4 x5, gen(stata_pred) ntrees(500) seed(42) cluster(cluster_id) equalizeclusterweights
keep stata_pred r_pred
export delimited "/tmp/grf_stata/tests/fidelity_reports/01_regression/test12_stata.csv", replace

* ---- Test 13: nomia ----
import delimited "/tmp/grf_stata/tests/fidelity_reports/01_regression/test13_data.csv", clear
grf_regression_forest y x1 x2 x3 x4 x5, gen(stata_pred) ntrees(500) seed(42) nomia
keep stata_pred r_pred
export delimited "/tmp/grf_stata/tests/fidelity_reports/01_regression/test13_stata.csv", replace

* ---- Test 14: cigroupsize=2 ----
import delimited "/tmp/grf_stata/tests/fidelity_reports/01_regression/test14_data.csv", clear
grf_regression_forest y x1 x2 x3 x4 x5, gen(stata_pred) ntrees(500) seed(42) cigroupsize(2)
keep stata_pred r_pred
export delimited "/tmp/grf_stata/tests/fidelity_reports/01_regression/test14_stata.csv", replace

* ---- Test 15: nohonesty + mtry=3 + minnodesize=15 ----
import delimited "/tmp/grf_stata/tests/fidelity_reports/01_regression/test15_data.csv", clear
grf_regression_forest y x1 x2 x3 x4 x5, gen(stata_pred) ntrees(500) seed(42) nohonesty mtry(3) minnodesize(15)
keep stata_pred r_pred
export delimited "/tmp/grf_stata/tests/fidelity_reports/01_regression/test15_stata.csv", replace

* ---- Test 16: cluster + weights + estimatevariance ----
import delimited "/tmp/grf_stata/tests/fidelity_reports/01_regression/test16_data.csv", clear
grf_regression_forest y x1 x2 x3 x4 x5, gen(stata_pred) ntrees(500) seed(42) cluster(cluster_id) weights(wts) estimatevariance vargenerate(stata_var) cigroupsize(2)
keep stata_pred stata_var r_pred r_var
export delimited "/tmp/grf_stata/tests/fidelity_reports/01_regression/test16_stata.csv", replace

* ---- Test 17: large p=20 ----
import delimited "/tmp/grf_stata/tests/fidelity_reports/01_regression/test17_data.csv", clear
grf_regression_forest y x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14 x15 x16 x17 x18 x19 x20, gen(stata_pred) ntrees(500) seed(42)
keep stata_pred r_pred
export delimited "/tmp/grf_stata/tests/fidelity_reports/01_regression/test17_stata.csv", replace

* ---- Test 18a: ntrees=100 ----
import delimited "/tmp/grf_stata/tests/fidelity_reports/01_regression/test18a_data.csv", clear
grf_regression_forest y x1 x2 x3 x4 x5, gen(stata_pred) ntrees(100) seed(42)
keep stata_pred r_pred
export delimited "/tmp/grf_stata/tests/fidelity_reports/01_regression/test18a_stata.csv", replace

* ---- Test 18b: ntrees=2000 ----
import delimited "/tmp/grf_stata/tests/fidelity_reports/01_regression/test18b_data.csv", clear
grf_regression_forest y x1 x2 x3 x4 x5, gen(stata_pred) ntrees(2000) seed(42)
keep stata_pred r_pred
export delimited "/tmp/grf_stata/tests/fidelity_reports/01_regression/test18b_stata.csv", replace

* ---- Test 19: missing data (MIA on, default) ----
import delimited "/tmp/grf_stata/tests/fidelity_reports/01_regression/test19_data.csv", clear
grf_regression_forest y x1 x2 x3 x4 x5, gen(stata_pred) ntrees(500) seed(42)
keep stata_pred r_pred
export delimited "/tmp/grf_stata/tests/fidelity_reports/01_regression/test19_stata.csv", replace

* ---- Test 20: seed reproducibility ----
import delimited "/tmp/grf_stata/tests/fidelity_reports/01_regression/test20_data.csv", clear
grf_regression_forest y x1 x2 x3 x4 x5, gen(stata_pred1) ntrees(500) seed(123)
grf_regression_forest y x1 x2 x3 x4 x5, gen(stata_pred2) ntrees(500) seed(123)
keep stata_pred1 stata_pred2 r_pred
export delimited "/tmp/grf_stata/tests/fidelity_reports/01_regression/test20_stata.csv", replace
