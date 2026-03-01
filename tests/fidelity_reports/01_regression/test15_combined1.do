* Test 15: Combined: nohonesty + mtry=3 + minnodesize=15
adopath ++ /tmp/grf_stata
import delimited "/tmp/grf_stata/tests/fidelity_reports/01_regression/test15_data.csv", clear
grf_regression_forest y x1 x2 x3 x4 x5, gen(stata_pred) ntrees(500) seed(42) nohonesty mtry(3) minnodesize(15)
keep stata_pred r_pred
export delimited "/tmp/grf_stata/tests/fidelity_reports/01_regression/test15_stata.csv", replace
