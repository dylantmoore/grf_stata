* Test 07: alpha=0.15
adopath ++ /tmp/grf_stata
import delimited "/tmp/grf_stata/tests/fidelity_reports/01_regression/test07_data.csv", clear
grf_regression_forest y x1 x2 x3 x4 x5, gen(stata_pred) ntrees(500) seed(42) alpha(0.15)
keep stata_pred r_pred
export delimited "/tmp/grf_stata/tests/fidelity_reports/01_regression/test07_stata.csv", replace
