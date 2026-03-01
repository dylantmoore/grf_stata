* Test 11: sample weights
adopath ++ /tmp/grf_stata
import delimited "/tmp/grf_stata/tests/fidelity_reports/01_regression/test11_data.csv", clear
grf_regression_forest y x1 x2 x3 x4 x5, gen(stata_pred) ntrees(500) seed(42) weights(wts)
keep stata_pred r_pred
export delimited "/tmp/grf_stata/tests/fidelity_reports/01_regression/test11_stata.csv", replace
