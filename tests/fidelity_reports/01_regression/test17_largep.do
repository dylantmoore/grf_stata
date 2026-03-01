* Test 17: Large p (p=20)
adopath ++ /tmp/grf_stata
import delimited "/tmp/grf_stata/tests/fidelity_reports/01_regression/test17_data.csv", clear
grf_regression_forest y x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14 x15 x16 x17 x18 x19 x20, gen(stata_pred) ntrees(500) seed(42)
keep stata_pred r_pred
export delimited "/tmp/grf_stata/tests/fidelity_reports/01_regression/test17_stata.csv", replace
