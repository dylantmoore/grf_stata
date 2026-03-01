* Test 09: estimate.variance
adopath ++ /tmp/grf_stata
import delimited "/tmp/grf_stata/tests/fidelity_reports/01_regression/test09_data.csv", clear
grf_regression_forest y x1 x2 x3 x4 x5, gen(stata_pred) ntrees(500) seed(42) estimatevariance vargenerate(stata_var) cigroupsize(2)
keep stata_pred stata_var r_pred r_var
export delimited "/tmp/grf_stata/tests/fidelity_reports/01_regression/test09_stata.csv", replace
