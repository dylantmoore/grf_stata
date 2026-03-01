* Test 05: sample.fraction=0.3
adopath ++ /tmp/grf_stata
import delimited "/tmp/grf_stata/tests/fidelity_reports/01_regression/test05_data.csv", clear
grf_regression_forest y x1 x2 x3 x4 x5, gen(stata_pred) ntrees(500) seed(42) samplefrac(0.3)
keep stata_pred r_pred
export delimited "/tmp/grf_stata/tests/fidelity_reports/01_regression/test05_stata.csv", replace
