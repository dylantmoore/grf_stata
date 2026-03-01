* Test 20: Seed reproducibility
adopath ++ /tmp/grf_stata
import delimited "/tmp/grf_stata/tests/fidelity_reports/01_regression/test20_data.csv", clear
* Run twice with same seed
grf_regression_forest y x1 x2 x3 x4 x5, gen(stata_pred1) ntrees(500) seed(123)
grf_regression_forest y x1 x2 x3 x4 x5, gen(stata_pred2) ntrees(500) seed(123)
* Compare them
keep stata_pred1 stata_pred2 r_pred
export delimited "/tmp/grf_stata/tests/fidelity_reports/01_regression/test20_stata.csv", replace
