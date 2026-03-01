* Run tests 19 and 20
adopath ++ /tmp/grf_stata

* ---- Test 19: missing data (MIA on, default) ----
import delimited "/tmp/grf_stata/tests/fidelity_reports/01_regression/test19_data.csv", clear
* Force numeric conversion where needed
destring x1 x2 x3 x4 x5 y r_pred, replace force
grf_regression_forest y x1 x2 x3 x4 x5, gen(stata_pred) ntrees(500) seed(42)
keep stata_pred r_pred
export delimited "/tmp/grf_stata/tests/fidelity_reports/01_regression/test19_stata.csv", replace

* ---- Test 20: seed reproducibility ----
import delimited "/tmp/grf_stata/tests/fidelity_reports/01_regression/test20_data.csv", clear
grf_regression_forest y x1 x2 x3 x4 x5, gen(stata_pred1) ntrees(500) seed(123)
grf_regression_forest y x1 x2 x3 x4 x5, gen(stata_pred2) ntrees(500) seed(123)
keep stata_pred1 stata_pred2 r_pred
export delimited "/tmp/grf_stata/tests/fidelity_reports/01_regression/test20_stata.csv", replace
