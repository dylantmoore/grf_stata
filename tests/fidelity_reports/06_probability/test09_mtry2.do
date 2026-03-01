* Test 09: mtry=2
adopath ++ /tmp/grf_stata
import delimited "/tmp/grf_stata/tests/fidelity_reports/06_probability/test09_data.csv", clear
grf_probability_forest y x1 x2 x3 x4 x5, gen(stata_pred) ntrees(500) seed(42) mtry(2)
keep stata_pred_c0 stata_pred_c1 r_pred_c0 r_pred_c1
export delimited "/tmp/grf_stata/tests/fidelity_reports/06_probability/test09_stata.csv", replace
