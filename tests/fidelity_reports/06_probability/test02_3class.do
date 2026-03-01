* Test 02: 3-class multinomial
adopath ++ /tmp/grf_stata
import delimited "/tmp/grf_stata/tests/fidelity_reports/06_probability/test02_data.csv", clear
grf_probability_forest y x1 x2 x3 x4 x5, gen(stata_pred) ntrees(500) seed(42)
keep stata_pred_c0 stata_pred_c1 stata_pred_c2 r_pred_c0 r_pred_c1 r_pred_c2
export delimited "/tmp/grf_stata/tests/fidelity_reports/06_probability/test02_stata.csv", replace
