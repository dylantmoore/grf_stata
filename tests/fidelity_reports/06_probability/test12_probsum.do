* Test 12: Probabilities sum to 1 check (3-class)
adopath ++ /tmp/grf_stata
import delimited "/tmp/grf_stata/tests/fidelity_reports/06_probability/test12_data.csv", clear
grf_probability_forest y x1 x2 x3 x4 x5, gen(stata_pred) ntrees(500) seed(42)

* Compute sum of class probabilities
gen stata_sum = stata_pred_c0 + stata_pred_c1 + stata_pred_c2
gen stata_sum_dev = abs(stata_sum - 1)
summarize stata_sum_dev

keep stata_pred_c0 stata_pred_c1 stata_pred_c2 stata_sum stata_sum_dev r_pred_c0 r_pred_c1 r_pred_c2 r_sum
export delimited "/tmp/grf_stata/tests/fidelity_reports/06_probability/test12_stata.csv", replace
