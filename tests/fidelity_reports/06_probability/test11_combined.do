* Test 11: Combined options: cluster + weights + nohonesty
adopath ++ /tmp/grf_stata
import delimited "/tmp/grf_stata/tests/fidelity_reports/06_probability/test11_data.csv", clear
grf_probability_forest y x1 x2 x3 x4 x5, gen(stata_pred) ntrees(500) seed(42) ///
    cluster(cluster_id) weights(wt) nohonesty
keep stata_pred_c0 stata_pred_c1 r_pred_c0 r_pred_c1
export delimited "/tmp/grf_stata/tests/fidelity_reports/06_probability/test11_stata.csv", replace
