* Test 10: With clusters
adopath ++ /tmp/grf_stata
import delimited "/tmp/grf_stata/tests/fidelity_reports/10_boosted/test10_data.csv", clear
grf_boosted_regression_forest y x1 x2 x3 x4 x5, gen(stata_pred) ntrees(500) seed(42) booststeps(2) cluster(cluster_id)
display "Boost steps: " e(boost_steps)
keep stata_pred r_pred true_mu cluster_id
export delimited "/tmp/grf_stata/tests/fidelity_reports/10_boosted/test10_stata.csv", replace
