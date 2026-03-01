* Test 07: boosterrorreduction(0.99) - less aggressive stopping
adopath ++ /tmp/grf_stata
import delimited "/tmp/grf_stata/tests/fidelity_reports/10_boosted/test07_data.csv", clear
grf_boosted_regression_forest y x1 x2 x3 x4 x5, gen(stata_pred) ntrees(500) seed(42) boosterrorreduction(0.99)
display "Boost steps: " e(boost_steps)
keep stata_pred r_pred true_mu
export delimited "/tmp/grf_stata/tests/fidelity_reports/10_boosted/test07_stata.csv", replace
