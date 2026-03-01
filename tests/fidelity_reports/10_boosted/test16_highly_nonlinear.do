* Test 16: Highly nonlinear data - boosting should help
adopath ++ /tmp/grf_stata
import delimited "/tmp/grf_stata/tests/fidelity_reports/10_boosted/test16_data.csv", clear
grf_boosted_regression_forest y x1 x2 x3 x4 x5, gen(stata_pred) ntrees(500) seed(42)
display "Boost steps used: " e(boost_steps)
correlate stata_pred r_pred
display "Stata-vs-R corr: " r(rho)
gen sq_err = (stata_pred - true_mu)^2
quietly summarize sq_err
display "Stata BRF MSE vs true: " r(mean)
keep stata_pred r_pred true_mu
export delimited "/tmp/grf_stata/tests/fidelity_reports/10_boosted/test16_stata.csv", replace
