* Test 14: Compare boosted vs plain regression_forest in OOB MSE
adopath ++ /tmp/grf_stata
import delimited "/tmp/grf_stata/tests/fidelity_reports/10_boosted/test14_data.csv", clear

* Run boosted regression forest
grf_boosted_regression_forest y x1 x2 x3 x4 x5, gen(brf_pred) ntrees(500) seed(42)
local brf_steps = e(boost_steps)

* Run plain regression forest
grf_regression_forest y x1 x2 x3 x4 x5, gen(rf_pred) ntrees(500) seed(42)

* Compute MSE vs true mu
gen brf_sq_err = (brf_pred - true_mu)^2
gen rf_sq_err  = (rf_pred  - true_mu)^2
quietly summarize brf_sq_err
local brf_mse = r(mean)
quietly summarize rf_sq_err
local rf_mse = r(mean)

display "Regression forest MSE vs true mu: " %8.6f `rf_mse'
display "Boosted forest MSE vs true mu:    " %8.6f `brf_mse'
display "Boost steps used:                  " `brf_steps'
display "Improvement: " %6.2f ((`rf_mse' - `brf_mse') / `rf_mse' * 100) "%"

* Correlation with R predictions
correlate brf_pred r_pred_brf
local brf_corr = r(rho)
correlate rf_pred r_pred_rf
local rf_corr = r(rho)
display "Boosted: Stata-vs-R corr = " %8.6f `brf_corr'
display "Plain RF: Stata-vs-R corr = " %8.6f `rf_corr'

keep brf_pred rf_pred r_pred_brf r_pred_rf true_mu brf_sq_err rf_sq_err
export delimited "/tmp/grf_stata/tests/fidelity_reports/10_boosted/test14_stata.csv", replace
