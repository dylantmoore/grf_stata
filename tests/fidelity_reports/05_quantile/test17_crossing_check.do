* Test 17: Quantile crossing check — q10 <= q50 <= q90 for all obs?
adopath ++ /tmp/grf_stata
import delimited "/tmp/grf_stata/tests/fidelity_reports/05_quantile/test17_data.csv", clear
grf_quantile_forest y x1 x2 x3 x4 x5, gen(stata_pred) quantiles(0.1 0.5 0.9) ntrees(500) seed(42)
* Count violations
gen cross_10_50 = (stata_pred_q10 > stata_pred_q50)
gen cross_50_90 = (stata_pred_q50 > stata_pred_q90)
quietly summarize cross_10_50
local n_cross_10_50 = r(sum)
quietly summarize cross_50_90
local n_cross_50_90 = r(sum)
display "Stata q10>q50 violations: " `n_cross_10_50'
display "Stata q50>q90 violations: " `n_cross_50_90'
keep stata_pred_q10 stata_pred_q50 stata_pred_q90 r_q10 r_q50 r_q90 cross_10_50 cross_50_90
export delimited "/tmp/grf_stata/tests/fidelity_reports/05_quantile/test17_stata.csv", replace
