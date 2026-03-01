* Test 02: Single quantile (0.5 - median only)
adopath ++ /tmp/grf_stata
import delimited "/tmp/grf_stata/tests/fidelity_reports/05_quantile/test02_data.csv", clear
grf_quantile_forest y x1 x2 x3 x4 x5, gen(stata_pred) quantiles(0.5) ntrees(500) seed(42)
keep stata_pred_q50 r_q50
export delimited "/tmp/grf_stata/tests/fidelity_reports/05_quantile/test02_stata.csv", replace
