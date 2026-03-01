* Test 04: Extreme quantiles (0.01, 0.99)
adopath ++ /tmp/grf_stata
import delimited "/tmp/grf_stata/tests/fidelity_reports/05_quantile/test04_data.csv", clear
grf_quantile_forest y x1 x2 x3 x4 x5, gen(stata_pred) quantiles(0.01 0.99) ntrees(500) seed(42)
keep stata_pred_q1 stata_pred_q99 r_q1 r_q99
export delimited "/tmp/grf_stata/tests/fidelity_reports/05_quantile/test04_stata.csv", replace
