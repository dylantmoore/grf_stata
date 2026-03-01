* Test 12: mtry=2 — restricted splitting
adopath ++ /tmp/grf_stata
import delimited "/tmp/grf_stata/tests/fidelity_reports/05_quantile/test12_data.csv", clear
grf_quantile_forest y x1 x2 x3 x4 x5, gen(stata_pred) quantiles(0.1 0.5 0.9) ntrees(500) seed(42) mtry(2)
keep stata_pred_q10 stata_pred_q50 stata_pred_q90 r_q10 r_q50 r_q90
export delimited "/tmp/grf_stata/tests/fidelity_reports/05_quantile/test12_stata.csv", replace
