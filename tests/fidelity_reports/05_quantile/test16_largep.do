* Test 16: Large p (p=20) — many predictors
adopath ++ /tmp/grf_stata
import delimited "/tmp/grf_stata/tests/fidelity_reports/05_quantile/test16_data.csv", clear
grf_quantile_forest y x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14 x15 x16 x17 x18 x19 x20, ///
    gen(stata_pred) quantiles(0.1 0.5 0.9) ntrees(500) seed(42)
keep stata_pred_q10 stata_pred_q50 stata_pred_q90 r_q10 r_q50 r_q90
export delimited "/tmp/grf_stata/tests/fidelity_reports/05_quantile/test16_stata.csv", replace
