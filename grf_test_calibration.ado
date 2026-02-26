*! grf_test_calibration.ado -- Calibration test for causal forest
*! Version 0.1.0
*! Implements Chernozhukov, Demirer, Duflo, Fernandez-Val (2020) calibration test

program define grf_test_calibration, rclass
    version 14.0

    syntax [if] [in]

    /* ---- Verify causal forest results exist ---- */
    if "`e(cmd)'" != "grf_causal_forest" {
        display as error "grf_test_calibration requires prior estimation by grf_causal_forest"
        exit 301
    }

    local depvar   "`e(depvar)'"
    local treatvar "`e(treatvar)'"
    local tauvar   "`e(predict_var)'"
    local yhatvar  "`e(yhat_var)'"
    local whatvar  "`e(what_var)'"

    /* ---- Confirm required variables exist ---- */
    foreach v in depvar treatvar tauvar yhatvar whatvar {
        confirm numeric variable ``v''
    }

    /* ---- Mark sample ---- */
    marksample touse
    markout `touse' `depvar' `treatvar' `tauvar' `yhatvar' `whatvar'
    quietly count if `touse'
    local n_use = r(N)

    if `n_use' < 2 {
        display as error "need at least 2 non-missing observations"
        exit 2000
    }

    /* ---- Compute calibration test ----
     *
     * Following grf's test_calibration:
     *   1. Compute residuals: Y_resid = Y - Y_hat, W_resid = W - W_hat
     *   2. Compute mean of tau_hat
     *   3. Construct:
     *        LHS = Y_resid
     *        RHS_1 = W_resid * 1           (mean forest prediction)
     *        RHS_2 = W_resid * (tau_hat - mean(tau_hat))  (differential prediction)
     *   4. Regress LHS on RHS_1 RHS_2 (no constant)
     *
     * Interpretation:
     *   - Coefficient on RHS_1: should be close to the ATE if forest is correct on average
     *   - Coefficient on RHS_2: should be close to 1 if heterogeneity is well-calibrated
     */

    /* Compute mean of tau_hat */
    quietly summarize `tauvar' if `touse'
    local tau_mean = r(mean)

    /* Generate regression variables */
    tempvar y_resid w_resid cal_lhs cal_rhs1 cal_rhs2

    quietly {
        gen double `y_resid'  = `depvar' - `yhatvar' if `touse'
        gen double `w_resid'  = `treatvar' - `whatvar' if `touse'
        gen double `cal_lhs'  = `y_resid' if `touse'
        gen double `cal_rhs1' = `w_resid' if `touse'
        gen double `cal_rhs2' = `w_resid' * (`tauvar' - `tau_mean') if `touse'
    }

    /* ---- Display header ---- */
    display as text ""
    display as text "Calibration Test for Causal Forest"
    display as text "(Chernozhukov, Demirer, Duflo, Fernandez-Val 2020)"
    display as text "{hline 65}"
    display as text "Outcome variable:      " as result "`depvar'"
    display as text "Treatment variable:    " as result "`treatvar'"
    display as text "CATE variable:         " as result "`tauvar'"
    display as text "Observations:          " as result `n_use'
    display as text "Mean tau_hat:          " as result %9.6f `tau_mean'
    display as text "{hline 65}"
    display as text ""
    display as text "Regression: (Y - Y_hat) on (W - W_hat)*{1, tau_hat - mean(tau_hat)}"
    display as text "  noconstant"
    display as text ""

    /* ---- Run calibration regression (no constant, HC3 robust SEs) ----
     * NOTE: regress overwrites e() results from grf_causal_forest.
     * We save the needed causal forest e() macros into locals before regress,
     * so downstream post-estimation commands that depend on e() will still
     * need to re-fit the causal forest (or call this last).
     */
    quietly regress `cal_lhs' `cal_rhs1' `cal_rhs2' if `touse', noconstant vce(hc3)

    /* Extract coefficients */
    local b_mean  = _b[`cal_rhs1']
    local b_diff  = _b[`cal_rhs2']
    local se_mean = _se[`cal_rhs1']
    local se_diff = _se[`cal_rhs2']
    local t_mean  = `b_mean' / `se_mean'
    local t_diff  = `b_diff' / `se_diff'
    local p_mean  = 2 * (1 - normal(abs(`t_mean')))
    local p_diff  = 2 * (1 - normal(abs(`t_diff')))
    local n_reg   = e(N)

    /* ---- Display results table ---- */
    display as text ""
    display as text %~25s " " %~12s "Coef." %~12s "Std. Err." ///
        %~10s "t" %~10s "P>|t|"
    display as text "{hline 65}"
    display as text "mean.forest.prediction" _col(26) ///
        as result %10.6f `b_mean' ///
        _col(39) as result %10.6f `se_mean' ///
        _col(51) as result %8.3f `t_mean' ///
        _col(61) as result %6.4f `p_mean'
    display as text "differential.forest.prediction" _col(26) ///
        as result %10.6f `b_diff' ///
        _col(39) as result %10.6f `se_diff' ///
        _col(51) as result %8.3f `t_diff' ///
        _col(61) as result %6.4f `p_diff'
    display as text "{hline 65}"
    display as text ""
    display as text "Note: A well-calibrated forest has differential.forest.prediction"
    display as text "      coefficient close to 1."
    display as text ""

    /* ---- Store results ---- */
    return scalar b_mean   = `b_mean'
    return scalar se_mean  = `se_mean'
    return scalar t_mean   = `t_mean'
    return scalar p_mean   = `p_mean'
    return scalar b_diff   = `b_diff'
    return scalar se_diff  = `se_diff'
    return scalar t_diff   = `t_diff'
    return scalar p_diff   = `p_diff'
    return scalar N        = `n_reg'
end
