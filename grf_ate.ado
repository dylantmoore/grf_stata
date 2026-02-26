*! grf_ate.ado -- Average Treatment Effect (AIPW) from causal forest
*! Version 0.1.0
*! Computes doubly-robust ATE from grf_causal_forest ereturn results

program define grf_ate, rclass
    version 14.0

    syntax [if] [in]

    /* ---- Verify causal forest results exist ---- */
    if "`e(cmd)'" != "grf_causal_forest" {
        display as error "grf_ate requires prior estimation by grf_causal_forest"
        exit 301
    }

    local depvar    "`e(depvar)'"
    local treatvar  "`e(treatvar)'"
    local tauvar    "`e(predict_var)'"
    local yhatvar   "`e(yhat_var)'"
    local whatvar   "`e(what_var)'"

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

    /* ---- Compute AIPW doubly-robust scores ----
     *
     * DR_score_i = tau_hat_i
     *   + (W_i - W_hat_i) / var(W_hat) * (Y_i - Y_hat_i - tau_hat_i * (W_i - W_hat_i))
     *
     * ATE = mean(DR_score_i)
     * SE  = sd(DR_score_i) / sqrt(n)
     */

    /* Compute var(W - W_hat) = variance of treatment residuals */
    tempvar w_resid y_resid dr_score
    quietly gen double `w_resid' = `treatvar' - `whatvar' if `touse'

    quietly summarize `w_resid' if `touse'
    local w_resid_var = r(Var)

    if `w_resid_var' < 1e-12 {
        display as error "variance of treatment residuals is near zero; cannot compute ATE"
        exit 498
    }

    /* Generate DR scores: tau_hat + (W-W_hat)/Var(W-W_hat) * (Y-Y_hat-tau_hat*(W-W_hat)) */
    quietly {
        gen double `y_resid'  = `depvar' - `yhatvar' - `tauvar' * `w_resid' if `touse'
        gen double `dr_score' = `tauvar' + (`w_resid' / `w_resid_var') * `y_resid' if `touse'
    }

    /* Compute ATE and SE */
    quietly summarize `dr_score' if `touse'
    local ate    = r(mean)
    local dr_sd  = r(sd)
    local se     = `dr_sd' / sqrt(`n_use')

    /* 95% CI and p-value (normal approximation) */
    local ci_lower = `ate' - invnormal(0.975) * `se'
    local ci_upper = `ate' + invnormal(0.975) * `se'
    local tstat    = `ate' / `se'
    local pvalue   = 2 * (1 - normal(abs(`tstat')))

    /* ---- Display results ---- */
    display as text ""
    display as text "Average Treatment Effect (AIPW)"
    display as text "{hline 60}"
    display as text "Outcome variable:      " as result "`depvar'"
    display as text "Treatment variable:    " as result "`treatvar'"
    display as text "Observations:          " as result `n_use'
    display as text "{hline 60}"
    display as text ""
    display as text %~14s "Estimate" %~14s "Std. Err." ///
        %~14s "[95% Conf." %~14s "Interval]"
    display as text "{hline 60}"
    display as text "ATE" _col(10) ///
        as result %12.6f `ate' ///
        _col(24) as result %12.6f `se' ///
        _col(38) as result %12.6f `ci_lower' ///
        _col(52) as result %12.6f `ci_upper'
    display as text "{hline 60}"
    display as text "t-statistic:   " as result %8.4f `tstat'
    display as text "p-value:       " as result %8.4f `pvalue'
    display as text ""

    /* ---- Store results ---- */
    return scalar ate      = `ate'
    return scalar se       = `se'
    return scalar ci_lower = `ci_lower'
    return scalar ci_upper = `ci_upper'
    return scalar pvalue   = `pvalue'
    return scalar N        = `n_use'
end
