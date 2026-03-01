*! grf_average_partial_effect.ado -- Average Partial Effect from causal forest
*! Version 0.2.0
*! Computes doubly-robust APE from grf_causal_forest ereturn results
*! For continuous treatment, APE estimates the average marginal effect of W on Y
*! Uses AIPW DR scores with optional debiasing weights
*!
*! Note: R's grf has deprecated average_partial_effect() in favor of
*! average_treatment_effect() with continuous treatment. Stata retains
*! grf_average_partial_effect as a standalone command for backward compatibility.
*! For new work, users may use grf_ate with continuous treatment.

program define grf_average_partial_effect, rclass
    version 14.0

    syntax [if] [in] [, DEBIASINGweights(varname numeric) ///
        NUMTreesvariance(integer 500) noCALIBRATE]

    /* ---- Verify causal forest results exist ---- */
    if "`e(cmd)'" != "grf_causal_forest" {
        display as error ///
            "grf_average_partial_effect requires prior estimation by grf_causal_forest"
        exit 301
    }

    local depvar    "`e(depvar)'"
    local treatvar  "`e(treatvar)'"
    local tauvar    "`e(predict_var)'"
    local yhatvar   "`e(yhat_var)'"
    local whatvar   "`e(what_var)'"
    local indepvars "`e(indepvars)'"

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
     *   + (W_i - W_hat_i) / Var(W - W_hat)
     *     * (Y_i - Y_hat_i - tau_hat_i * (W_i - W_hat_i))
     *
     * APE = weighted mean of DR scores (with debiasing weights)
     */

    /* Compute treatment residuals and their variance */
    tempvar w_resid y_resid dr_score
    quietly gen double `w_resid' = `treatvar' - `whatvar' if `touse'

    quietly summarize `w_resid' if `touse'
    local w_resid_var = r(Var)

    if `w_resid_var' < 1e-12 {
        display as error ///
            "variance of treatment residuals is near zero; cannot compute APE"
        exit 498
    }

    /* Generate DR scores */
    quietly {
        gen double `y_resid'  = `depvar' - `yhatvar' ///
            - `tauvar' * `w_resid' if `touse'
        gen double `dr_score' = `tauvar' ///
            + (`w_resid' / `w_resid_var') * `y_resid' if `touse'
    }

    /* ---- Compute debiasing weights ---- */
    tempvar dbweight
    if "`debiasweights'" != "" {
        /* User-supplied debiasing weights */
        confirm numeric variable `debiasweights'
        quietly gen double `dbweight' = `debiasweights' if `touse'
    }
    else {
        /* Default: estimate Var[W|X] from squared treatment residuals,
         * then use 1/Var[W|X] as debiasing weights.
         * For binary W, w_resid^2 approximates e(X)*(1-e(X)).
         * For continuous W, w_resid^2 estimates Var[W|X] at each obs. */
        quietly gen double `dbweight' = 1 / max(`w_resid'^2, 1e-6) if `touse'
    }

    /* Drop observations with non-positive weights */
    quietly replace `touse' = 0 if `dbweight' <= 0 & `touse'
    quietly count if `touse'
    local n_use = r(N)

    if `n_use' < 2 {
        display as error "fewer than 2 observations with positive debiasing weights"
        exit 2000
    }

    /* ---- Calibrate weights ---- */
    if "`calibrate'" != "nocalibrate" {
        /* Normalize weights so they sum to n */
        quietly summarize `dbweight' if `touse'
        local wt_sum = r(sum)
        if `wt_sum' < 1e-12 {
            display as error "sum of debiasing weights is zero"
            exit 498
        }
        quietly replace `dbweight' = `dbweight' * `n_use' / `wt_sum' if `touse'
    }

    /* ---- Compute weighted APE and SE ---- */
    /* APE = sum(w_i * DR_i) / sum(w_i) */
    quietly summarize `dbweight' if `touse'
    local sum_wt = r(sum)

    if `sum_wt' < 1e-12 {
        display as error "sum of weights is zero"
        exit 498
    }

    tempvar wt_dr
    quietly gen double `wt_dr' = `dbweight' * `dr_score' if `touse'
    quietly summarize `wt_dr' if `touse'
    local ape = r(sum) / `sum_wt'

    /* SE: sqrt( sum(w_i^2 * (DR_i - APE)^2) ) / sum(w_i) */
    tempvar wt_dev_sq
    quietly gen double `wt_dev_sq' = ///
        `dbweight'^2 * (`dr_score' - `ape')^2 if `touse'
    quietly summarize `wt_dev_sq' if `touse'
    local se = sqrt(r(sum)) / `sum_wt'

    /* 95% CI and p-value (normal approximation) */
    local ci_lower = `ape' - invnormal(0.975) * `se'
    local ci_upper = `ape' + invnormal(0.975) * `se'
    local tstat    = `ape' / `se'
    local pvalue   = 2 * (1 - normal(abs(`tstat')))

    /* ---- Display results ---- */
    display as text ""
    display as text "Average Partial Effect (AIPW)"
    display as text "{hline 60}"
    display as text "Outcome variable:      " as result "`depvar'"
    display as text "Treatment variable:    " as result "`treatvar'"
    if "`debiasweights'" != "" {
        display as text "Debiasing weights:     " as result "`debiasweights'"
    }
    else {
        display as text "Debiasing weights:     " as result "(estimated)"
    }
    if "`calibrate'" == "nocalibrate" {
        display as text "Calibrated:            " as result "no"
    }
    else {
        display as text "Calibrated:            " as result "yes"
    }
    display as text "Observations:          " as result `n_use'
    display as text "{hline 60}"
    display as text ""
    display as text %~14s "Estimate" %~14s "Std. Err." ///
        %~14s "[95% Conf." %~14s "Interval]"
    display as text "{hline 60}"
    display as text "APE" _col(10) ///
        as result %12.6f `ape' ///
        _col(24) as result %12.6f `se' ///
        _col(38) as result %12.6f `ci_lower' ///
        _col(52) as result %12.6f `ci_upper'
    display as text "{hline 60}"
    display as text "t-statistic:   " as result %8.4f `tstat'
    display as text "p-value:       " as result %8.4f `pvalue'
    display as text ""

    /* ---- Store results ---- */
    return scalar estimate  = `ape'
    return scalar std_err   = `se'
    return scalar ci_lower  = `ci_lower'
    return scalar ci_upper  = `ci_upper'
    return scalar pvalue    = `pvalue'
    return scalar N         = `n_use'
end
