*! grf_best_linear_projection.ado -- Best Linear Projection of CATE
*! Version 0.2.0
*! Projects AIPW DR scores onto covariates via OLS with HC3 robust SEs
*! Matches grf::best_linear_projection() which uses DR scores, not raw tau_hat
*!
*! NOTE: This is an eclass command. Running it replaces the causal forest
*! e() results with regression output (same as R's best_linear_projection).
*! Run grf_ate/grf_test_calibration BEFORE this command if needed.

program define grf_best_linear_projection, eclass
    version 14.0

    syntax [varlist(numeric default=none)] [if] [in]

    /* ---- Verify causal forest results exist ---- */
    if "`e(cmd)'" != "grf_causal_forest" {
        display as error "grf_best_linear_projection requires prior estimation by grf_causal_forest"
        exit 301
    }

    local depvar    "`e(depvar)'"
    local treatvar  "`e(treatvar)'"
    local tauvar    "`e(predict_var)'"
    local yhatvar   "`e(yhat_var)'"
    local whatvar   "`e(what_var)'"
    local indepvars_default "`e(indepvars)'"

    /* ---- Determine projection variables ---- */
    if "`varlist'" != "" {
        local proj_vars "`varlist'"
    }
    else {
        local proj_vars "`indepvars_default'"
    }

    if "`proj_vars'" == "" {
        display as error "no covariates specified and e(indepvars) is empty"
        exit 198
    }

    /* ---- Confirm required variables exist ---- */
    foreach v in depvar treatvar tauvar yhatvar whatvar {
        confirm numeric variable ``v''
    }
    foreach v of local proj_vars {
        confirm numeric variable `v'
    }

    /* ---- Mark sample ---- */
    marksample touse
    markout `touse' `depvar' `treatvar' `tauvar' `yhatvar' `whatvar' `proj_vars'
    quietly count if `touse'
    local n_use = r(N)

    if `n_use' < 2 {
        display as error "need at least 2 non-missing observations"
        exit 2000
    }

    /* ---- Compute AIPW doubly-robust scores ----
     * Same formula as grf_ate.ado (from grf's get_scores.R):
     *   DR_score_i = tau_hat_i
     *     + (W_i - W_hat_i) / Var(W - W_hat) * (Y_i - Y_hat_i - tau_hat_i * (W_i - W_hat_i))
     *
     * BLP then regresses DR_score on covariates with HC3 robust SEs.
     */
    tempvar w_resid y_resid dr_score
    quietly gen double `w_resid' = `treatvar' - `whatvar' if `touse'

    quietly summarize `w_resid' if `touse'
    local w_resid_var = r(Var)

    if `w_resid_var' < 1e-12 {
        display as error "variance of treatment residuals is near zero; cannot compute BLP"
        exit 498
    }

    quietly {
        gen double `y_resid'  = `depvar' - `yhatvar' - `tauvar' * `w_resid' if `touse'
        gen double `dr_score' = `tauvar' + (`w_resid' / `w_resid_var') * `y_resid' if `touse'
    }

    /* ---- Display header ---- */
    display as text ""
    display as text "Best Linear Projection of CATE (AIPW)"
    display as text "{hline 55}"
    display as text "CATE variable:         " as result "`tauvar'"
    display as text "Projection covariates: " as result "`proj_vars'"
    display as text "Observations:          " as result `n_use'
    display as text "{hline 55}"
    display as text ""
    display as text "Regression: DR_score ~ covariates, vce(hc3)"
    display as text ""

    /* ---- Run OLS: DR_score ~ covariates with HC3 robust SEs ---- */
    regress `dr_score' `proj_vars' if `touse', vce(hc3)

end
