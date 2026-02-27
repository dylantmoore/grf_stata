*! grf_get_scores.ado -- Extract doubly-robust scores from causal/instrumental forest
*! Version 0.1.0
*! Computes AIPW doubly-robust scores Gamma_i for post-estimation analysis

program define grf_get_scores, rclass
    version 14.0

    syntax , GENerate(name) [REPlace]

    /* ---- Verify prior estimation ---- */
    local forest_type "`e(forest_type)'"
    if "`forest_type'" != "causal" & "`forest_type'" != "instrumental" {
        display as error "grf_get_scores requires prior estimation by" ///
            " grf_causal_forest or grf_instrumental_forest"
        exit 301
    }

    /* ---- Handle replace ---- */
    if "`replace'" != "" {
        capture drop `generate'
    }
    confirm new variable `generate'

    /* ---- Read e() results ---- */
    local depvar   "`e(depvar)'"
    local tauvar   "`e(predict_var)'"

    if "`forest_type'" == "causal" {
        local treatvar "`e(treatvar)'"
        local yhatvar  "`e(yhat_var)'"
        local whatvar  "`e(what_var)'"
    }
    else if "`forest_type'" == "instrumental" {
        local treatvar "`e(treatment)'"
        /* Instrumental forest stores nuisance vars differently;
         * the centering is done internally. We need Y, W, Y.hat, W.hat.
         * The instrumental forest stores predict_var as LATE estimates.
         * For DR scores, we use the same AIPW formula with the LATE as tau_hat. */
        local yhatvar  "`e(yhat_var)'"
        local whatvar  "`e(what_var)'"

        /* Check if nuisance variables are available */
        if "`yhatvar'" == "" | "`whatvar'" == "" {
            /* Instrumental forest may not store nuisance vars externally.
             * Check for internal _grf_ variables */
            capture confirm numeric variable _grf_yhat
            if !_rc {
                local yhatvar "_grf_yhat"
            }
            capture confirm numeric variable _grf_what
            if !_rc {
                local whatvar "_grf_what"
            }
        }
    }

    /* ---- Confirm required variables exist ---- */
    foreach v in depvar treatvar tauvar yhatvar whatvar {
        if "``v''" == "" {
            display as error "required variable not found: `v'"
            display as error "the prior estimation must store Y.hat and W.hat"
            display as error "(for causal forest, use the default settings)"
            exit 111
        }
        confirm numeric variable ``v''
    }

    /* ---- Mark sample ---- */
    tempvar touse
    quietly gen byte `touse' = 1
    markout `touse' `depvar' `treatvar' `tauvar' `yhatvar' `whatvar'
    quietly count if `touse'
    local n_use = r(N)

    if `n_use' < 2 {
        display as error "need at least 2 non-missing observations"
        exit 2000
    }

    /* ---- Compute doubly-robust scores ----
     *
     * Gamma_i = tau_hat_i
     *   + (W_i - W_hat_i) / Var(W - W_hat)
     *     * (Y_i - Y_hat_i - tau_hat_i * (W_i - W_hat_i))
     *
     * This is the AIPW doubly-robust score. Under correct specification
     * of either the outcome model or the propensity model, the mean of
     * these scores is a consistent estimator of the ATE.
     */

    /* Compute treatment residuals and their variance */
    tempvar w_resid y_resid
    quietly gen double `w_resid' = `treatvar' - `whatvar' if `touse'

    quietly summarize `w_resid' if `touse'
    local w_resid_var = r(Var)

    if `w_resid_var' < 1e-12 {
        display as error "variance of treatment residuals (W - W.hat) is near zero"
        display as error "cannot compute doubly-robust scores"
        exit 498
    }

    /* Compute outcome residual: Y - Y.hat - tau_hat * (W - W.hat) */
    quietly gen double `y_resid' = `depvar' - `yhatvar' ///
        - `tauvar' * `w_resid' if `touse'

    /* Compute DR scores */
    quietly gen double `generate' = `tauvar' ///
        + (`w_resid' / `w_resid_var') * `y_resid' if `touse'

    label variable `generate' "Doubly-robust scores from `forest_type' forest"

    /* ---- Summary statistics ---- */
    quietly summarize `generate' if `touse'
    local n_scores = r(N)
    local score_mean = r(mean)
    local score_sd = r(sd)
    local score_min = r(min)
    local score_max = r(max)
    local score_se = `score_sd' / sqrt(`n_scores')

    /* ---- Display results ---- */
    display as text ""
    display as text "Doubly-Robust Scores"
    display as text "{hline 55}"
    display as text "Forest type:           " as result "`forest_type'"
    display as text "Outcome variable:      " as result "`depvar'"
    display as text "Treatment variable:    " as result "`treatvar'"
    display as text "CATE variable:         " as result "`tauvar'"
    display as text "Observations:          " as result `n_scores'
    display as text "{hline 55}"
    display as text ""
    display as text "Summary of DR scores (`generate'):"
    display as text "  Mean:         " as result %12.6f `score_mean'
    display as text "  Std. Dev.:    " as result %12.6f `score_sd'
    display as text "  Std. Err.:    " as result %12.6f `score_se'
    display as text "  Min:          " as result %12.6f `score_min'
    display as text "  Max:          " as result %12.6f `score_max'
    display as text ""
    display as text "Note: mean(DR scores) = AIPW estimate of ATE = " ///
        as result %9.6f `score_mean'
    display as text ""

    /* ---- Store results ---- */
    return scalar N        = `n_scores'
    return scalar mean     = `score_mean'
    return scalar sd       = `score_sd'
    return scalar se       = `score_se'
    return scalar min      = `score_min'
    return scalar max      = `score_max'
    return local  generate   "`generate'"
    return local  forest_type "`forest_type'"
end
