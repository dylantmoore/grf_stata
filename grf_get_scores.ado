*! grf_get_scores.ado -- Extract DR/IPCW scores from supported GRF causal forests
*! Version 0.3.0
*! Computes AIPW-style or IPCW-moment-corrected scores for post-estimation analysis

program define grf_get_scores, rclass
    version 14.0

    syntax , GENerate(name) [REPlace]

    /* ---- Verify prior estimation ---- */
    local forest_type "`e(forest_type)'"
    local cluster_var "`e(cluster_var)'"
    if "`forest_type'" != "causal" & "`forest_type'" != "instrumental" ///
        & "`forest_type'" != "multi_causal" & "`forest_type'" != "causal_survival" {
        display as error "grf_get_scores requires prior estimation by" ///
            " grf_causal_forest, grf_instrumental_forest," ///
            " grf_multi_arm_causal_forest, or grf_causal_survival_forest"
        exit 301
    }

    /* ====================================================================
     * Multi-arm causal forest: per-arm DR scores
     * ==================================================================== */
    if "`forest_type'" == "multi_causal" {
        local n_treat = e(n_treat)
        local depvar   "`e(depvar)'"
        local yhatvar  "`e(yhat_var)'"
        local predict_stub "`e(predict_stub)'"

        /* Handle replace for all arm-specific score variables */
        forvalues j = 1/`n_treat' {
            if "`replace'" != "" {
                capture drop `generate'_t`j'
            }
            confirm new variable `generate'_t`j'
        }

        /* Mark sample */
        tempvar touse
        quietly gen byte `touse' = 1
        markout `touse' `depvar' `yhatvar'
        forvalues j = 1/`n_treat' {
            local tau_j "`predict_stub'_t`j'"
            local what_j "`e(what_var_`j')'"
            markout `touse' `tau_j' `what_j'
        }
        quietly count if `touse'
        local n_use = r(N)

        if `n_use' < 2 {
            display as error "need at least 2 non-missing observations"
            exit 2000
        }

        /* Compute per-arm DR scores (joint Y residual, per R's get_scores):
         *
         * Step 1: For each arm j, compute w_resid_j = W_j - W_hat_j
         * Step 2: Compute joint Y residual:
         *   Y_resid = Y - Y_hat - SUM_k(tau_hat_k * w_resid_k)
         * Step 3: For each arm j:
         *   DR_j = tau_hat_j + w_resid_j / Var(w_resid_j) * Y_resid
         */
        local output_vars ""
        local treatvars "`e(treatvars)'"

        /* Step 1: Compute all treatment residuals and their variances */
        forvalues j = 1/`n_treat' {
            local what_j "`e(what_var_`j')'"
            local tv : word `j' of `treatvars'
            tempvar w_resid_`j'
            quietly gen double `w_resid_`j'' = `tv' - `what_j' if `touse'
            quietly summarize `w_resid_`j'' if `touse'
            local w_resid_var_`j' = r(Var)
            if `w_resid_var_`j'' < 1e-12 {
                local w_resid_var_`j' = 1e-12
            }
        }

        /* Step 2: Compute joint Y residual subtracting ALL arms' effects */
        tempvar y_resid_joint
        quietly gen double `y_resid_joint' = `depvar' - `yhatvar' if `touse'
        forvalues k = 1/`n_treat' {
            local tau_k "`predict_stub'_t`k'"
            quietly replace `y_resid_joint' = `y_resid_joint' ///
                - `tau_k' * `w_resid_`k'' if `touse'
        }

        /* Step 3: Compute per-arm DR scores using joint Y residual */
        forvalues j = 1/`n_treat' {
            local tau_j "`predict_stub'_t`j'"
            quietly gen double `generate'_t`j' = `tau_j' ///
                + (`w_resid_`j'' / `w_resid_var_`j'') * `y_resid_joint' if `touse'
            label variable `generate'_t`j' "DR scores arm `j' from multi_causal forest"
            local output_vars "`output_vars' `generate'_t`j'"
        }

        /* Display results */
        display as text ""
        display as text "Multi-Arm Doubly-Robust Scores"
        display as text "{hline 55}"
        display as text "Forest type:           " as result "multi_causal"
        display as text "Outcome variable:      " as result "`depvar'"
        display as text "Treatment arms:        " as result `n_treat'
        display as text "Observations:          " as result `n_use'
        display as text "{hline 55}"

        forvalues j = 1/`n_treat' {
            quietly summarize `generate'_t`j' if `touse'
            display as text ""
            display as text "Arm `j' DR scores (`generate'_t`j'):"
            display as text "  Mean (ATE_`j'): " as result %12.6f r(mean)
            display as text "  Std. Dev.:      " as result %12.6f r(sd)
        }
        display as text ""

        /* Store results */
        return scalar N         = `n_use'
        return scalar n_treat   = `n_treat'
        forvalues j = 1/`n_treat' {
            quietly summarize `generate'_t`j' if `touse'
            return scalar mean_`j' = r(mean)
            return scalar sd_`j'   = r(sd)
            return scalar se_`j'   = r(sd) / sqrt(r(N))
        }
        return local  generate    "`output_vars'"
        return local  forest_type "multi_causal"
        exit
    }

    /* ====================================================================
     * Causal survival forest: IPCW/DR scores from nuisance moments
     * ==================================================================== */
    if "`forest_type'" == "causal_survival" {
        local tauvar      "`e(predict_var)'"
        local numervar    "`e(numer_var)'"
        local denomvar    "`e(denom_var)'"
        local nuisance_mode "`e(nuisance_mode)'"
        local cluster_var "`e(cluster_var)'"

        /* Handle replace */
        if "`replace'" != "" {
            capture drop `generate'
        }
        confirm new variable `generate'

        /* Fall back to canonical nuisance variable names if e() macros are absent */
        if "`numervar'" == "" {
            capture confirm numeric variable _grf_cs_numer
            if !_rc local numervar "_grf_cs_numer"
        }
        if "`denomvar'" == "" {
            capture confirm numeric variable _grf_cs_denom
            if !_rc local denomvar "_grf_cs_denom"
        }

        /* Confirm required variables */
        foreach v in tauvar numervar denomvar {
            if "``v''" == "" {
                display as error "required variable not found: `v'"
                exit 111
            }
            confirm numeric variable ``v''
        }

        /* Mark sample */
        tempvar touse
        quietly gen byte `touse' = 1
        markout `touse' `tauvar' `numervar' `denomvar'
        quietly count if `touse'
        local n_use = r(N)

        if `n_use' < 2 {
            display as error "need at least 2 non-missing observations"
            exit 2000
        }

        /* IPCW/DR moment correction for ratio estimands:
         *
         *   tau_hat_i = forest prediction
         *   numer_i, denom_i = causal-survival nuisance moments
         *
         * Influence-style score:
         *   score_i = tau_hat_i + (numer_i - denom_i * tau_hat_i) / E[denom]
         *
         * This uses the same orthogonal moments that define the causal-survival
         * forest objective, and reduces to plug-in tau only when moment residuals
         * are identically zero.
         */
        quietly summarize `denomvar' if `touse', meanonly
        local denom_mean = r(mean)
        if abs(`denom_mean') < 1e-12 {
            display as error "mean IPCW denominator is near zero; cannot compute scores"
            exit 498
        }

        tempvar moment_resid
        quietly gen double `moment_resid' = `numervar' - `denomvar' * `tauvar' if `touse'
        quietly gen double `generate' = `tauvar' + (`moment_resid' / `denom_mean') if `touse'
        label variable `generate' "IPCW doubly-robust scores from causal_survival forest"

        /* Summary statistics */
        quietly summarize `generate' if `touse', detail
        local n_scores = r(N)
        local score_mean = r(mean)
        local score_sd = r(sd)
        local score_min = r(min)
        local score_max = r(max)

        if "`cluster_var'" != "" {
            confirm numeric variable `cluster_var'
            tempvar score_dev cl_sum cl_tag cl_sum_sq
            quietly gen double `score_dev' = `generate' - `score_mean' if `touse'
            quietly bysort `cluster_var': egen double `cl_sum' = total(`score_dev') if `touse'
            quietly bysort `cluster_var': gen byte `cl_tag' = (_n == 1) if `touse'
            quietly gen double `cl_sum_sq' = `cl_sum'^2 if `cl_tag' & `touse'
            quietly summarize `cl_sum_sq' if `cl_tag' & `touse', meanonly
            local score_se = sqrt(r(sum)) / `n_scores'
        }
        else {
            local score_se = `score_sd' / sqrt(`n_scores')
        }

        /* Display results */
        display as text ""
        display as text "Causal Survival IPCW Doubly-Robust Scores"
        display as text "{hline 55}"
        display as text "Forest type:           " as result "causal_survival"
        display as text "Nuisance mode:         " as result cond("`nuisance_mode'"=="", "unspecified", "`nuisance_mode'")
        display as text "Observations:          " as result `n_scores'
        display as text "E[denominator]:        " as result %12.6f `denom_mean'
        display as text "{hline 55}"
        display as text ""
        display as text "Summary of IPCW DR scores (`generate'):"
        display as text "  Mean:         " as result %12.6f `score_mean'
        display as text "  Std. Dev.:    " as result %12.6f `score_sd'
        display as text "  Std. Err.:    " as result %12.6f `score_se'
        display as text "  Min:          " as result %12.6f `score_min'
        display as text "  Max:          " as result %12.6f `score_max'
        display as text ""
        display as text "Note: scores use stored IPCW nuisance moments (numerator/denominator)."
        display as text ""

        /* Store results */
        return scalar N        = `n_scores'
        return scalar mean     = `score_mean'
        return scalar sd       = `score_sd'
        return scalar se       = `score_se'
        return scalar min      = `score_min'
        return scalar max      = `score_max'
        return scalar denom_mean = `denom_mean'
        return local  generate   "`generate'"
        return local  forest_type "causal_survival"
        exit
    }

    /* ====================================================================
     * Causal and instrumental forest: standard DR scores
     * ==================================================================== */

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
        local yhatvar  "`e(yhat_var)'"
        local whatvar  "`e(what_var)'"

        /* Check if nuisance variables are available */
        if "`yhatvar'" == "" | "`whatvar'" == "" {
            capture confirm numeric variable _grf_if_yhat
            if !_rc {
                local yhatvar "_grf_if_yhat"
            }
            else {
                capture confirm numeric variable _grf_yhat
                if !_rc {
                    local yhatvar "_grf_yhat"
                }
            }
            capture confirm numeric variable _grf_if_what
            if !_rc {
                local whatvar "_grf_if_what"
            }
            else {
                capture confirm numeric variable _grf_what
                if !_rc {
                    local whatvar "_grf_what"
                }
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

    if "`cluster_var'" != "" {
        /* Cluster-robust SE of the mean */
        confirm numeric variable `cluster_var'
        tempvar score_dev cl_sum cl_tag cl_sum_sq
        quietly gen double `score_dev' = `generate' - `score_mean' if `touse'
        quietly bysort `cluster_var': egen double `cl_sum' = total(`score_dev') if `touse'
        quietly bysort `cluster_var': gen byte `cl_tag' = (_n == 1) if `touse'
        quietly gen double `cl_sum_sq' = `cl_sum'^2 if `cl_tag' & `touse'
        quietly summarize `cl_sum_sq' if `cl_tag' & `touse'
        local score_se = sqrt(r(sum)) / `n_scores'
    }
    else {
        local score_se = `score_sd' / sqrt(`n_scores')
    }

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
