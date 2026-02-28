*! grf_ate.ado -- Average Treatment Effect (AIPW) from causal forest
*! Version 0.2.0
*! Computes doubly-robust ATE from grf_causal_forest ereturn results
*! Supports target.sample: all, treated, control, overlap

program define grf_ate, rclass
    version 14.0

    syntax [if] [in] [, TARGETsample(string) DEBIASINGweights(varname numeric) METHod(string)]

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
    local cluster_var "`e(cluster_var)'"

    /* ---- Default and validate method ---- */
    if "`method'" == "" {
        local method "AIPW"
    }
    local method = upper("`method'")
    if "`method'" != "AIPW" & "`method'" != "TMLE" {
        display as error "method() must be AIPW or TMLE"
        exit 198
    }

    /* ---- Default and validate target.sample ---- */
    if "`targetsample'" == "" {
        local targetsample "all"
    }
    if "`targetsample'" != "all" & "`targetsample'" != "treated" ///
        & "`targetsample'" != "control" & "`targetsample'" != "overlap" {
        display as error ///
            "target.sample must be one of: all, treated, control, overlap"
        exit 198
    }

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

    if "`method'" == "AIPW" {

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

    /* ---- Apply debiasing weights if specified ---- */
    if "`debiasingweights'" != "" {
        confirm numeric variable `debiasingweights'
        markout `touse' `debiasingweights'
        quietly count if `touse'
        local n_use = r(N)
        tempvar dw_dr
        quietly gen double `dw_dr' = `dr_score' * `debiasingweights' if `touse'
        quietly replace `dr_score' = `dw_dr' if `touse'
    }

    /* ---- Compute weights based on target.sample ---- */
    tempvar tsweight
    if "`targetsample'" == "all" {
        quietly gen double `tsweight' = 1 if `touse'
    }
    else if "`targetsample'" == "treated" {
        quietly gen double `tsweight' = `treatvar' if `touse'
    }
    else if "`targetsample'" == "control" {
        quietly gen double `tsweight' = (1 - `treatvar') if `touse'
    }
    else if "`targetsample'" == "overlap" {
        quietly gen double `tsweight' = `whatvar' * (1 - `whatvar') if `touse'
    }

    /* Check that weights have positive sum */
    quietly summarize `tsweight' if `touse'
    local sum_wt = r(sum)
    if `sum_wt' < 1e-12 {
        display as error "sum of target.sample weights is zero"
        exit 498
    }

    /* Compute weighted ATE and SE */
    tempvar wt_dr
    quietly gen double `wt_dr' = `tsweight' * `dr_score' if `touse'
    quietly summarize `wt_dr' if `touse'
    local ate = r(sum) / `sum_wt'

    /* SE: depends on whether clustering is active */
    if "`cluster_var'" != "" {
        /* Cluster-robust SE:
         * 1. Sum weighted deviations within clusters
         * 2. Square cluster-level sums
         * 3. SE = sqrt(sum of squared cluster sums) / sum(w)
         */
        confirm numeric variable `cluster_var'
        tempvar wt_dev cl_sum cl_tag cl_sum_sq
        quietly gen double `wt_dev' = `tsweight' * (`dr_score' - `ate') if `touse'
        quietly bysort `cluster_var': egen double `cl_sum' = total(`wt_dev') if `touse'
        quietly bysort `cluster_var': gen byte `cl_tag' = (_n == 1) if `touse'
        quietly gen double `cl_sum_sq' = `cl_sum'^2 if `cl_tag' & `touse'
        quietly summarize `cl_sum_sq' if `cl_tag' & `touse'
        local se = sqrt(r(sum)) / `sum_wt'
    }
    else {
        /* iid SE: sqrt( sum(w_i^2 * (Gamma_i - ATE)^2) ) / sum(w_i) */
        tempvar wt_dev_sq
        quietly gen double `wt_dev_sq' = ///
            `tsweight'^2 * (`dr_score' - `ate')^2 if `touse'
        quietly summarize `wt_dev_sq' if `touse'
        local se = sqrt(r(sum)) / `sum_wt'
    }

    }
    else if "`method'" == "TMLE" {

    /* ---- Compute TMLE (Targeted Maximum Likelihood Estimation) ----
     *
     * TMLE corrects the initial CATE estimates via targeted bias correction.
     * For target.sample=="overlap", TMLE is not applicable; fall back to AIPW.
     *
     * Note: All regressions are computed manually (OLS formula) to avoid
     * clobbering e() results from grf_causal_forest.
     */

    if "`targetsample'" == "overlap" {
        display as text "Note: TMLE not applicable for overlap weights; using AIPW."
        local method "AIPW"
    }

    if "`method'" == "AIPW" {
        /* Overlap fallback: replicate AIPW logic inline */
        tempvar w_resid y_resid dr_score
        quietly gen double `w_resid' = `treatvar' - `whatvar' if `touse'
        quietly summarize `w_resid' if `touse'
        local w_resid_var = r(Var)
        if `w_resid_var' < 1e-12 {
            display as error "variance of treatment residuals is near zero"
            exit 498
        }
        quietly {
            gen double `y_resid' = `depvar' - `yhatvar' - `tauvar' * `w_resid' if `touse'
            gen double `dr_score' = `tauvar' + (`w_resid' / `w_resid_var') * `y_resid' if `touse'
        }
        if "`debiasingweights'" != "" {
            confirm numeric variable `debiasingweights'
            markout `touse' `debiasingweights'
            quietly count if `touse'
            local n_use = r(N)
            tempvar dw_dr
            quietly gen double `dw_dr' = `dr_score' * `debiasingweights' if `touse'
            quietly replace `dr_score' = `dw_dr' if `touse'
        }
        tempvar tsweight
        quietly gen double `tsweight' = `whatvar' * (1 - `whatvar') if `touse'
        quietly summarize `tsweight' if `touse'
        local sum_wt = r(sum)
        if `sum_wt' < 1e-12 {
            display as error "sum of target.sample weights is zero"
            exit 498
        }
        tempvar wt_dr
        quietly gen double `wt_dr' = `tsweight' * `dr_score' if `touse'
        quietly summarize `wt_dr' if `touse'
        local ate = r(sum) / `sum_wt'
        if "`cluster_var'" != "" {
            confirm numeric variable `cluster_var'
            tempvar wt_dev cl_sum cl_tag cl_sum_sq
            quietly gen double `wt_dev' = `tsweight' * (`dr_score' - `ate') if `touse'
            quietly bysort `cluster_var': egen double `cl_sum' = total(`wt_dev') if `touse'
            quietly bysort `cluster_var': gen byte `cl_tag' = (_n == 1) if `touse'
            quietly gen double `cl_sum_sq' = `cl_sum'^2 if `cl_tag' & `touse'
            quietly summarize `cl_sum_sq' if `cl_tag' & `touse'
            local se = sqrt(r(sum)) / `sum_wt'
        }
        else {
            tempvar wt_dev_sq
            quietly gen double `wt_dev_sq' = ///
                `tsweight'^2 * (`dr_score' - `ate')^2 if `touse'
            quietly summarize `wt_dev_sq' if `touse'
            local se = sqrt(r(sum)) / `sum_wt'
        }
    }
    else {
        /* ---- True TMLE path ---- */

        /* Step 1: Initial potential outcome estimates */
        tempvar yhat0 yhat1
        quietly gen double `yhat0' = `yhatvar' - `tauvar' * `whatvar' if `touse'
        quietly gen double `yhat1' = `yhatvar' + `tauvar' * (1 - `whatvar') if `touse'

        /* Step 2: Average initial CATE */
        quietly summarize `tauvar' if `touse'
        local tau_avg = r(mean)

        if "`targetsample'" == "all" {
            /* ---- TMLE for ATE (all) ----
             * Controls (W=0): regress (Y - Y.hat.0) on 1/(1-W.hat), no intercept => epsilon0
             * Treated (W=1): regress (Y - Y.hat.1) on 1/W.hat, no intercept => epsilon1
             */
            tempvar resid0 clever0 resid1 clever1
            quietly gen double `resid0' = `depvar' - `yhat0' if `touse' & `treatvar' == 0
            quietly gen double `clever0' = 1 / (1 - `whatvar') if `touse' & `treatvar' == 0
            quietly gen double `resid1' = `depvar' - `yhat1' if `touse' & `treatvar' == 1
            quietly gen double `clever1' = 1 / `whatvar' if `touse' & `treatvar' == 1

            /* Manual OLS (no intercept): epsilon = sum(X*Y) / sum(X^2) */
            tempvar _xy0 _xx0 _xy1 _xx1
            quietly gen double `_xy0' = `clever0' * `resid0' if `touse' & `treatvar' == 0
            quietly gen double `_xx0' = `clever0'^2 if `touse' & `treatvar' == 0
            quietly summarize `_xy0' if `touse' & `treatvar' == 0
            local _sum_xy0 = r(sum)
            quietly summarize `_xx0' if `touse' & `treatvar' == 0
            local epsilon0 = `_sum_xy0' / r(sum)

            quietly gen double `_xy1' = `clever1' * `resid1' if `touse' & `treatvar' == 1
            quietly gen double `_xx1' = `clever1'^2 if `touse' & `treatvar' == 1
            quietly summarize `_xy1' if `touse' & `treatvar' == 1
            local _sum_xy1 = r(sum)
            quietly summarize `_xx1' if `touse' & `treatvar' == 1
            local epsilon1 = `_sum_xy1' / r(sum)

            /* Step 5: DR correction */
            quietly summarize `clever1' if `touse' & `treatvar' == 1
            local mean_clever1 = r(mean)
            quietly count if `touse' & `treatvar' == 1
            local n1 = r(N)

            quietly summarize `clever0' if `touse' & `treatvar' == 0
            local mean_clever0 = r(mean)
            quietly count if `touse' & `treatvar' == 0
            local n0 = r(N)

            local dr_correction = `epsilon1' * `mean_clever1' - `epsilon0' * `mean_clever0'

            /* Step 6: ATE */
            local ate = `tau_avg' + `dr_correction'

            /* Step 7: Sandwich SE */
            /* Influence function: psi_i for each obs */
            tempvar psi
            quietly gen double `psi' = . if `touse'
            /* Treated: psi = tau_hat - ate + clever1 * (Y - yhat1 - epsilon1 * clever1) ... simplified:
             * Use the standard TMLE influence function:
             * psi_i = (yhat1* - yhat0*) - ate + W*clever1*(Y-yhat1*) - (1-W)*clever0*(Y-yhat0*)
             * where yhat1* = yhat1 + epsilon1/W.hat, yhat0* = yhat0 + epsilon0/(1-W.hat)
             */
            tempvar yhat0_star yhat1_star
            quietly gen double `yhat0_star' = `yhat0' + `epsilon0' / (1 - `whatvar') if `touse'
            quietly gen double `yhat1_star' = `yhat1' + `epsilon1' / `whatvar' if `touse'
            quietly replace `psi' = (`yhat1_star' - `yhat0_star') - `ate' ///
                + `treatvar' / `whatvar' * (`depvar' - `yhat1_star') ///
                - (1 - `treatvar') / (1 - `whatvar') * (`depvar' - `yhat0_star') if `touse'

            if "`cluster_var'" != "" {
                confirm numeric variable `cluster_var'
                tempvar cl_psi cl_tag_t cl_psi_sq
                quietly bysort `cluster_var': egen double `cl_psi' = total(`psi') if `touse'
                quietly bysort `cluster_var': gen byte `cl_tag_t' = (_n == 1) if `touse'
                quietly gen double `cl_psi_sq' = `cl_psi'^2 if `cl_tag_t' & `touse'
                quietly summarize `cl_psi_sq' if `cl_tag_t' & `touse'
                local se = sqrt(r(sum)) / `n_use'
            }
            else {
                quietly summarize `psi' if `touse'
                local se = r(sd) / sqrt(`n_use')
            }
        }
        else if "`targetsample'" == "treated" {
            /* TMLE for ATT: fit control regression with weights W.hat/(1-W.hat) */
            tempvar resid0 clever0 ipw_wt
            quietly gen double `resid0' = `depvar' - `yhat0' if `touse' & `treatvar' == 0
            quietly gen double `clever0' = 1 / (1 - `whatvar') if `touse' & `treatvar' == 0
            quietly gen double `ipw_wt' = `whatvar' / (1 - `whatvar') if `touse' & `treatvar' == 0

            /* Weighted OLS (no intercept): epsilon = sum(w*X*Y) / sum(w*X^2) */
            tempvar _wxy0 _wxx0
            quietly gen double `_wxy0' = `ipw_wt' * `clever0' * `resid0' if `touse' & `treatvar' == 0
            quietly gen double `_wxx0' = `ipw_wt' * `clever0'^2 if `touse' & `treatvar' == 0
            quietly summarize `_wxy0' if `touse' & `treatvar' == 0
            local _sum_wxy0 = r(sum)
            quietly summarize `_wxx0' if `touse' & `treatvar' == 0
            local epsilon0 = `_sum_wxy0' / r(sum)

            /* Treated side: simple mean correction */
            tempvar resid1 clever1
            quietly gen double `resid1' = `depvar' - `yhat1' if `touse' & `treatvar' == 1
            quietly gen double `clever1' = 1 / `whatvar' if `touse' & `treatvar' == 1

            /* Manual OLS (no intercept) */
            tempvar _xy1 _xx1
            quietly gen double `_xy1' = `clever1' * `resid1' if `touse' & `treatvar' == 1
            quietly gen double `_xx1' = `clever1'^2 if `touse' & `treatvar' == 1
            quietly summarize `_xy1' if `touse' & `treatvar' == 1
            local _sum_xy1 = r(sum)
            quietly summarize `_xx1' if `touse' & `treatvar' == 1
            local epsilon1 = `_sum_xy1' / r(sum)

            /* Corrected potential outcomes */
            tempvar yhat0_star yhat1_star
            quietly gen double `yhat0_star' = `yhat0' + `epsilon0' / (1 - `whatvar') if `touse'
            quietly gen double `yhat1_star' = `yhat1' + `epsilon1' / `whatvar' if `touse'

            /* ATT = mean among treated of (yhat1* - yhat0*) + correction */
            quietly count if `touse' & `treatvar' == 1
            local n1 = r(N)

            tempvar psi
            quietly gen double `psi' = . if `touse'
            quietly replace `psi' = (`yhat1_star' - `yhat0_star') ///
                + `treatvar' / `whatvar' * (`depvar' - `yhat1_star') ///
                - (1 - `treatvar') * `whatvar' / (1 - `whatvar') / `whatvar' * (`depvar' - `yhat0_star') if `touse'

            /* Weight by treatment indicator for ATT */
            tempvar tsweight wt_psi
            quietly gen double `tsweight' = `treatvar' if `touse'
            quietly gen double `wt_psi' = `tsweight' * `psi' if `touse'
            quietly summarize `wt_psi' if `touse'
            local sum_wt = `n1'
            local ate = r(sum) / `sum_wt'

            if "`cluster_var'" != "" {
                confirm numeric variable `cluster_var'
                tempvar wt_dev cl_sum cl_tag_t cl_sum_sq
                quietly gen double `wt_dev' = `tsweight' * (`psi' - `ate') if `touse'
                quietly bysort `cluster_var': egen double `cl_sum' = total(`wt_dev') if `touse'
                quietly bysort `cluster_var': gen byte `cl_tag_t' = (_n == 1) if `touse'
                quietly gen double `cl_sum_sq' = `cl_sum'^2 if `cl_tag_t' & `touse'
                quietly summarize `cl_sum_sq' if `cl_tag_t' & `touse'
                local se = sqrt(r(sum)) / `sum_wt'
            }
            else {
                tempvar wt_dev_sq
                quietly gen double `wt_dev_sq' = `tsweight'^2 * (`psi' - `ate')^2 if `touse'
                quietly summarize `wt_dev_sq' if `touse'
                local se = sqrt(r(sum)) / `sum_wt'
            }
        }
        else if "`targetsample'" == "control" {
            /* TMLE for ATC: fit treated regression with weights (1-W.hat)/W.hat */
            tempvar resid1 clever1 ipw_wt
            quietly gen double `resid1' = `depvar' - `yhat1' if `touse' & `treatvar' == 1
            quietly gen double `clever1' = 1 / `whatvar' if `touse' & `treatvar' == 1
            quietly gen double `ipw_wt' = (1 - `whatvar') / `whatvar' if `touse' & `treatvar' == 1

            /* Weighted OLS (no intercept): epsilon = sum(w*X*Y) / sum(w*X^2) */
            tempvar _wxy1 _wxx1
            quietly gen double `_wxy1' = `ipw_wt' * `clever1' * `resid1' if `touse' & `treatvar' == 1
            quietly gen double `_wxx1' = `ipw_wt' * `clever1'^2 if `touse' & `treatvar' == 1
            quietly summarize `_wxy1' if `touse' & `treatvar' == 1
            local _sum_wxy1 = r(sum)
            quietly summarize `_wxx1' if `touse' & `treatvar' == 1
            local epsilon1 = `_sum_wxy1' / r(sum)

            /* Control side: simple OLS */
            tempvar resid0 clever0
            quietly gen double `resid0' = `depvar' - `yhat0' if `touse' & `treatvar' == 0
            quietly gen double `clever0' = 1 / (1 - `whatvar') if `touse' & `treatvar' == 0

            /* Manual OLS (no intercept) */
            tempvar _xy0 _xx0
            quietly gen double `_xy0' = `clever0' * `resid0' if `touse' & `treatvar' == 0
            quietly gen double `_xx0' = `clever0'^2 if `touse' & `treatvar' == 0
            quietly summarize `_xy0' if `touse' & `treatvar' == 0
            local _sum_xy0 = r(sum)
            quietly summarize `_xx0' if `touse' & `treatvar' == 0
            local epsilon0 = `_sum_xy0' / r(sum)

            /* Corrected potential outcomes */
            tempvar yhat0_star yhat1_star
            quietly gen double `yhat0_star' = `yhat0' + `epsilon0' / (1 - `whatvar') if `touse'
            quietly gen double `yhat1_star' = `yhat1' + `epsilon1' / `whatvar' if `touse'

            /* ATC = mean among controls of (yhat1* - yhat0*) + correction */
            quietly count if `touse' & `treatvar' == 0
            local n0 = r(N)

            tempvar psi
            quietly gen double `psi' = . if `touse'
            quietly replace `psi' = (`yhat1_star' - `yhat0_star') ///
                + `treatvar' * (1 - `whatvar') / `whatvar' / (1 - `whatvar') * (`depvar' - `yhat1_star') ///
                - (1 - `treatvar') / (1 - `whatvar') * (`depvar' - `yhat0_star') if `touse'

            /* Weight by control indicator for ATC */
            tempvar tsweight wt_psi
            quietly gen double `tsweight' = (1 - `treatvar') if `touse'
            quietly gen double `wt_psi' = `tsweight' * `psi' if `touse'
            quietly summarize `wt_psi' if `touse'
            local sum_wt = `n0'
            local ate = r(sum) / `sum_wt'

            if "`cluster_var'" != "" {
                confirm numeric variable `cluster_var'
                tempvar wt_dev cl_sum cl_tag_t cl_sum_sq
                quietly gen double `wt_dev' = `tsweight' * (`psi' - `ate') if `touse'
                quietly bysort `cluster_var': egen double `cl_sum' = total(`wt_dev') if `touse'
                quietly bysort `cluster_var': gen byte `cl_tag_t' = (_n == 1) if `touse'
                quietly gen double `cl_sum_sq' = `cl_sum'^2 if `cl_tag_t' & `touse'
                quietly summarize `cl_sum_sq' if `cl_tag_t' & `touse'
                local se = sqrt(r(sum)) / `sum_wt'
            }
            else {
                tempvar wt_dev_sq
                quietly gen double `wt_dev_sq' = `tsweight'^2 * (`psi' - `ate')^2 if `touse'
                quietly summarize `wt_dev_sq' if `touse'
                local se = sqrt(r(sum)) / `sum_wt'
            }
        }
    }

    }

    /* 95% CI and p-value (normal approximation) */
    local ci_lower = `ate' - invnormal(0.975) * `se'
    local ci_upper = `ate' + invnormal(0.975) * `se'
    local tstat    = `ate' / `se'
    local pvalue   = 2 * (1 - normal(abs(`tstat')))

    /* ---- Display results ---- */
    display as text ""
    display as text "Average Treatment Effect (`method')"
    display as text "{hline 60}"
    display as text "Outcome variable:      " as result "`depvar'"
    display as text "Treatment variable:    " as result "`treatvar'"
    display as text "Method:                " as result "`method'"
    display as text "Target sample:         " as result "`targetsample'"
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
    return local  target_sample "`targetsample'"
    return local  method "`method'"
end
