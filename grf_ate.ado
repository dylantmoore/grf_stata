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

    if "`method'" == "TMLE" & "`debiasingweights'" != "" {
        display as error "method(TMLE) does not support debiasingweights()"
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

    /* ---- TMLE + overlap: redirect to AIPW ---- */
    if "`method'" == "TMLE" & "`targetsample'" == "overlap" {
        display as text "Note: TMLE not applicable for overlap weights; using AIPW."
        local method "AIPW"
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
        quietly egen byte `cl_tag' = tag(`cluster_var') if `touse'
        quietly gen double `cl_sum_sq' = `cl_sum'^2 if `cl_tag' == 1 & `touse'
        quietly summarize `cl_sum_sq' if `cl_tag' == 1 & `touse'
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
    else {

    /* ---- Compute TMLE (Targeted Maximum Likelihood Estimation) ----
     *
     * TMLE corrects the initial CATE estimates via targeted bias correction.
     * Note: All regressions are computed manually (OLS formula) to avoid
     * clobbering e() results from grf_causal_forest.
     */

    /* Step 1: Initial potential outcome estimates */
    tempvar yhat0 yhat1
    quietly gen double `yhat0' = `yhatvar' - `tauvar' * `whatvar' if `touse'
    quietly gen double `yhat1' = `yhatvar' + `tauvar' * (1 - `whatvar') if `touse'

    if "`targetsample'" == "all" {
        /* ---- TMLE for ATE (all) ----
         * Controls (W=0): regress (Y - Y.hat.0) on 1/(1-W.hat), no intercept => epsilon0
         * Treated (W=1): regress (Y - Y.hat.1) on 1/W.hat, no intercept => epsilon1
         * DR correction uses full-sample means of clever covariates (matching R's grf).
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

        /* Step 5: DR correction -- full-sample means (R: mean(1/e), mean(1/(1-e))) */
        quietly summarize `tauvar' if `touse'
        local tau_avg = r(mean)

        tempvar inv_e inv_1me
        quietly gen double `inv_e' = 1 / `whatvar' if `touse'
        quietly gen double `inv_1me' = 1 / (1 - `whatvar') if `touse'
        quietly summarize `inv_e' if `touse'
        local mean_clever1 = r(mean)
        quietly summarize `inv_1me' if `touse'
        local mean_clever0 = r(mean)

        local dr_correction = `epsilon1' * `mean_clever1' - `epsilon0' * `mean_clever0'

        /* Step 6: ATE */
        local ate = `tau_avg' + `dr_correction'

        /* Step 7: Sandwich SE (standard TMLE influence function) */
        tempvar psi yhat0_star yhat1_star
        quietly gen double `yhat0_star' = `yhat0' + `epsilon0' / (1 - `whatvar') if `touse'
        quietly gen double `yhat1_star' = `yhat1' + `epsilon1' / `whatvar' if `touse'
        quietly gen double `psi' = (`yhat1_star' - `yhat0_star') - `ate' ///
            + `treatvar' / `whatvar' * (`depvar' - `yhat1_star') ///
            - (1 - `treatvar') / (1 - `whatvar') * (`depvar' - `yhat0_star') if `touse'

        if "`cluster_var'" != "" {
            confirm numeric variable `cluster_var'
            tempvar cl_psi cl_tag_t cl_psi_sq
            quietly bysort `cluster_var': egen double `cl_psi' = total(`psi') if `touse'
            quietly egen byte `cl_tag_t' = tag(`cluster_var') if `touse'
            quietly gen double `cl_psi_sq' = `cl_psi'^2 if `cl_tag_t' == 1 & `touse'
            quietly summarize `cl_psi_sq' if `cl_tag_t' == 1 & `touse'
            local se = sqrt(r(sum)) / `n_use'
        }
        else {
            quietly summarize `psi' if `touse'
            local se = r(sd) / sqrt(`n_use')
        }
    }
    else if "`targetsample'" == "treated" {
        /* ---- TMLE for ATT ----
         * Control-side OLS (unweighted, no intercept):
         *   clever0 = e/(1-e), resid0 = Y - Y.hat.0, among W=0
         *   epsilon0 = sum(clever0 * resid0) / sum(clever0^2)
         * Point estimate: tau_att_raw - epsilon0 * new_center
         * SE: HC1 sandwich from control OLS + treated variance
         */
        quietly count if `touse' & `treatvar' == 1
        local n1 = r(N)
        quietly count if `touse' & `treatvar' == 0
        local n0 = r(N)

        /* Control-side regression */
        tempvar resid0 clever0
        quietly gen double `resid0' = `depvar' - `yhat0' if `touse' & `treatvar' == 0
        quietly gen double `clever0' = `whatvar' / (1 - `whatvar') if `touse' & `treatvar' == 0

        /* Unweighted OLS (no intercept): epsilon0 = sum(clever0*resid0)/sum(clever0^2) */
        tempvar _xy0 _xx0
        quietly gen double `_xy0' = `clever0' * `resid0' if `touse' & `treatvar' == 0
        quietly gen double `_xx0' = `clever0'^2 if `touse' & `treatvar' == 0
        quietly summarize `_xy0' if `touse' & `treatvar' == 0
        local _sum_xy0 = r(sum)
        quietly summarize `_xx0' if `touse' & `treatvar' == 0
        local _sum_xx0 = r(sum)
        local epsilon0 = `_sum_xy0' / `_sum_xx0'

        /* Point estimate */
        quietly summarize `tauvar' if `touse' & `treatvar' == 1
        local tau_att_raw = r(mean)

        tempvar att_clever_t
        quietly gen double `att_clever_t' = `whatvar' / (1 - `whatvar') if `touse' & `treatvar' == 1
        quietly summarize `att_clever_t' if `touse' & `treatvar' == 1
        local new_center = r(mean)

        local ate = `tau_att_raw' - `epsilon0' * `new_center'

        /* SE: HC1 sandwich */
        tempvar ols_resid
        quietly gen double `ols_resid' = `resid0' - `epsilon0' * `clever0' if `touse' & `treatvar' == 0

        if "`cluster_var'" != "" {
            confirm numeric variable `cluster_var'
            /* Clustered V_eps */
            tempvar score0 cl_score0 cl_tag0 cl_score0_sq
            quietly gen double `score0' = `ols_resid' * `clever0' if `touse' & `treatvar' == 0
            quietly bysort `cluster_var': egen double `cl_score0' = total(`score0') if `touse' & `treatvar' == 0
            quietly egen byte `cl_tag0' = tag(`cluster_var') if `touse' & `treatvar' == 0
            quietly gen double `cl_score0_sq' = `cl_score0'^2 if `cl_tag0' == 1 & `touse' & `treatvar' == 0
            quietly summarize `cl_score0_sq' if `cl_tag0' == 1 & `touse' & `treatvar' == 0
            local _sum_cl_score0_sq = r(sum)
            quietly count if `cl_tag0' == 1 & `touse' & `treatvar' == 0
            local G0 = r(N)
            local V_eps = (`G0' / (`G0' - 1)) * `_sum_cl_score0_sq' / (`_sum_xx0'^2)

            /* Clustered V_treat */
            tempvar tresid cl_tresid cl_tag1 cl_nobs1 cl_tresid_dev cl_tresid_sq
            quietly gen double `tresid' = `depvar' - `yhat1' if `touse' & `treatvar' == 1
            quietly bysort `cluster_var': egen double `cl_tresid' = total(`tresid') if `touse' & `treatvar' == 1
            quietly bysort `cluster_var': egen long `cl_nobs1' = count(`tresid') if `touse' & `treatvar' == 1
            quietly egen byte `cl_tag1' = tag(`cluster_var') if `touse' & `treatvar' == 1
            quietly summarize `tresid' if `touse' & `treatvar' == 1
            local tresid_mean = r(mean)
            quietly gen double `cl_tresid_dev' = `cl_tresid' - `cl_nobs1' * `tresid_mean' if `cl_tag1' == 1 & `touse' & `treatvar' == 1
            quietly gen double `cl_tresid_sq' = `cl_tresid_dev'^2 if `cl_tag1' == 1 & `touse' & `treatvar' == 1
            quietly summarize `cl_tresid_sq' if `cl_tag1' == 1 & `touse' & `treatvar' == 1
            quietly count if `cl_tag1' == 1 & `touse' & `treatvar' == 1
            local G1 = r(N)
            local V_treat = (`G1' / (`G1' - 1)) * r(sum) / (`n1'^2)

            local sigma2 = `V_eps' * `new_center'^2 + `V_treat'
            local se = sqrt(`sigma2')
        }
        else {
            /* iid SE */
            tempvar r_sq_x_sq
            quietly gen double `r_sq_x_sq' = `ols_resid'^2 * `clever0'^2 if `touse' & `treatvar' == 0
            quietly summarize `r_sq_x_sq' if `touse' & `treatvar' == 0
            local V_eps = (`n0' / (`n0' - 1)) * r(sum) / (`_sum_xx0'^2)

            tempvar tresid
            quietly gen double `tresid' = `depvar' - `yhat1' if `touse' & `treatvar' == 1
            quietly summarize `tresid' if `touse' & `treatvar' == 1
            local V_treat = r(Var) / `n1'

            local sigma2 = `V_eps' * `new_center'^2 + `V_treat'
            local se = sqrt(`sigma2')
        }
    }
    else if "`targetsample'" == "control" {
        /* ---- TMLE for ATC ----
         * Treated-side OLS (unweighted, no intercept):
         *   clever1 = (1-e)/e, resid1 = Y - Y.hat.1, among W=1
         *   epsilon1 = sum(clever1 * resid1) / sum(clever1^2)
         * Point estimate: tau_atc_raw + epsilon1 * new_center
         * SE: HC1 sandwich from treated OLS + control variance
         */
        quietly count if `touse' & `treatvar' == 1
        local n1 = r(N)
        quietly count if `touse' & `treatvar' == 0
        local n0 = r(N)

        /* Treated-side regression */
        tempvar resid1 clever1
        quietly gen double `resid1' = `depvar' - `yhat1' if `touse' & `treatvar' == 1
        quietly gen double `clever1' = (1 - `whatvar') / `whatvar' if `touse' & `treatvar' == 1

        /* Unweighted OLS (no intercept): epsilon1 = sum(clever1*resid1)/sum(clever1^2) */
        tempvar _xy1 _xx1
        quietly gen double `_xy1' = `clever1' * `resid1' if `touse' & `treatvar' == 1
        quietly gen double `_xx1' = `clever1'^2 if `touse' & `treatvar' == 1
        quietly summarize `_xy1' if `touse' & `treatvar' == 1
        local _sum_xy1 = r(sum)
        quietly summarize `_xx1' if `touse' & `treatvar' == 1
        local _sum_xx1 = r(sum)
        local epsilon1 = `_sum_xy1' / `_sum_xx1'

        /* Point estimate */
        quietly summarize `tauvar' if `touse' & `treatvar' == 0
        local tau_atc_raw = r(mean)

        tempvar atc_clever_c
        quietly gen double `atc_clever_c' = (1 - `whatvar') / `whatvar' if `touse' & `treatvar' == 0
        quietly summarize `atc_clever_c' if `touse' & `treatvar' == 0
        local new_center = r(mean)

        local ate = `tau_atc_raw' + `epsilon1' * `new_center'

        /* SE: HC1 sandwich */
        tempvar ols_resid
        quietly gen double `ols_resid' = `resid1' - `epsilon1' * `clever1' if `touse' & `treatvar' == 1

        if "`cluster_var'" != "" {
            confirm numeric variable `cluster_var'
            /* Clustered V_eps */
            tempvar score1 cl_score1 cl_tag1 cl_score1_sq
            quietly gen double `score1' = `ols_resid' * `clever1' if `touse' & `treatvar' == 1
            quietly bysort `cluster_var': egen double `cl_score1' = total(`score1') if `touse' & `treatvar' == 1
            quietly egen byte `cl_tag1' = tag(`cluster_var') if `touse' & `treatvar' == 1
            quietly gen double `cl_score1_sq' = `cl_score1'^2 if `cl_tag1' == 1 & `touse' & `treatvar' == 1
            quietly summarize `cl_score1_sq' if `cl_tag1' == 1 & `touse' & `treatvar' == 1
            local _sum_cl_score1_sq = r(sum)
            quietly count if `cl_tag1' == 1 & `touse' & `treatvar' == 1
            local G1 = r(N)
            local V_eps = (`G1' / (`G1' - 1)) * `_sum_cl_score1_sq' / (`_sum_xx1'^2)

            /* Clustered V_ctrl */
            tempvar cresid cl_cresid cl_tag0 cl_nobs0 cl_cresid_dev cl_cresid_sq
            quietly gen double `cresid' = `depvar' - `yhat0' if `touse' & `treatvar' == 0
            quietly bysort `cluster_var': egen double `cl_cresid' = total(`cresid') if `touse' & `treatvar' == 0
            quietly bysort `cluster_var': egen long `cl_nobs0' = count(`cresid') if `touse' & `treatvar' == 0
            quietly egen byte `cl_tag0' = tag(`cluster_var') if `touse' & `treatvar' == 0
            quietly summarize `cresid' if `touse' & `treatvar' == 0
            local cresid_mean = r(mean)
            quietly gen double `cl_cresid_dev' = `cl_cresid' - `cl_nobs0' * `cresid_mean' if `cl_tag0' == 1 & `touse' & `treatvar' == 0
            quietly gen double `cl_cresid_sq' = `cl_cresid_dev'^2 if `cl_tag0' == 1 & `touse' & `treatvar' == 0
            quietly summarize `cl_cresid_sq' if `cl_tag0' == 1 & `touse' & `treatvar' == 0
            quietly count if `cl_tag0' == 1 & `touse' & `treatvar' == 0
            local G0 = r(N)
            local V_ctrl = (`G0' / (`G0' - 1)) * r(sum) / (`n0'^2)

            local sigma2 = `V_eps' * `new_center'^2 + `V_ctrl'
            local se = sqrt(`sigma2')
        }
        else {
            /* iid SE */
            tempvar r_sq_x_sq
            quietly gen double `r_sq_x_sq' = `ols_resid'^2 * `clever1'^2 if `touse' & `treatvar' == 1
            quietly summarize `r_sq_x_sq' if `touse' & `treatvar' == 1
            local V_eps = (`n1' / (`n1' - 1)) * r(sum) / (`_sum_xx1'^2)

            tempvar cresid
            quietly gen double `cresid' = `depvar' - `yhat0' if `touse' & `treatvar' == 0
            quietly summarize `cresid' if `touse' & `treatvar' == 0
            local V_ctrl = r(Var) / `n0'

            local sigma2 = `V_eps' * `new_center'^2 + `V_ctrl'
            local se = sqrt(`sigma2')
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
