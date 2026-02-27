*! grf_best_linear_projection.ado -- Best Linear Projection of CATE
*! Version 0.3.0
*! Projects AIPW DR scores onto covariates via OLS with robust SEs
*! Matches grf::best_linear_projection() which uses DR scores, not raw tau_hat
*! Supports vcov.type (HC0/HC1/HC2/HC3) and target.sample (all/overlap)
*!
*! NOTE: This is an eclass command. Running it replaces the causal forest
*! e() results with regression output (same as R's best_linear_projection).
*! Run grf_ate/grf_test_calibration BEFORE this command if needed.

program define grf_best_linear_projection, eclass
    version 14.0

    syntax [varlist(numeric default=none)] [if] [in] ///
        [, VCOVtype(string) TARGETsample(string)]

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

    /* ---- Default and validate vcov.type ---- */
    if "`vcovtype'" == "" {
        local vcovtype "HC3"
    }
    local vcovtype = upper("`vcovtype'")
    if "`vcovtype'" != "HC0" & "`vcovtype'" != "HC1" ///
        & "`vcovtype'" != "HC2" & "`vcovtype'" != "HC3" {
        display as error "vcov.type must be one of: HC0, HC1, HC2, HC3"
        exit 198
    }

    /* ---- Default and validate target.sample ---- */
    if "`targetsample'" == "" {
        local targetsample "all"
    }
    if "`targetsample'" != "all" & "`targetsample'" != "overlap" {
        display as error "target.sample must be one of: all, overlap"
        exit 198
    }

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

    /* ---- Compute observation weights for target.sample ---- */
    tempvar obsweight
    if "`targetsample'" == "overlap" {
        quietly gen double `obsweight' = `whatvar' * (1 - `whatvar') if `touse'
        /* Drop observations with near-zero weight */
        quietly replace `touse' = 0 if `obsweight' < 1e-12 & `touse'
        quietly count if `touse'
        local n_use = r(N)
        if `n_use' < 2 {
            display as error "too few observations with positive overlap weights"
            exit 2000
        }
    }
    else {
        quietly gen double `obsweight' = 1 if `touse'
    }

    /* ---- Display header ---- */
    display as text ""
    display as text "Best Linear Projection of CATE (AIPW)"
    display as text "{hline 55}"
    display as text "CATE variable:         " as result "`tauvar'"
    display as text "Projection covariates: " as result "`proj_vars'"
    display as text "Target sample:         " as result "`targetsample'"
    display as text "VCOV type:             " as result "`vcovtype'"
    display as text "Observations:          " as result `n_use'
    display as text "{hline 55}"
    display as text ""
    display as text "Regression: DR_score ~ covariates, vce(`vcovtype')"
    display as text ""

    /* ---- Run OLS with specified vcov.type via Mata ---- */
    /* For target.sample="all" with HC3, use Stata's built-in for simplicity */
    if "`targetsample'" == "all" & "`vcovtype'" == "HC3" {
        regress `dr_score' `proj_vars' if `touse', vce(hc3)
    }
    else if "`targetsample'" == "all" & "`vcovtype'" == "HC2" {
        regress `dr_score' `proj_vars' if `touse', vce(hc2)
    }
    else if "`targetsample'" == "all" & "`vcovtype'" == "HC1" {
        regress `dr_score' `proj_vars' if `touse', vce(robust)
    }
    else if "`targetsample'" == "all" & "`vcovtype'" == "HC0" {
        /* HC0 needs manual computation; run OLS first for coefficients */
        quietly regress `dr_score' `proj_vars' if `touse'
        local k = e(df_m) + 1
        /* Repost with HC0 vcov via Mata */
        tempname b_ols V_hc0
        matrix `b_ols' = e(b)
        mata: _grf_blp_hc0("`dr_score'", "`proj_vars'", "`touse'", ///
            "`V_hc0'")
        ereturn repost V = `V_hc0'
        ereturn display
    }
    else {
        /* target.sample="overlap": weighted regression via Mata */
        /* Normalize weights to sum to n for proper regression */
        quietly summarize `obsweight' if `touse'
        local wt_sum = r(sum)
        tempvar nweight
        quietly gen double `nweight' = `obsweight' * `n_use' / `wt_sum' if `touse'

        /* Run weighted regression, then compute specified HC vcov */
        quietly regress `dr_score' `proj_vars' [aweight=`nweight'] if `touse'
        local k = e(df_m) + 1

        tempname V_wt
        mata: _grf_blp_weighted("`dr_score'", "`proj_vars'", ///
            "`nweight'", "`touse'", "`vcovtype'", "`V_wt'")
        ereturn repost V = `V_wt'
        ereturn display
    }

    /* Store additional info */
    ereturn local vcov_type     "`vcovtype'"
    ereturn local target_sample "`targetsample'"

end

/* ====================================================================
 * Mata helper: HC0 sandwich variance (unweighted)
 * ==================================================================== */
mata:
void _grf_blp_hc0(string scalar yvar, string scalar xvars,
                   string scalar tousevar, string scalar vname)
{
    real matrix X, XpX_inv, V
    real colvector y, e
    real scalar n, k, i

    st_view(y, ., yvar, tousevar)
    st_view(X, ., tokens(xvars), tousevar)
    /* Add constant */
    X = X, J(rows(X), 1, 1)
    n = rows(X)
    k = cols(X)

    XpX_inv = invsym(cross(X, X))
    e = y - X * (XpX_inv * cross(X, y))

    /* HC0: (X'X)^{-1} X' diag(e^2) X (X'X)^{-1} */
    V = XpX_inv * cross(X, e:^2, X) * XpX_inv

    st_matrix(vname, V)
    st_matrixcolstripe(vname, st_matrixcolstripe("e(V)"))
    st_matrixrowstripe(vname, st_matrixrowstripe("e(V)"))
}
end

/* ====================================================================
 * Mata helper: HC0/HC1/HC2/HC3 sandwich variance (weighted)
 * ==================================================================== */
mata:
void _grf_blp_weighted(string scalar yvar, string scalar xvars,
                       string scalar wtvar, string scalar tousevar,
                       string scalar hctype, string scalar vname)
{
    real matrix X, XpX_inv, H, V, meat
    real colvector y, w, sw, e, h
    real scalar n, k, i

    st_view(y, ., yvar, tousevar)
    st_view(X, ., tokens(xvars), tousevar)
    st_view(w, ., wtvar, tousevar)
    /* Add constant */
    X = X, J(rows(X), 1, 1)
    n = rows(X)
    k = cols(X)

    /* Weighted OLS: (X'WX)^{-1} X'Wy */
    sw = sqrt(w)
    XpX_inv = invsym(cross(sw:*X, sw:*X))
    e = y - X * (XpX_inv * cross(sw:*X, sw:*y))

    if (hctype == "HC0") {
        meat = cross(sw:*X, (w:*(e:^2)), X)
        V = XpX_inv * meat * XpX_inv
    }
    else if (hctype == "HC1") {
        meat = cross(sw:*X, (w:*(e:^2)), X)
        V = (n / (n - k)) * XpX_inv * meat * XpX_inv
    }
    else if (hctype == "HC2") {
        /* h_ii = x_i' (X'WX)^{-1} x_i * w_i */
        h = J(n, 1, 0)
        for (i = 1; i <= n; i++) {
            h[i] = (sw[i] * X[i,.]) * XpX_inv * (sw[i] * X[i,.])'
        }
        meat = cross(sw:*X, (w:*(e:^2 :/ (1 :- h))), X)
        V = XpX_inv * meat * XpX_inv
    }
    else {
        /* HC3 */
        h = J(n, 1, 0)
        for (i = 1; i <= n; i++) {
            h[i] = (sw[i] * X[i,.]) * XpX_inv * (sw[i] * X[i,.])'
        }
        meat = cross(sw:*X, (w:*(e:^2 :/ (1 :- h):^2)), X)
        V = XpX_inv * meat * XpX_inv
    }

    st_matrix(vname, V)
    st_matrixcolstripe(vname, st_matrixcolstripe("e(V)"))
    st_matrixrowstripe(vname, st_matrixrowstripe("e(V)"))
}
end
