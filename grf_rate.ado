*! grf_rate.ado -- RATE (Rank-Weighted Average Treatment Effect)
*! Version 0.2.0
*! Evaluates treatment prioritization rules
*! Implements AUTOC and QINI from Yadlowsky, Fleming, Shah, Brunskill, Wager (2021)
*! Supports compliance.score for instrumental forest RATE evaluation
*! Note: Full R name is rank_average_treatment_effect; shortened for Stata's
*!       32-character command name limit.

program define grf_rate, rclass
    version 14.0

    syntax varname(numeric) [if] [in],  ///
        [                               ///
            TARget(string)              ///
            Quantiles(numlist >0 <=1)   ///
            BOOTstrap(integer 200)      ///
            CATEvar(varname numeric)    ///
            COMPliancescore(varname numeric) ///
            SEED(integer -1)            ///
        ]

    /* ---- Defaults ---- */
    if "`target'" == "" {
        local target "AUTOC"
    }

    /* Validate target */
    local target = upper("`target'")
    if "`target'" != "AUTOC" & "`target'" != "QINI" {
        display as error "target() must be AUTOC or QINI"
        exit 198
    }

    /* Default quantile grid */
    if "`quantiles'" == "" {
        local quantiles "0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0"
    }

    /* Validate bootstrap reps */
    if `bootstrap' < 2 {
        display as error "bootstrap() must be at least 2"
        exit 198
    }

    /* ---- Determine CATE variable ---- */
    local priorities "`varlist'"

    if "`catevar'" != "" {
        local tauvar "`catevar'"
    }
    else {
        /* Read from e() results of grf_causal_forest */
        if "`e(cmd)'" != "grf_causal_forest" {
            display as error "grf_rate requires prior estimation" ///
                " by grf_causal_forest, or specify catevar()"
            exit 301
        }
        local tauvar "`e(predict_var)'"
        if "`tauvar'" == "" {
            display as error "no CATE variable found in e() results; specify catevar()"
            exit 301
        }
    }

    /* ---- Confirm required variables exist ---- */
    confirm numeric variable `tauvar'
    confirm numeric variable `priorities'

    /* ---- Validate compliance score if specified ---- */
    if "`compliancescore'" != "" {
        confirm numeric variable `compliancescore'
    }

    /* ---- Compute AIPW doubly-robust scores ----
     *
     * R's rank_average_treatment_effect uses DR scores, not raw CATE.
     * DR_i = tau_hat_i + (W_i - W_hat_i) / Var(W-W_hat) *
     *        (Y_i - Y_hat_i - tau_hat_i * (W_i - W_hat_i))
     *
     * This requires access to Y, W, Y.hat, W.hat from the causal forest.
     * If available, we use DR scores. Otherwise, fall back to raw CATE.
     */
    local use_dr 0
    local depvar    "`e(depvar)'"
    local treatvar  "`e(treatvar)'"
    local yhatvar   "`e(yhat_var)'"
    local whatvar   "`e(what_var)'"

    if "`depvar'" != "" & "`treatvar'" != "" & "`yhatvar'" != "" & "`whatvar'" != "" {
        capture confirm numeric variable `depvar'
        if !_rc {
            capture confirm numeric variable `treatvar'
            if !_rc {
                capture confirm numeric variable `yhatvar'
                if !_rc {
                    capture confirm numeric variable `whatvar'
                    if !_rc {
                        local use_dr 1
                    }
                }
            }
        }
    }

    /* ---- Mark sample ---- */
    marksample touse
    markout `touse' `tauvar' `priorities'
    if `use_dr' {
        markout `touse' `depvar' `treatvar' `yhatvar' `whatvar'
    }
    if "`compliancescore'" != "" {
        markout `touse' `compliancescore'
    }
    quietly count if `touse'
    local n_use = r(N)

    if `n_use' < 10 {
        display as error "need at least 10 non-missing observations"
        exit 2000
    }

    /* ---- Compute DR scores or use raw CATE ---- */
    if `use_dr' {
        tempvar w_resid y_resid dr_score
        quietly gen double `w_resid' = `treatvar' - `whatvar' if `touse'
        quietly summarize `w_resid' if `touse'
        local w_resid_var = r(Var)

        if `w_resid_var' < 1e-12 {
            display as text "Warning: Var(W-W.hat) near zero; falling back to raw CATE"
            local use_dr 0
        }
        else {
            quietly {
                gen double `y_resid' = `depvar' - `yhatvar' - `tauvar' * `w_resid' if `touse'
                gen double `dr_score' = `tauvar' + (`w_resid' / `w_resid_var') * `y_resid' if `touse'
            }
        }
    }

    /* The scoring variable: DR scores if available, raw CATE otherwise */
    if `use_dr' {
        local scorevar "`dr_score'"
        local score_label "DR scores"
    }
    else {
        local scorevar "`tauvar'"
        local score_label "raw CATE"
    }

    /* ---- Apply compliance score weighting if specified ---- */
    local use_compliance 0
    if "`compliancescore'" != "" {
        local use_compliance 1
        tempvar cscore_wt
        quietly gen double `cscore_wt' = `scorevar' * `compliancescore' if `touse'
        local scorevar "`cscore_wt'"
        local score_label "`score_label' (compliance-weighted)"
    }

    /* ---- Set seed if requested ---- */
    if `seed' >= 0 {
        set seed `seed'
    }

    /* ---- Count quantiles ---- */
    local n_quantiles 0
    foreach q of local quantiles {
        local n_quantiles = `n_quantiles' + 1
    }

    /* ---- Display header ---- */
    display as text ""
    display as text "Rank-Weighted Average Treatment Effect"
    display as text "{hline 55}"
    display as text "Target:                " as result "`target'"
    display as text "Priorities variable:   " as result "`priorities'"
    display as text "CATE variable:         " as result "`tauvar'"
    display as text "Scoring:               " as result "`score_label'"
    if "`compliancescore'" != "" {
        display as text "Compliance score:      " as result "`compliancescore'"
    }
    display as text "Observations:          " as result `n_use'
    display as text "Bootstrap reps:        " as result `bootstrap'
    display as text "Quantile grid points:  " as result `n_quantiles'
    display as text "{hline 55}"
    display as text ""

    /* ====================================================================
     * RATE Computation
     *
     * The TOC (Targeting Operator Characteristic) at fraction q is:
     *   TOC(q) = mean(Score_i | rank(S_i) in top q-fraction) - mean(Score_i)
     *
     * where Score = DR scores (AIPW) if available, raw CATE otherwise.
     *
     * AUTOC = integral_0^1 TOC(q) dq
     * QINI  = integral_0^1 q * TOC(q) dq
     *
     * We approximate the integrals using a trapezoidal rule on the
     * quantile grid.
     * ==================================================================== */

    /* ---- Compute point estimate ---- */
    display as text "Computing point estimate ..."

    /* Compute overall mean score */
    quietly summarize `scorevar' if `touse'
    local tau_mean = r(mean)

    /* Create rank variable based on priorities (descending) */
    tempvar rank_var fractional_rank
    quietly egen `rank_var' = rank(`priorities') if `touse', field
    /* field gives rank 1 to highest value */
    quietly gen double `fractional_rank' = `rank_var' / `n_use' if `touse'

    /* Compute TOC at each quantile point */
    /* Store TOC values in a temporary variable indexed by quantile */
    local toc_values ""
    local prev_q 0
    local qi 0

    foreach q of local quantiles {
        local qi = `qi' + 1

        /* Mean score among top q-fraction */
        quietly summarize `scorevar' if `touse' & `fractional_rank' <= `q'
        local toc_`qi' = r(mean) - `tau_mean'

        local toc_values "`toc_values' `toc_`qi''"
    }

    /* Compute integral using trapezoidal rule */
    local rate_est = 0
    local prev_q   = 0
    local prev_toc = 0
    local qi = 0

    foreach q of local quantiles {
        local qi = `qi' + 1
        local cur_toc = `toc_`qi''
        local dq = `q' - `prev_q'

        if "`target'" == "AUTOC" {
            /* Trapezoidal: 0.5 * (f(a) + f(b)) * (b - a) */
            local rate_est = `rate_est' + 0.5 * (`prev_toc' + `cur_toc') * `dq'
        }
        else {
            /* QINI: integral of q * TOC(q) dq */
            /* Trapezoidal: 0.5 * (a*f(a) + b*f(b)) * (b - a) */
            local rate_est = `rate_est' + 0.5 * (`prev_q' * `prev_toc' + `q' * `cur_toc') * `dq'
        }

        local prev_q   = `q'
        local prev_toc = `cur_toc'
    }

    /* ---- Bootstrap standard errors ---- */
    display as text "Computing bootstrap standard errors (`bootstrap' replications) ..."

    /* We store bootstrap estimates in a temporary variable */
    tempname bs_mat
    matrix `bs_mat' = J(`bootstrap', 1, .)

    /* Preserve current data state for resampling */
    tempvar orig_order bs_weight

    /* Generate original order for restoring sort */
    quietly gen long `orig_order' = _n

    forvalues b = 1/`bootstrap' {

        /* Generate bootstrap weights (Poisson bootstrap for speed) */
        /* Each obs gets a Poisson(1) weight -- equivalent to multinomial resampling */
        quietly {
            capture drop `bs_weight'
            gen long `bs_weight' = rpoisson(1) if `touse'

            /* Compute weighted overall mean score */
            summarize `scorevar' [fw=`bs_weight'] if `touse' & `bs_weight' > 0
            local bs_tau_mean = r(mean)
            local bs_n = r(sum_w)

            if `bs_n' < 2 {
                matrix `bs_mat'[`b', 1] = .
                continue
            }

            /* Rank by priorities within bootstrap sample (weighted) */
            /* For weighted ranks, we use the priorities directly and
               compute the fraction based on cumulative weights */

            /* Sort by priorities descending to compute cumulative weight fractions */
            tempvar sort_key cum_wt frac_wt
            gen double `sort_key' = -`priorities' if `touse' & `bs_weight' > 0
            sort `sort_key'

            gen double `cum_wt' = sum(`bs_weight') if `touse' & `bs_weight' > 0
            gen double `frac_wt' = `cum_wt' / `bs_n' if `touse' & `bs_weight' > 0
        }

        /* Compute TOC at each quantile for this bootstrap sample */
        local bs_prev_q   = 0
        local bs_prev_toc = 0
        local bs_rate     = 0
        local bqi = 0

        foreach q of local quantiles {
            local bqi = `bqi' + 1

            quietly summarize `scorevar' [fw=`bs_weight'] ///
                if `touse' & `bs_weight' > 0 & `frac_wt' <= `q'
            if r(N) > 0 {
                local bs_toc = r(mean) - `bs_tau_mean'
            }
            else {
                local bs_toc = 0
            }

            local dq = `q' - `bs_prev_q'

            if "`target'" == "AUTOC" {
                local bs_rate = `bs_rate' + 0.5 * (`bs_prev_toc' + `bs_toc') * `dq'
            }
            else {
                local bs_rate = `bs_rate' + 0.5 * (`bs_prev_q' * `bs_prev_toc' + `q' * `bs_toc') * `dq'
            }

            local bs_prev_q   = `q'
            local bs_prev_toc = `bs_toc'
        }

        matrix `bs_mat'[`b', 1] = `bs_rate'

        /* Clean up temp vars from this iteration */
        quietly {
            capture drop `sort_key'
            capture drop `cum_wt'
            capture drop `frac_wt'
        }

        /* Progress indicator every 50 reps */
        if mod(`b', 50) == 0 {
            display as text "  ... completed `b'/`bootstrap' bootstrap replications"
        }
    }

    /* Restore original sort order */
    sort `orig_order'

    /* Clean up bootstrap weight */
    capture drop `bs_weight'

    /* ---- Compute bootstrap SE ---- */
    /* Extract bootstrap estimates into a temp variable for summarize */
    tempvar bs_vals
    quietly gen double `bs_vals' = .

    local n_valid_bs = 0
    forvalues b = 1/`bootstrap' {
        local bval = `bs_mat'[`b', 1]
        if !missing(`bval') {
            local n_valid_bs = `n_valid_bs' + 1
            quietly replace `bs_vals' = `bval' in `n_valid_bs'
        }
    }

    if `n_valid_bs' < 2 {
        display as error "too few valid bootstrap replications to compute SE"
        exit 498
    }

    quietly summarize `bs_vals' in 1/`n_valid_bs'
    local bs_se = r(sd)

    /* ---- Compute z-statistic and p-value ---- */
    local z_stat = `rate_est' / `bs_se'
    local p_value = 2 * (1 - normal(abs(`z_stat')))

    /* ---- Display results ---- */
    display as text ""
    display as text "Rank-Weighted Average Treatment Effect"
    display as text "{hline 55}"
    display as text "Target:                " as result "`target'"
    display as text "Observations:          " as result `n_use'
    display as text "Bootstrap reps:        " as result `n_valid_bs'
    display as text "{hline 55}"
    display as text ""
    display as text %~14s " " %~12s "Estimate" %~12s "Std.Err." ///
        %~8s "z" %~10s "P>|z|"
    display as text "{hline 55}"
    display as text "`target'" _col(10) ///
        as result %12.6f `rate_est' ///
        _col(24) as result %12.6f `bs_se' ///
        _col(38) as result %8.2f `z_stat' ///
        _col(49) as result %8.4f `p_value'
    display as text "{hline 55}"
    display as text ""

    /* ---- Store results ---- */
    return scalar estimate    = `rate_est'
    return scalar std_err     = `bs_se'
    return scalar z_stat      = `z_stat'
    return scalar p_value     = `p_value'
    return scalar n           = `n_use'
    return scalar n_bootstrap = `n_valid_bs'
    return local  target        "`target'"
    return local  priorities    "`priorities'"
    return local  catevar       "`tauvar'"
    if "`compliancescore'" != "" {
        return local compliance_score_var "`compliancescore'"
    }
end
