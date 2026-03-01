*! grf_lm_forest.ado -- Linear Model Forest via grf C++ library
*! Version 0.1.0
*! Implements Friedberg et al. (2020) lm_forest()
*!
*! Estimates Y = c(x) + h_1(x)*W_1 + ... + h_K(x)*W_K where
*! h_k(x) are heterogeneous coefficients estimated via GRF.

program define grf_lm_forest, eclass
    version 14.0

    syntax varlist(min=2 numeric) [if] [in],  ///
        GENerate(name)                         ///
        Xvars(varlist numeric)                 ///
        [                                      ///
            NTrees(integer 2000)               ///
            NUISancetrees(integer 500)         ///
            SEED(integer 42)                   ///
            MTRY(integer 0)                    ///
            MINNodesize(integer 5)             ///
            SAMPLEfrac(real 0.5)               ///
            noHONesty                          ///
            HONestyfrac(real 0.5)              ///
            noHONestyprune                     ///
            STABilizesplits                    ///
            ALPha(real 0.05)                   ///
            IMBalancepenalty(real 0.0)          ///
            CIGroupsize(integer 1)             ///
            NUMThreads(integer 0)              ///
            ESTIMATEVariance                   ///
            REPlace                            ///
            noMIA                              ///
            YHATinput(varname numeric)          ///
            WHATinput(varlist numeric)          ///
            YHATgenerate(name)                 ///
            WHATgenerate(namelist)              ///
            GRADientweights(numlist)            ///
            CLuster(varname numeric)           ///
            WEIghts(varname numeric)           ///
            EQUALizeclusterweights             ///
            TUNEParameters(string)             ///
            TUNENumtrees(integer 200)          ///
            TUNENumreps(integer 50)            ///
        ]

    /* ---- Parse honesty ---- */
    local do_honesty 1
    if "`honesty'" == "nohonesty" {
        local do_honesty 0
    }

    /* ---- Parse honesty prune ---- */
    local do_honesty_prune 1
    if "`honestyprune'" == "nohonestyprune" {
        local do_honesty_prune 0
    }

    /* ---- Parse stabilize_splits ---- */
    /* R's lm_forest defaults stabilize.splits = FALSE.
     * STABilizesplits is opt-in (default OFF). */
    local do_stabilize 0
    if "`stabilizesplits'" != "" {
        local do_stabilize 1
    }

    /* ---- Parse gradient_weights ---- */
    /* R's lm_forest accepts gradient.weights as a vector of weights for
     * gradient computation. If not specified, an empty vector is passed
     * (C++ default behavior). */
    local gradient_weights_str ""
    if "`gradientweights'" != "" {
        /* Validate: should have exactly n_regressors * n_y weights */
        local n_gw : word count `gradientweights'
        /* Convert numlist to space-separated string for argv */
        local gradient_weights_str "`gradientweights'"
    }

    /* ---- Parse estimate_variance ---- */
    local do_est_var 0
    if "`estimatevariance'" != "" {
        local do_est_var 1
        if `cigroupsize' < 2 {
            local cigroupsize 2
        }
    }

    /* ---- Parse MIA ---- */
    local allow_missing_x 1
    if "`mia'" == "nomia" {
        local allow_missing_x 0
    }

    /* ---- Parse cluster ---- */
    local cluster_col_idx 0
    if "`cluster'" != "" {
        local cluster_var `cluster'
    }

    /* ---- Parse weights ---- */
    local weight_col_idx 0
    if "`weights'" != "" {
        local weight_var `weights'
    }

    /* ---- Parse varlist: depvar regressor1 [regressor2 ...] ---- */
    gettoken depvar regvars : varlist
    local n_regressors : word count `regvars'

    if `n_regressors' < 1 {
        display as error "need at least 1 regressor variable (W)"
        exit 198
    }

    /* ---- Parse xvars ---- */
    local xvarlist `xvars'
    local nindep : word count `xvarlist'

    if `nindep' < 1 {
        display as error "xvars() must contain at least 1 covariate"
        exit 198
    }

    /* ---- Handle replace ---- */
    if "`replace'" != "" {
        forvalues j = 1/`n_regressors' {
            capture drop `generate'_`j'
            if `do_est_var' {
                capture drop `generate'_`j'_var
            }
        }
    }

    /* ---- Mark sample ---- */
    if `allow_missing_x' {
        marksample touse, novarlist
        markout `touse' `depvar'
        foreach wv of local regvars {
            markout `touse' `wv'
        }
        if "`cluster_var'" != "" {
            markout `touse' `cluster_var'
        }
        if "`weight_var'" != "" {
            markout `touse' `weight_var'
        }
    }
    else {
        marksample touse
        markout `touse' `xvarlist'
    }
    quietly count if `touse'
    local n_use = r(N)

    if `n_use' < 2 {
        display as error "need at least 2 non-missing observations"
        exit 2000
    }

    /* ---- Parse equalize cluster weights ---- */
    if "`equalizeclusterweights'" != "" {
        if "`cluster'" == "" {
            display as error "equalizeclusterweights requires cluster() option"
            exit 198
        }
        /* Compute 1/cluster_size for each observation */
        tempvar _eq_clsize _eq_wt
        quietly bysort `cluster': gen long `_eq_clsize' = _N if `touse'
        quietly gen double `_eq_wt' = 1.0 / `_eq_clsize' if `touse'
        /* Combine with existing weights if any */
        if "`weight_var'" != "" {
            quietly replace `_eq_wt' = `_eq_wt' * `weight_var' if `touse'
        }
        /* Use equalized weights as the weight variable */
        local weight_var `_eq_wt'
    }

    /* ---- Create output variables ---- */
    /* One coefficient estimate per regressor */
    local n_output `n_regressors'
    if `do_est_var' {
        local n_output = `n_regressors' * 2
    }

    forvalues j = 1/`n_regressors' {
        confirm new variable `generate'_`j'
        quietly gen double `generate'_`j' = .
        if `do_est_var' {
            confirm new variable `generate'_`j'_var
            quietly gen double `generate'_`j'_var = .
        }
    }

        /* ---- Inline tuning ---- */
    if `"`tuneparameters'"' != "" {
        display as text ""
        display as text "Running inline parameter tuning..."
        display as text "  Parameters: `tuneparameters'"
        display as text "  Tune trees: `tunenumtrees'  Tune reps: `tunenumreps'"

        /* Call grf_tune to find best parameters */
        grf_tune `varlist' if `touse', foresttype(lm_forest) ///
            numreps(`tunenumreps') tunetrees(`tunenumtrees') seed(`seed') ///
            numthreads(`numthreads') xvars(`xvars')

        /* Override specified parameters with tuned values */
        foreach _tp of local tuneparameters {
            if "`_tp'" == "mtry" {
                local mtry = r(best_mtry)
                display as text "  Tuned mtry: `mtry'"
            }
            else if "`_tp'" == "minnodesize" {
                local minnodesize = r(best_min_node_size)
                display as text "  Tuned min_node_size: `minnodesize'"
            }
            else if "`_tp'" == "samplefrac" {
                local samplefrac = r(best_sample_fraction)
                display as text "  Tuned sample_fraction: `samplefrac'"
            }
            else if "`_tp'" == "honestyfrac" {
                local honestyfrac = r(best_honesty_fraction)
                display as text "  Tuned honesty_fraction: `honestyfrac'"
            }
            else if "`_tp'" == "alpha" {
                local alpha = r(best_alpha)
                display as text "  Tuned alpha: `alpha'"
            }
            else if "`_tp'" == "imbalancepenalty" {
                local imbalancepenalty = r(best_imbalance_penalty)
                display as text "  Tuned imbalance_penalty: `imbalancepenalty'"
            }
        }
        display as text ""
    }

    /* ---- Display header ---- */
    display as text ""
    display as text "Generalized Random Forest: Linear Model Forest"
    display as text "{hline 60}"
    display as text "Outcome variable:      " as result "`depvar'"
    display as text "Regressor variables:   " as result "`regvars'"
    display as text "  Number of W:         " as result `n_regressors'
    display as text "Splitting covariates:  " as result "`xvarlist'"
    display as text "Observations:          " as result `n_use'
    display as text "Trees:                 " as result `ntrees'
    display as text "Nuisance trees:        " as result `nuisancetrees'
    display as text "Honesty:               " as result cond(`do_honesty', "yes", "no")
    display as text "Stabilize splits:      " as result cond(`do_stabilize', "yes", "no")
    if `do_est_var' {
        display as text "Variance estimation:   " as result "yes (ci_group_size=`cigroupsize')"
    }
    display as text "{hline 60}"
    display as text ""

    /* ---- Load plugin ---- */
    if ( inlist("`c(os)'", "MacOSX") | strpos("`c(machine_type)'", "Mac") ) local c_os_ macosx
    else local c_os_: di lower("`c(os)'")

    capture program grf_plugin, plugin using("grf_plugin_`c_os_'.plugin")

    /* ==== Nuisance estimation pipeline ====
     *
     * Linear model forests require orthogonalization:
     *   1. Fit Y ~ X to get Y.hat (outcome model)
     *   2. For each W_k, fit W_k ~ X to get W_k.hat
     *   3. Center: Y.c = Y - Y.hat, W_k.c = W_k - W_k.hat
     *   4. Fit lm_forest on (X, Y.c, W1.c, ..., Wk.c)
     */

    local total_steps = `n_regressors' + 2

    /* Nuisance column indices (regression: X + Y = nindep + 1 data cols) */
    local _nuis_data_cols = `nindep' + 1
    local _nuis_cluster_idx 0
    local _nuis_weight_idx 0
    if "`cluster_var'" != "" {
        local _nuis_cluster_idx = `_nuis_data_cols' + 1
    }
    if "`weight_var'" != "" {
        local _nuis_offset = 0
        if "`cluster_var'" != "" {
            local _nuis_offset = 1
        }
        local _nuis_weight_idx = `_nuis_data_cols' + `_nuis_offset' + 1
    }

    /* Build extra_vars for nuisance calls */
    local _nuis_extra_vars ""
    if "`cluster_var'" != "" {
        local _nuis_extra_vars `_nuis_extra_vars' `cluster_var'
    }
    if "`weight_var'" != "" {
        local _nuis_extra_vars `_nuis_extra_vars' `weight_var'
    }

    if "`yhatinput'" != "" & "`whatinput'" != "" {
        /* ---- User-supplied nuisance estimates ---- */
        local n_whatinput : word count `whatinput'
        if `n_whatinput' != `n_regressors' {
            display as error "whatinput() must contain exactly `n_regressors' variables"
            exit 198
        }
        display as text "Using user-supplied nuisance estimates"
        tempvar yhat
        quietly gen double `yhat' = `yhatinput' if `touse'
        local w_centered_vars ""
        forvalues j = 1/`n_regressors' {
            local wv : word `j' of `regvars'
            local wiv : word `j' of `whatinput'
            tempvar what_`j'
            tempvar wc_`j'
            quietly gen double `what_`j'' = `wiv' if `touse'
            quietly gen double `wc_`j'' = `wv' - `what_`j'' if `touse'
            local w_centered_vars `w_centered_vars' `wc_`j''
        }
        tempvar y_centered
        quietly gen double `y_centered' = `depvar' - `yhat' if `touse'
    }
    else if "`yhatinput'" != "" | "`whatinput'" != "" {
        display as error "both yhatinput() and whatinput() must be specified together"
        exit 198
    }
    else {
        /* ---- Step 1: Fit Y.hat ---- */
        display as text "Step 1/`total_steps': Fitting nuisance model Y ~ X ..."
        tempvar yhat
        quietly gen double `yhat' = .
        plugin call grf_plugin `xvarlist' `depvar' `_nuis_extra_vars' `yhat' ///
            if `touse',                                    ///
            "regression"                                   ///
            "`nuisancetrees'"                              ///
            "`seed'"                                       ///
            "`mtry'"                                       ///
            "`minnodesize'"                                ///
            "`samplefrac'"                                 ///
            "`do_honesty'"                                 ///
            "`honestyfrac'"                                ///
            "`do_honesty_prune'"                           ///
            "`alpha'"                                      ///
            "`imbalancepenalty'"                            ///
            "`cigroupsize'"                                ///
            "`numthreads'"                                 ///
            "0"                                            ///
            "0"                                            ///
            "`nindep'"                                     ///
            "1"                                            ///
            "0"                                            ///
            "0"                                            ///
            "1"                                            ///
            "`allow_missing_x'"                            ///
            "`_nuis_cluster_idx'"                          ///
            "`_nuis_weight_idx'"

        /* ---- Step 2: Fit W_k.hat for each regressor ---- */
        local w_centered_vars ""
        forvalues j = 1/`n_regressors' {
            local step = `j' + 1
            local wv : word `j' of `regvars'
            display as text "Step `step'/`total_steps': Fitting nuisance model `wv' ~ X ..."

            tempvar what_`j'
            quietly gen double `what_`j'' = .
            plugin call grf_plugin `xvarlist' `wv' `_nuis_extra_vars' `what_`j'' ///
                if `touse',                                     ///
                "regression"                                    ///
                "`nuisancetrees'"                               ///
                "`seed'"                                        ///
                "`mtry'"                                        ///
                "`minnodesize'"                                 ///
                "`samplefrac'"                                  ///
                "`do_honesty'"                                  ///
                "`honestyfrac'"                                 ///
                "`do_honesty_prune'"                            ///
                "`alpha'"                                       ///
                "`imbalancepenalty'"                             ///
                "`cigroupsize'"                                 ///
                "`numthreads'"                                  ///
                "0"                                             ///
                "0"                                             ///
                "`nindep'"                                      ///
                "1"                                             ///
                "0"                                             ///
                "0"                                             ///
                "1"                                             ///
                "`allow_missing_x'"                             ///
                "`_nuis_cluster_idx'"                           ///
                "`_nuis_weight_idx'"

            tempvar wc_`j'
            quietly gen double `wc_`j'' = `wv' - `what_`j'' if `touse'
            local w_centered_vars `w_centered_vars' `wc_`j''
        }

        /* ---- Center Y ---- */
        tempvar y_centered
        quietly gen double `y_centered' = `depvar' - `yhat' if `touse'
    }

    /* ---- Save nuisance estimates ---- */
    capture drop _grf_lm_yhat
    quietly gen double _grf_lm_yhat = `yhat'
    label variable _grf_lm_yhat "Y.hat from nuisance regression (Y ~ X)"

    forvalues j = 1/`n_regressors' {
        local wv : word `j' of `regvars'
        capture drop _grf_lm_what`j'
        quietly gen double _grf_lm_what`j' = `what_`j''
        label variable _grf_lm_what`j' "W`j'.hat from nuisance regression (`wv' ~ X)"
    }

    /* ---- Save user-specified generate variables ---- */
    if "`yhatgenerate'" != "" {
        if "`replace'" != "" {
            capture drop `yhatgenerate'
        }
        confirm new variable `yhatgenerate'
        quietly gen double `yhatgenerate' = `yhat'
        label variable `yhatgenerate' "Y.hat from nuisance regression (Y ~ X)"
    }
    if "`whatgenerate'" != "" {
        local n_wg : word count `whatgenerate'
        if `n_wg' != `n_regressors' {
            display as error "whatgenerate() must contain exactly `n_regressors' names"
            exit 198
        }
        forvalues j = 1/`n_regressors' {
            local wgn : word `j' of `whatgenerate'
            if "`replace'" != "" {
                capture drop `wgn'
            }
            confirm new variable `wgn'
            quietly gen double `wgn' = `what_`j''
            local wv : word `j' of `regvars'
            label variable `wgn' "W`j'.hat from nuisance regression (`wv' ~ X)"
        }
    }

    /* ---- Build output varlist ----
     * Plugin writes contiguously: all coefficient predictions first, then all variances.
     * So output_vars must be: coef_1 coef_2 ... var_1 var_2 ...
     */
    local output_vars ""
    forvalues j = 1/`n_regressors' {
        local output_vars `output_vars' `generate'_`j'
    }
    if `do_est_var' {
        forvalues j = 1/`n_regressors' {
            local output_vars `output_vars' `generate'_`j'_var
        }
    }

    /* ---- Build extra vars for cluster/weight ---- */
    local extra_vars ""
    local n_data_before = `nindep' + 1 + `n_regressors'
    if "`cluster_var'" != "" {
        local extra_vars `extra_vars' `cluster_var'
        local cluster_col_idx = `n_data_before' + 1
        local n_data_before = `n_data_before' + 1
    }
    if "`weight_var'" != "" {
        local extra_vars `extra_vars' `weight_var'
        local weight_col_idx = `n_data_before' + 1
    }

    /* ---- Call plugin for linear model forest ----
     *
     * Variable order: X1..Xp Y.c W1.c W2.c ... [cluster] [weight] out_1 [out_1_var] out_2 [out_2_var] ...
     *   n_x = nindep
     *   n_y = 1 (centered outcome)
     *   n_w = n_regressors (centered regressor variables)
     *   n_z = 0
     *   n_output = n_regressors (or n_regressors*2 with variance)
     *
     * argv[20-22]: allow_missing_x, cluster_col_idx, weight_col_idx
     * argv[23]: stabilize_splits (0 by default for lm_forest)
     */
    local final_step = `n_regressors' + 2
    display as text "Step `final_step'/`final_step': Fitting linear model forest on centered data ..."
    plugin call grf_plugin `xvarlist' `y_centered' `w_centered_vars' ///
        `extra_vars' `output_vars'                                    ///
        if `touse',                                                   ///
        "lm_forest"                                                  ///
        "`ntrees'"                                                   ///
        "`seed'"                                                     ///
        "`mtry'"                                                     ///
        "`minnodesize'"                                              ///
        "`samplefrac'"                                               ///
        "`do_honesty'"                                               ///
        "`honestyfrac'"                                              ///
        "`do_honesty_prune'"                                         ///
        "`alpha'"                                                    ///
        "`imbalancepenalty'"                                          ///
        "`cigroupsize'"                                              ///
        "`numthreads'"                                               ///
        "`do_est_var'"                                               ///
        "0"                                                          ///
        "`nindep'"                                                   ///
        "1"                                                          ///
        "`n_regressors'"                                             ///
        "0"                                                          ///
        "`n_output'"                                                 ///
        "`allow_missing_x'"                                          ///
        "`cluster_col_idx'"                                          ///
        "`weight_col_idx'"                                           ///
        "`do_stabilize'"                                                 ///
        "`gradient_weights_str'"

    /* ---- Store results ---- */
    ereturn clear
    ereturn scalar N           = `n_use'
    ereturn scalar n_trees     = `ntrees'
    ereturn scalar seed        = `seed'
    ereturn scalar mtry        = `mtry'
    ereturn scalar min_node    = `minnodesize'
    ereturn scalar alpha       = `alpha'
    ereturn scalar honesty     = `do_honesty'
    ereturn scalar honesty_prune = `do_honesty_prune'
    ereturn scalar sample_fraction    = `samplefrac'
    ereturn scalar honesty_fraction   = `honestyfrac'
    ereturn scalar imbalance_penalty  = `imbalancepenalty'
    ereturn scalar ci_group_size      = `cigroupsize'
    ereturn scalar stabilize   = `do_stabilize'
    ereturn scalar n_regressors = `n_regressors'
    ereturn local  cmd           "grf_lm_forest"
    ereturn scalar allow_missing_x = `allow_missing_x'
    ereturn local  forest_type   "lm_forest"
    ereturn local  depvar        "`depvar'"
    ereturn local  regvars       "`regvars'"
    ereturn local  indepvars     "`xvarlist'"
    ereturn local  predict_stub  "`generate'"
    ereturn local  yhat_var      "_grf_lm_yhat"
    if "`cluster_var'" != "" {
        ereturn local cluster_var "`cluster_var'"
    }
    if "`weight_var'" != "" {
        ereturn local weight_var "`weight_var'"
    }

    /* Per-regressor results */
    forvalues j = 1/`n_regressors' {
        local wv : word `j' of `regvars'
        quietly summarize `generate'_`j' if `touse'
        ereturn scalar coef_mean_`j' = r(mean)
        ereturn scalar coef_sd_`j'   = r(sd)
        ereturn local  what_var_`j'    "_grf_lm_what`j'"
        ereturn local  coef_var_`j'    "`generate'_`j'"
    }

    /* ---- Summary stats ---- */
    display as text ""
    display as text "Linear Model Forest Results"
    display as text "{hline 60}"

    forvalues j = 1/`n_regressors' {
        local wv : word `j' of `regvars'
        quietly summarize `generate'_`j' if `touse'
        local n_pred = r(N)
        local pred_mean = r(mean)
        local pred_sd = r(sd)

        display as text ""
        display as text "Coefficient h_`j'(x) for `wv':"
        display as text "  Variable:     " as result "`generate'_`j'"
        display as text "  Non-missing:  " as result `n_pred'
        display as text "  Mean:         " as result %9.4f `pred_mean'
        display as text "  SD:           " as result %9.4f `pred_sd'
        if `do_est_var' {
            quietly summarize `generate'_`j'_var if `touse'
            display as text "  Mean variance: " as result %9.6f r(mean)
        }
    }

    display as text ""
    display as text "{hline 60}"
    display as text ""
end
