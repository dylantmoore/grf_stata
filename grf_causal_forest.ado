*! grf_causal_forest.ado -- Causal Forest via grf C++ library
*! Version 0.1.0
*! Implements Athey, Tibshirani, Wager (2019) causal_forest()

program define grf_causal_forest, eclass
    version 14.0

    syntax varlist(min=3 numeric) [if] [in],  ///
        GENerate(name)                         ///
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
            ALPha(real 0.05)                   ///
            IMBalancepenalty(real 0.0)          ///
            CIGroupsize(integer 1)             ///
            NUMThreads(integer 0)              ///
            ESTIMATEVariance                   ///
            noSTABilizesplits                  ///
            YHATGenerate(name)                 ///
            WHATGenerate(name)                 ///
            YHATinput(varname numeric)         ///
            WHATinput(varname numeric)         ///
            REPlace                            ///
            VARGenerate(name)                  ///
            noMIA                              ///
            CLuster(varname numeric)           ///
            WEIghts(varname numeric)           ///
            EQUALizeclusterweights             ///
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
    local do_stabilize 1
    if "`stabilizesplits'" == "nostabilizesplits" {
        local do_stabilize 0
    }

    /* ---- Parse estimate_variance ---- */
    local do_est_var 0
    if "`estimatevariance'" != "" {
        local do_est_var 1
        if `cigroupsize' < 2 {
            local cigroupsize 2
        }
    }

    /* ---- Handle replace ---- */
    if "`replace'" != "" {
        capture drop `generate'
        if `do_est_var' & "`vargenerate'" != "" {
            capture drop `vargenerate'
        }
        if "`yhatgenerate'" != "" {
            capture drop `yhatgenerate'
        }
        if "`whatgenerate'" != "" {
            capture drop `whatgenerate'
        }
    }
    confirm new variable `generate'

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

    /* ---- Parse varlist: depvar treatvar indepvars ---- */
    gettoken depvar rest : varlist
    gettoken treatvar indepvars : rest
    local nindep : word count `indepvars'

    if `nindep' < 1 {
        display as error "need at least 1 predictor variable"
        exit 198
    }

    /* ---- Mark sample ---- */
    if `allow_missing_x' {
        marksample touse, novarlist
        markout `touse' `depvar'
        markout `touse' `treatvar'
        if "`cluster_var'" != "" {
            markout `touse' `cluster_var'
        }
        if "`weight_var'" != "" {
            markout `touse' `weight_var'
        }
    }
    else {
        marksample touse
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

    /* ---- Validate treatment variable ---- */
    quietly summarize `treatvar' if `touse'
    local w_min = r(min)
    local w_max = r(max)
    if `w_min' == `w_max' {
        display as error "treatment variable `treatvar' has no variation"
        exit 198
    }

    /* ---- Display header ---- */
    display as text ""
    display as text "Generalized Random Forest: Causal Forest"
    display as text "{hline 55}"
    display as text "Outcome variable:      " as result "`depvar'"
    display as text "Treatment variable:    " as result "`treatvar'"
    display as text "Predictors:            " as result "`indepvars'"
    display as text "Observations:          " as result `n_use'
    display as text "Trees:                 " as result `ntrees'
    display as text "Honesty:               " as result cond(`do_honesty', "yes", "no")
    display as text "Stabilize splits:      " as result cond(`do_stabilize', "yes", "no")
    if `do_est_var' {
        display as text "Variance estimation:   " as result "yes (ci_group_size=`cigroupsize')"
    }
    display as text "{hline 55}"
    display as text ""

    /* ---- Load plugin ---- */
    if ( inlist("`c(os)'", "MacOSX") | strpos("`c(machine_type)'", "Mac") ) local c_os_ macosx
    else local c_os_: di lower("`c(os)'")

    capture program grf_plugin, plugin using("grf_plugin_`c_os_'.plugin")

    /* ==== Nuisance estimation pipeline ====
     *
     * Causal forests require orthogonalization (Robinson, 1988):
     *   1. Fit Y ~ X to get Y.hat (outcome model)
     *   2. Fit W ~ X to get W.hat (propensity model)
     *   3. Center: Y.centered = Y - Y.hat, W.centered = W - W.hat
     *   4. Fit causal forest on (X, Y.centered, W.centered)
     */

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

    /* Check user-supplied nuisance estimates */
    if "`yhatinput'" != "" & "`whatinput'" != "" {
        /* Use user-supplied nuisance estimates */
        display as text "Steps 1-2: Using user-supplied nuisance estimates"
        display as text "  Y.hat from: `yhatinput'"
        display as text "  W.hat from: `whatinput'"
        tempvar yhat what
        quietly gen double `yhat' = `yhatinput' if `touse'
        quietly gen double `what' = `whatinput' if `touse'
    }
    else if "`yhatinput'" != "" | "`whatinput'" != "" {
        display as error "both yhatinput() and whatinput() must be specified together"
        exit 198
    }
    else {
        display as text "Step 1/3: Fitting nuisance model Y ~ X ..."
        tempvar yhat
        quietly gen double `yhat' = .
        plugin call grf_plugin `indepvars' `depvar' `extra_vars' `yhat' ///
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

        display as text "Step 2/3: Fitting nuisance model W ~ X ..."
        tempvar what
        quietly gen double `what' = .
        plugin call grf_plugin `indepvars' `treatvar' `extra_vars' `what' ///
            if `touse',                                      ///
            "regression"                                     ///
            "`nuisancetrees'"                                ///
            "`seed'"                                         ///
            "`mtry'"                                         ///
            "`minnodesize'"                                  ///
            "`samplefrac'"                                   ///
            "`do_honesty'"                                   ///
            "`honestyfrac'"                                  ///
            "`do_honesty_prune'"                             ///
            "`alpha'"                                        ///
            "`imbalancepenalty'"                              ///
            "`cigroupsize'"                                  ///
            "`numthreads'"                                   ///
            "0"                                              ///
            "0"                                              ///
            "`nindep'"                                       ///
            "1"                                              ///
            "0"                                              ///
            "0"                                              ///
            "1"                                              ///
            "`allow_missing_x'"                              ///
            "`_nuis_cluster_idx'"                            ///
            "`_nuis_weight_idx'"
    }

    /* ---- Center Y and W ---- */
    tempvar y_centered w_centered
    quietly gen double `y_centered' = `depvar' - `yhat' if `touse'
    quietly gen double `w_centered' = `treatvar' - `what' if `touse'

    /* ---- Save nuisance estimates (always, for post-estimation) ---- */
    /* Internal names used by grf_ate and grf_test_calibration */
    capture drop _grf_yhat
    capture drop _grf_what
    quietly gen double _grf_yhat = `yhat'
    quietly gen double _grf_what = `what'
    label variable _grf_yhat "Y.hat from nuisance regression (Y ~ X)"
    label variable _grf_what "W.hat from nuisance regression (W ~ X)"

    /* Optionally save with user-specified names too */
    if "`yhatgenerate'" != "" {
        if "`replace'" != "" {
            capture drop `yhatgenerate'
        }
        confirm new variable `yhatgenerate'
        quietly gen double `yhatgenerate' = `yhat'
        label variable `yhatgenerate' "Y.hat from nuisance regression (Y ~ X)"
    }
    if "`whatgenerate'" != "" {
        if "`replace'" != "" {
            capture drop `whatgenerate'
        }
        confirm new variable `whatgenerate'
        quietly gen double `whatgenerate' = `what'
        label variable `whatgenerate' "W.hat from nuisance regression (W ~ X)"
    }

    /* ---- Create output variable(s) ---- */
    quietly gen double `generate' = .
    local n_output 1

    if `do_est_var' {
        if "`vargenerate'" == "" {
            local vargenerate `generate'_var
        }
        if "`replace'" != "" {
            capture drop `vargenerate'
        }
        confirm new variable `vargenerate'
        quietly gen double `vargenerate' = .
        local n_output 2
    }

    /* ---- Build output varlist ---- */
    local output_vars `generate'
    if `do_est_var' {
        local output_vars `generate' `vargenerate'
    }

    /* ---- Build plugin varlist with optional cluster/weight columns ---- */
    local extra_vars ""
    if "`cluster_var'" != "" {
        local extra_vars `extra_vars' `cluster_var'
    }
    if "`weight_var'" != "" {
        local extra_vars `extra_vars' `weight_var'
    }

    /* Compute cluster/weight column indices in the plugin varlist */
    /* Causal forest data vars: X1..Xp Y.c W.c => nindep + 2 columns */
    local _data_col_count = `nindep' + 2
    if "`cluster_var'" != "" {
        local cluster_col_idx = `_data_col_count' + 1
    }
    if "`weight_var'" != "" {
        local _offset = 0
        if "`cluster_var'" != "" {
            local _offset = 1
        }
        local weight_col_idx = `_data_col_count' + `_offset' + 1
    }

    /* ---- Call plugin for causal forest ----
     *
     * Variable order: X1..Xp Y.centered W.centered [cluster] [weight] out1 [out2]
     * argv: forest_type num_trees seed mtry min_node_size sample_fraction
     *       honesty honesty_fraction honesty_prune alpha imbalance_penalty
     *       ci_group_size num_threads estimate_variance compute_oob
     *       n_x n_y n_w n_z n_output
     *       allow_missing_x cluster_col_idx weight_col_idx
     *       stabilize_splits
     */
    display as text "Step 3/3: Fitting causal forest on centered data ..."
    plugin call grf_plugin `indepvars' `y_centered' `w_centered' `extra_vars' `output_vars' ///
        if `touse',                                                             ///
        "causal"                                                                ///
        "`ntrees'"                                                              ///
        "`seed'"                                                                ///
        "`mtry'"                                                                ///
        "`minnodesize'"                                                         ///
        "`samplefrac'"                                                          ///
        "`do_honesty'"                                                          ///
        "`honestyfrac'"                                                         ///
        "`do_honesty_prune'"                                                    ///
        "`alpha'"                                                               ///
        "`imbalancepenalty'"                                                     ///
        "`cigroupsize'"                                                         ///
        "`numthreads'"                                                          ///
        "`do_est_var'"                                                          ///
        "0"                                                                     ///
        "`nindep'"                                                              ///
        "1"                                                                     ///
        "1"                                                                     ///
        "0"                                                                     ///
        "`n_output'"                                                            ///
        "`allow_missing_x'"                                                     ///
        "`cluster_col_idx'"                                                     ///
        "`weight_col_idx'"                                                      ///
        "`do_stabilize'"

    /* ---- Compute ATE ---- */
    quietly summarize `generate' if `touse'
    local ate = r(mean)
    local ate_se = r(sd) / sqrt(r(N))

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
    ereturn scalar ate         = `ate'
    ereturn scalar ate_se      = `ate_se'
    ereturn local  cmd           "grf_causal_forest"
    ereturn scalar allow_missing_x = `allow_missing_x'
    ereturn local  forest_type   "causal"
    ereturn local  depvar        "`depvar'"
    ereturn local  treatvar      "`treatvar'"
    ereturn local  indepvars     "`indepvars'"
    ereturn local  predict_var   "`generate'"
    if `do_est_var' {
        ereturn local variance_var "`vargenerate'"
    }
    ereturn local yhat_var   "_grf_yhat"
    ereturn local what_var   "_grf_what"
    if "`cluster_var'" != "" {
        ereturn local cluster_var "`cluster_var'"
    }
    if "`weight_var'" != "" {
        ereturn local weight_var "`weight_var'"
    }

    /* ---- Summary stats ---- */
    quietly summarize `generate' if `touse'
    local n_pred = r(N)
    local pred_mean = r(mean)
    local pred_sd = r(sd)

    display as text ""
    display as text "Causal Forest Results"
    display as text "{hline 55}"
    display as text "CATE predictions:       " as result "`generate'"
    display as text "  Non-missing:  " as result `n_pred'
    display as text "  Mean (ATE):   " as result %9.4f `pred_mean'
    display as text "  SD:           " as result %9.4f `pred_sd'
    display as text "  ATE s.e.:     " as result %9.4f `ate_se'
    if `do_est_var' {
        quietly summarize `vargenerate' if `touse'
        display as text "Variance estimates:     " as result "`vargenerate'"
        display as text "  Mean variance: " as result %9.6f r(mean)
    }
    display as text "{hline 55}"
    display as text ""
end
