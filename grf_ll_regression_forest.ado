*! grf_ll_regression_forest.ado -- Local Linear Regression Forest via grf C++ library
*! Version 0.2.0
*! Implements Friedberg et al. (2021) ll_regression_forest()
*!
*! Like regression forest but with local linear correction at prediction time.
*!
*! Known limitation: R's ll.split.variables accepts a vector of variable
*! indices to restrict which variables are used for splitting. Stata's
*! llsplit option is a boolean toggle (splits on all variables when enabled).

program define grf_ll_regression_forest, eclass
    version 14.0

    syntax varlist(min=2 numeric) [if] [in],  ///
        GENerate(name)                         ///
        [                                      ///
            NTrees(integer 2000)               ///
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
            REPlace                            ///
            VARGenerate(name)                  ///
            LLSplit                            ///
            LLSplitVars(varlist numeric)       ///
            LLLambda(real 0.1)                 ///
            LLWeightpenalty                    ///
            LLCutoff(integer 0)               ///
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

    /* ---- Parse estimate_variance ---- */
    local do_est_var 0
    if "`estimatevariance'" != "" {
        local do_est_var 1
        if `cigroupsize' < 2 {
            local cigroupsize 2
        }
    }

    /* ---- Parse local linear options ---- */
    local enable_ll_split 0
    if "`llsplit'" != "" {
        local enable_ll_split 1
    }

    local ll_weight_penalty 0
    if "`llweightpenalty'" != "" {
        local ll_weight_penalty 1
    }

    local ll_split_cutoff `llcutoff'

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

    /* ---- Handle replace ---- */
    if "`replace'" != "" {
        capture drop `generate'
        if `do_est_var' & "`vargenerate'" != "" {
            capture drop `vargenerate'
        }
    }
    confirm new variable `generate'

    /* ---- Parse varlist ---- */
    gettoken depvar indepvars : varlist
    local nindep : word count `indepvars'

    if `nindep' < 1 {
        display as error "need at least 1 predictor variable"
        exit 198
    }

    /* ---- Parse ll.split.variables ----
     * LLSplitVars(varlist) specifies which X variables are used for LL splitting.
     * Convert variable names to 0-indexed column positions in the X matrix.
     * If LLSplit is specified without LLSplitVars, all variables are used (default). */
    local ll_split_vars_str ""
    if "`llsplitvars'" != "" {
        foreach sv of local llsplitvars {
            local found 0
            local pos 0
            foreach xv of local indepvars {
                if "`sv'" == "`xv'" {
                    local found 1
                    local ll_split_vars_str "`ll_split_vars_str' `pos'"
                }
                local pos = `pos' + 1
            }
            if !`found' {
                display as error "llsplitvars: variable `sv' not found in predictor list"
                exit 198
            }
        }
        local ll_split_vars_str = trim("`ll_split_vars_str'")
        /* Enable LL split if splitvars are specified */
        if "`llsplit'" == "" {
            local enable_ll_split 1
        }
    }

    /* ---- Mark sample ---- */
    if `allow_missing_x' {
        marksample touse, novarlist
        markout `touse' `depvar'
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

    /* ---- Display header ---- */
    display as text ""
    display as text "Generalized Random Forest: Local Linear Regression Forest"
    display as text "{hline 60}"
    display as text "Dependent variable:    " as result "`depvar'"
    display as text "Predictors:            " as result "`indepvars'"
    display as text "Observations:          " as result `n_use'
    display as text "Trees:                 " as result `ntrees'
    display as text "Honesty:               " as result cond(`do_honesty', "yes", "no")
    display as text "LL lambda:             " as result `lllambda'
    display as text "LL split:              " as result cond(`enable_ll_split', "yes", "no")
    display as text "LL weight penalty:     " as result cond(`ll_weight_penalty', "yes", "no")
    display as text "LL split cutoff:       " as result cond(`ll_split_cutoff'==0, "sqrt(n)", "`ll_split_cutoff'")
    if `do_est_var' {
        display as text "Variance estimation:   " as result "yes (ci_group_size=`cigroupsize')"
    }
    display as text "{hline 60}"
    display as text ""

    /* ---- Load plugin ---- */
    if ( inlist("`c(os)'", "MacOSX") | strpos("`c(machine_type)'", "Mac") ) local c_os_ macosx
    else local c_os_: di lower("`c(os)'")

    capture program grf_plugin, plugin using("grf_plugin_`c_os_'.plugin")

    /* ---- Build output varlist ---- */
    local output_vars `generate'
    if `do_est_var' {
        local output_vars `generate' `vargenerate'
    }

    /* ---- Build extra vars for cluster/weight ---- */
    local extra_vars ""
    local n_data_before = `nindep' + 1
    if "`cluster_var'" != "" {
        local extra_vars `extra_vars' `cluster_var'
        local cluster_col_idx = `n_data_before' + 1
        local n_data_before = `n_data_before' + 1
    }
    if "`weight_var'" != "" {
        local extra_vars `extra_vars' `weight_var'
        local weight_col_idx = `n_data_before' + 1
    }

    /* ---- Call plugin ----
     *
     * Variable order: X1..Xp Y [cluster] [weight] out1 [out2]
     * argv: forest_type num_trees seed mtry min_node_size sample_fraction
     *       honesty honesty_fraction honesty_prune alpha imbalance_penalty
     *       ci_group_size num_threads estimate_variance compute_oob
     *       n_x n_y n_w n_z n_output
     *       allow_missing_x cluster_col_idx weight_col_idx
     *       enable_ll_split ll_lambda ll_weight_penalty ll_split_cutoff
     */
    plugin call grf_plugin `indepvars' `depvar' `extra_vars' `output_vars' ///
        if `touse',                                            ///
        "ll_regression"                                        ///
        "`ntrees'"                                             ///
        "`seed'"                                               ///
        "`mtry'"                                               ///
        "`minnodesize'"                                        ///
        "`samplefrac'"                                         ///
        "`do_honesty'"                                         ///
        "`honestyfrac'"                                        ///
        "`do_honesty_prune'"                                   ///
        "`alpha'"                                              ///
        "`imbalancepenalty'"                                    ///
        "`cigroupsize'"                                        ///
        "`numthreads'"                                         ///
        "`do_est_var'"                                         ///
        "0"                                                    ///
        "`nindep'"                                             ///
        "1"                                                    ///
        "0"                                                    ///
        "0"                                                    ///
        "`n_output'"                                           ///
        "`allow_missing_x'"                                    ///
        "`cluster_col_idx'"                                    ///
        "`weight_col_idx'"                                     ///
        "`enable_ll_split'"                                    ///
        "`lllambda'"                                           ///
        "`ll_weight_penalty'"                                  ///
        "`ll_split_cutoff'"                                    ///
        "`ll_split_vars_str'"

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
    ereturn scalar ll_lambda          = `lllambda'
    ereturn scalar ll_weight_penalty  = `ll_weight_penalty'
    ereturn scalar ll_split_cutoff    = `ll_split_cutoff'
    ereturn scalar enable_ll_split    = `enable_ll_split'
    ereturn local  cmd           "grf_ll_regression_forest"
    ereturn scalar allow_missing_x = `allow_missing_x'
    ereturn local  forest_type   "ll_regression"
    ereturn local  depvar        "`depvar'"
    ereturn local  indepvars     "`indepvars'"
    ereturn local  predict_var   "`generate'"
    if `do_est_var' {
        ereturn local variance_var "`vargenerate'"
    }
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
    display as text "Predictions written to: " as result "`generate'"
    display as text "  Non-missing:  " as result `n_pred'
    display as text "  Mean:         " as result %9.4f `pred_mean'
    display as text "  SD:           " as result %9.4f `pred_sd'
    if `do_est_var' {
        quietly summarize `vargenerate' if `touse'
        display as text "Variance estimates:     " as result "`vargenerate'"
        display as text "  Mean variance: " as result %9.6f r(mean)
    }
    display as text ""
end
