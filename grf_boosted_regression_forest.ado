*! grf_boosted_regression_forest.ado -- Boosted Regression Forest via grf C++ library
*! Version 0.1.0
*! Implements boosted_regression_forest()
*!
*! Iteratively fits regression forests to residuals with optional
*! cross-validation-based auto-tuning for the number of boosting steps.

program define grf_boosted_regression_forest, eclass
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
            noSTABilizesplits                  ///
            ALPha(real 0.05)                   ///
            IMBalancepenalty(real 0.0)          ///
            CIGroupsize(integer 1)             ///
            NUMThreads(integer 0)              ///
            ESTIMATEVariance                   ///
            REPlace                            ///
            VARGenerate(name)                  ///
            BOOSTSteps(integer 0)              ///
            BOOSTMaxsteps(integer 5)           ///
            BOOSTErrorreduction(real 0.97)     ///
            BOOSTTreestune(integer 10)         ///
            noMIA                              ///
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

    /* ---- Validate boost options ---- */
    if `booststeps' < 0 {
        display as error "booststeps() must be >= 0 (0 = auto-tune via CV)"
        exit 198
    }
    if `boostmaxsteps' < 1 {
        display as error "boostmaxsteps() must be >= 1"
        exit 198
    }
    if `boosterrorreduction' <= 0 | `boosterrorreduction' >= 1 {
        display as error "boosterrorreduction() must be in (0, 1)"
        exit 198
    }
    if `boosttreestune' < 1 {
        display as error "boosttreestune() must be >= 1"
        exit 198
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

        /* ---- Inline tuning ---- */
    if `"`tuneparameters'"' != "" {
        display as text ""
        display as text "Running inline parameter tuning..."
        display as text "  Parameters: `tuneparameters'"
        display as text "  Tune trees: `tunenumtrees'  Tune reps: `tunenumreps'"

        /* Call grf_tune to find best parameters */
        grf_tune `varlist' if `touse', foresttype(boosted_regression) ///
            numreps(`tunenumreps') tunetrees(`tunenumtrees') seed(`seed') ///
            numthreads(`numthreads')

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
    display as text "Generalized Random Forest: Boosted Regression Forest"
    display as text "{hline 60}"
    display as text "Dependent variable:    " as result "`depvar'"
    display as text "Predictors:            " as result "`indepvars'"
    display as text "Observations:          " as result `n_use'
    display as text "Trees per step:        " as result `ntrees'
    display as text "Honesty:               " as result cond(`do_honesty', "yes", "no")
    display as text "Stabilize splits:      " as result cond(`do_stabilize', "yes", "no")
    if `booststeps' == 0 {
        display as text "Boost steps:           " as result "auto-tune (max `boostmaxsteps')"
        display as text "  Error reduction:     " as result %6.4f `boosterrorreduction'
        display as text "  CV trees:            " as result `boosttreestune'
    }
    else {
        display as text "Boost steps:           " as result `booststeps'
    }
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
     *       boost_steps boost_error_reduction boost_max_steps boost_trees_tune
     *       stabilize_splits
     */
    plugin call grf_plugin `indepvars' `depvar' `extra_vars' `output_vars' ///
        if `touse',                                            ///
        "boosted_regression"                                   ///
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
        "`booststeps'"                                         ///
        "`boosterrorreduction'"                                ///
        "`boostmaxsteps'"                                      ///
        "`boosttreestune'"                                     ///
        "`do_stabilize'"

    /* ---- Read actual boost steps from plugin scalar ---- */
    local actual_boost_steps = _grf_boost_steps

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
    ereturn scalar boost_steps        = `actual_boost_steps'
    ereturn scalar boost_max_steps    = `boostmaxsteps'
    ereturn scalar boost_error_reduction = `boosterrorreduction'
    ereturn local  cmd           "grf_boosted_regression_forest"
    ereturn scalar allow_missing_x = `allow_missing_x'
    ereturn local  forest_type   "boosted_regression"
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
    display as text "Boosted Regression Forest Results"
    display as text "{hline 55}"
    display as text "Boost steps used:       " as result `actual_boost_steps'
    display as text "Predictions written to: " as result "`generate'"
    display as text "  Non-missing:  " as result `n_pred'
    display as text "  Mean:         " as result %9.4f `pred_mean'
    display as text "  SD:           " as result %9.4f `pred_sd'
    if `do_est_var' {
        quietly summarize `vargenerate' if `touse'
        display as text "Variance estimates:     " as result "`vargenerate'"
        display as text "  Mean variance: " as result %9.6f r(mean)
    }
    display as text "{hline 55}"
    display as text ""
end
