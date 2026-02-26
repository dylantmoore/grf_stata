*! grf_ll_regression_forest.ado -- Local Linear Regression Forest via grf C++ library
*! Version 0.1.0
*! Implements Friedberg et al. (2021) ll_regression_forest()
*!
*! Like regression forest but with local linear correction at prediction time.

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
            LLLambda(real 0.1)                 ///
            LLWeightpenalty                    ///
            LLCutoff(integer 0)               ///
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
    marksample touse
    quietly count if `touse'
    local n_use = r(N)

    if `n_use' < 2 {
        display as error "need at least 2 non-missing observations"
        exit 2000
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
    local plugin_loaded 0
    foreach plat in darwin-arm64 darwin-x86_64 linux-x86_64 windows-x86_64 {
        if !`plugin_loaded' {
            capture findfile grf_plugin.`plat'.plugin
            if _rc == 0 {
                capture program grf_plugin, plugin using("`r(fn)'")
                if _rc == 0 | _rc == 110 {
                    local plugin_loaded 1
                }
            }
        }
    }
    if !`plugin_loaded' {
        display as error "could not load grf_plugin"
        display as error "make sure the .plugin file is installed"
        exit 601
    }

    /* ---- Build output varlist ---- */
    local output_vars `generate'
    if `do_est_var' {
        local output_vars `generate' `vargenerate'
    }

    /* ---- Call plugin ----
     *
     * Variable order: X1..Xp Y out1 [out2]
     * argv: forest_type num_trees seed mtry min_node_size sample_fraction
     *       honesty honesty_fraction honesty_prune alpha imbalance_penalty
     *       ci_group_size num_threads estimate_variance compute_oob
     *       n_x n_y n_w n_z n_output
     *       enable_ll_split ll_lambda ll_weight_penalty ll_split_cutoff
     */
    plugin call grf_plugin `indepvars' `depvar' `output_vars' ///
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
        "`enable_ll_split'"                                    ///
        "`lllambda'"                                           ///
        "`ll_weight_penalty'"                                  ///
        "`ll_split_cutoff'"

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
    ereturn local  forest_type   "ll_regression"
    ereturn local  depvar        "`depvar'"
    ereturn local  indepvars     "`indepvars'"
    ereturn local  predict_var   "`generate'"
    if `do_est_var' {
        ereturn local variance_var "`vargenerate'"
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
