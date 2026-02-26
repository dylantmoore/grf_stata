*! grf_multi_arm_causal_forest.ado -- Multi-Arm Causal Forest via grf C++ library
*! Version 0.1.0
*! Implements Nie and Wager (2021) multi_arm_causal_forest()
*!
*! Estimates heterogeneous treatment effects with multiple discrete treatment arms.
*! Each treatment arm gets its own CATE estimate relative to the control (first level).

program define grf_multi_arm_causal_forest, eclass
    version 14.0

    /* ---- Syntax ----
     * varlist: depvar treatvar1 [treatvar2 ...] indepvar1 [indepvar2 ...]
     * We need to know where treatment vars end and indep vars begin.
     * Use ntreat() to specify how many of the front variables are treatment indicators.
     */
    syntax varlist(min=3 numeric) [if] [in],  ///
        GENerate(name)                         ///
        NTReat(integer)                        ///
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
            noSTABilizesplits                  ///
            REPlace                            ///
        ]

    /* ---- Validate ntreat ---- */
    if `ntreat' < 1 {
        display as error "ntreat() must be at least 1"
        exit 198
    }

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

    /* ---- Parse varlist: depvar treat1 treat2 ... indepvar1 indepvar2 ... ---- */
    gettoken depvar rest : varlist
    local treatvars ""
    forvalues j = 1/`ntreat' {
        gettoken tv rest : rest
        local treatvars `treatvars' `tv'
    }
    local indepvars `rest'
    local nindep : word count `indepvars'

    if `nindep' < 1 {
        display as error "need at least 1 predictor variable"
        exit 198
    }

    /* ---- Handle replace ---- */
    if "`replace'" != "" {
        forvalues j = 1/`ntreat' {
            capture drop `generate'_t`j'
            if `do_est_var' {
                capture drop `generate'_t`j'_var
            }
        }
    }

    /* ---- Mark sample ---- */
    marksample touse
    quietly count if `touse'
    local n_use = r(N)

    if `n_use' < 2 {
        display as error "need at least 2 non-missing observations"
        exit 2000
    }

    /* ---- Validate treatment variables ---- */
    foreach tv of local treatvars {
        quietly summarize `tv' if `touse'
        if r(min) == r(max) {
            display as error "treatment variable `tv' has no variation"
            exit 198
        }
    }

    /* ---- Create output variables ---- */
    /* One CATE estimate per treatment arm */
    local n_output `ntreat'
    if `do_est_var' {
        local n_output = `ntreat' * 2
    }

    forvalues j = 1/`ntreat' {
        confirm new variable `generate'_t`j'
        quietly gen double `generate'_t`j' = .
        if `do_est_var' {
            confirm new variable `generate'_t`j'_var
            quietly gen double `generate'_t`j'_var = .
        }
    }

    /* ---- Build output varlist ----
     * Plugin writes contiguously: all predictions first, then all variances.
     * So output_vars must be: pred_t1 pred_t2 ... var_t1 var_t2 ...
     */
    local output_vars ""
    forvalues j = 1/`ntreat' {
        local output_vars `output_vars' `generate'_t`j'
    }
    if `do_est_var' {
        forvalues j = 1/`ntreat' {
            local output_vars `output_vars' `generate'_t`j'_var
        }
    }

    /* ---- Display header ---- */
    display as text ""
    display as text "Generalized Random Forest: Multi-Arm Causal Forest"
    display as text "{hline 60}"
    display as text "Outcome variable:      " as result "`depvar'"
    display as text "Treatment arms:        " as result "`treatvars'"
    display as text "  Number of arms:      " as result `ntreat'
    display as text "Predictors:            " as result "`indepvars'"
    display as text "Observations:          " as result `n_use'
    display as text "Trees:                 " as result `ntrees'
    display as text "Honesty:               " as result cond(`do_honesty', "yes", "no")
    display as text "Stabilize splits:      " as result cond(`do_stabilize', "yes", "no")
    if `do_est_var' {
        display as text "Variance estimation:   " as result "yes (ci_group_size=`cigroupsize')"
    }
    display as text "{hline 60}"
    display as text ""

    /* ---- Load plugin ---- */
    capture program grf_plugin, plugin ///
        using("grf_plugin.darwin-arm64.plugin")
    if _rc & _rc != 110 {
        capture program grf_plugin, plugin ///
            using("grf_plugin.darwin-x86_64.plugin")
        if _rc & _rc != 110 {
            capture program grf_plugin, plugin ///
                using("grf_plugin.linux-x86_64.plugin")
            if _rc & _rc != 110 {
                capture program grf_plugin, plugin ///
                    using("grf_plugin.windows-x86_64.plugin")
                if _rc & _rc != 110 {
                    display as error "could not load grf_plugin"
                    display as error "make sure the .plugin file is installed"
                    exit 601
                }
            }
        }
    }

    /* ==== Nuisance estimation pipeline ====
     *
     * Multi-arm causal forests require orthogonalization:
     *   1. Fit Y ~ X to get Y.hat (outcome model)
     *   2. For each treatment arm k, fit W_k ~ X to get W_k.hat (propensity)
     *   3. Center: Y.c = Y - Y.hat, W_k.c = W_k - W_k.hat
     *   4. Fit multi-arm causal forest on (X, Y.c, W1.c, W2.c, ...)
     */

    /* ---- Step 1: Fit Y.hat ---- */
    display as text "Step 1/`=`ntreat'+2': Fitting nuisance model Y ~ X ..."
    tempvar yhat
    quietly gen double `yhat' = .
    plugin call grf_plugin `indepvars' `depvar' `yhat' ///
        if `touse',                                     ///
        "regression"                                    ///
        "`ntrees'"                                      ///
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
        "1"

    /* ---- Step 2: Fit W_k.hat for each treatment arm ---- */
    local w_centered_vars ""
    forvalues j = 1/`ntreat' {
        local step = `j' + 1
        local tv : word `j' of `treatvars'
        display as text "Step `step'/`=`ntreat'+2': Fitting nuisance model `tv' ~ X ..."

        tempvar what_`j'
        quietly gen double `what_`j'' = .
        plugin call grf_plugin `indepvars' `tv' `what_`j'' ///
            if `touse',                                      ///
            "regression"                                     ///
            "`ntrees'"                                       ///
            "`seed'"                                         ///
            "`mtry'"                                         ///
            "`minnodesize'"                                  ///
            "`samplefrac'"                                   ///
            "`do_honesty'"                                   ///
            "`honestyfrac'"                                  ///
            "`do_honesty_prune'"                             ///
            "`alpha'"                                        ///
            "`imbalancepenalty'"                               ///
            "`cigroupsize'"                                  ///
            "`numthreads'"                                   ///
            "0"                                              ///
            "0"                                              ///
            "`nindep'"                                       ///
            "1"                                              ///
            "0"                                              ///
            "0"                                              ///
            "1"

        tempvar wc_`j'
        quietly gen double `wc_`j'' = `tv' - `what_`j'' if `touse'
        local w_centered_vars `w_centered_vars' `wc_`j''
    }

    /* ---- Center Y ---- */
    tempvar y_centered
    quietly gen double `y_centered' = `depvar' - `yhat' if `touse'

    /* ---- Save nuisance estimates ---- */
    capture drop _grf_mac_yhat
    quietly gen double _grf_mac_yhat = `yhat'
    label variable _grf_mac_yhat "Y.hat from nuisance regression (Y ~ X)"

    forvalues j = 1/`ntreat' {
        local tv : word `j' of `treatvars'
        capture drop _grf_mac_what`j'
        quietly gen double _grf_mac_what`j' = `what_`j''
        label variable _grf_mac_what`j' "W`j'.hat from nuisance regression (`tv' ~ X)"
    }

    /* ---- Call plugin for multi-arm causal forest ----
     *
     * Variable order: X1..Xp Y.c W1.c W2.c ... out_t1 [out_t1_var] out_t2 [out_t2_var] ...
     *   n_x = nindep
     *   n_y = 1 (centered outcome)
     *   n_w = ntreat (centered treatment indicators)
     *   n_z = 0
     *   n_output = ntreat (or ntreat*2 with variance)
     *
     * argv[20]: stabilize_splits
     * argv[21]: num_treatments
     */
    local final_step = `ntreat' + 2
    display as text "Step `final_step'/`final_step': Fitting multi-arm causal forest ..."
    plugin call grf_plugin `indepvars' `y_centered' `w_centered_vars' `output_vars' ///
        if `touse',                                                                  ///
        "multi_arm_causal"                                                           ///
        "`ntrees'"                                                                   ///
        "`seed'"                                                                     ///
        "`mtry'"                                                                     ///
        "`minnodesize'"                                                              ///
        "`samplefrac'"                                                               ///
        "`do_honesty'"                                                               ///
        "`honestyfrac'"                                                              ///
        "`do_honesty_prune'"                                                         ///
        "`alpha'"                                                                    ///
        "`imbalancepenalty'"                                                          ///
        "`cigroupsize'"                                                              ///
        "`numthreads'"                                                               ///
        "`do_est_var'"                                                               ///
        "0"                                                                          ///
        "`nindep'"                                                                   ///
        "1"                                                                          ///
        "`ntreat'"                                                                   ///
        "0"                                                                          ///
        "`n_output'"                                                                 ///
        "`do_stabilize'"                                                             ///
        "`ntreat'"

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
    ereturn scalar n_treat     = `ntreat'
    ereturn local  cmd           "grf_multi_arm_causal_forest"
    ereturn local  forest_type   "multi_causal"
    ereturn local  depvar        "`depvar'"
    ereturn local  treatvars     "`treatvars'"
    ereturn local  indepvars     "`indepvars'"
    ereturn local  predict_stub  "`generate'"
    ereturn local  yhat_var      "_grf_mac_yhat"

    /* Per-arm ATE estimates */
    forvalues j = 1/`ntreat' {
        local tv : word `j' of `treatvars'
        quietly summarize `generate'_t`j' if `touse'
        ereturn scalar ate_`j'    = r(mean)
        ereturn scalar ate_se_`j' = r(sd) / sqrt(r(N))
        ereturn local  what_var_`j' "_grf_mac_what`j'"
    }

    /* ---- Summary stats ---- */
    display as text ""
    display as text "Multi-Arm Causal Forest Results"
    display as text "{hline 60}"

    forvalues j = 1/`ntreat' {
        local tv : word `j' of `treatvars'
        quietly summarize `generate'_t`j' if `touse'
        local n_pred = r(N)
        local pred_mean = r(mean)
        local pred_sd = r(sd)
        local ate_se_j = r(sd) / sqrt(r(N))

        display as text ""
        display as text "Treatment arm `j' (`tv'):"
        display as text "  Variable:     " as result "`generate'_t`j'"
        display as text "  Non-missing:  " as result `n_pred'
        display as text "  Mean (ATE):   " as result %9.4f `pred_mean'
        display as text "  SD:           " as result %9.4f `pred_sd'
        display as text "  ATE s.e.:     " as result %9.4f `ate_se_j'
        if `do_est_var' {
            quietly summarize `generate'_t`j'_var if `touse'
            display as text "  Mean variance: " as result %9.6f r(mean)
        }
    }

    display as text ""
    display as text "{hline 60}"
    display as text ""
end
