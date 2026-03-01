*! grf_instrumental_forest.ado -- Instrumental Forest via grf C++ library
*! Version 0.1.0
*! Implements Athey, Tibshirani, Wager (2019) instrumental_forest()

program define grf_instrumental_forest, eclass
    version 14.0

    syntax varlist(min=4 numeric) [if] [in],  ///
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
            REDucedformweight(real 0.0)        ///
            noSTABilizesplits                  ///
            NUISancetrees(integer 500)         ///
            YHATinput(varname numeric)         ///
            WHATinput(varname numeric)         ///
            ZHATinput(varname numeric)         ///
            YHATgenerate(name)                 ///
            WHATgenerate(name)                 ///
            ZHATgenerate(name)                 ///
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

    /* ---- Parse estimate_variance ---- */
    local do_est_var 0
    if "`estimatevariance'" != "" {
        local do_est_var 1
        if `cigroupsize' < 2 {
            local cigroupsize 2
        }
    }

    /* ---- Parse stabilize_splits ---- */
    /* R's instrumental_forest defaults stabilize.splits = TRUE.
     * noSTABilizesplits is opt-out (default ON). */
    local do_stabilize 1
    if "`stabilizesplits'" == "nostabilizesplits" {
        local do_stabilize 0
    }

    /* ---- Handle replace ---- */
    if "`replace'" != "" {
        capture drop `generate'
        if `do_est_var' & "`vargenerate'" != "" {
            capture drop `vargenerate'
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

    /* ---- Parse varlist: Y W Z X1..Xp ---- */
    gettoken depvar rest : varlist
    gettoken treatment rest : rest
    gettoken instrument indepvars : rest
    local nindep : word count `indepvars'

    if `nindep' < 1 {
        display as error "need at least 1 predictor variable"
        exit 198
    }

    /* ---- Mark sample ---- */
    if `allow_missing_x' {
        marksample touse, novarlist
        markout `touse' `depvar'
        markout `touse' `treatment'
        markout `touse' `instrument'
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
        grf_tune `varlist' if `touse', foresttype(instrumental) ///
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
    display as text "Generalized Random Forest: Instrumental Forest"
    display as text "{hline 55}"
    display as text "Dependent variable:    " as result "`depvar'"
    display as text "Treatment variable:    " as result "`treatment'"
    display as text "Instrument variable:   " as result "`instrument'"
    display as text "Predictors:            " as result "`indepvars'"
    display as text "Observations:          " as result `n_use'
    display as text "Trees:                 " as result `ntrees'
    display as text "Honesty:               " as result cond(`do_honesty', "yes", "no")
    display as text "Reduced form weight:   " as result `reducedformweight'
    if `do_stabilize' {
        display as text "Stabilize splits:      " as result "yes"
    }
    if `do_est_var' {
        display as text "Variance estimation:   " as result "yes (ci_group_size=`cigroupsize')"
    }
    display as text "{hline 55}"
    display as text ""

    /* ---- Nuisance pipeline: center Y, W, Z using regression forests ----
     *
     * Step 1: Y.hat = E[Y|X] via regression forest
     * Step 2: W.hat = E[W|X] via regression forest
     * Step 3: Z.hat = E[Z|X] via regression forest
     * Step 4: Y.c = Y - Y.hat, W.c = W - W.hat, Z.c = Z - Z.hat
     */

    /* ---- Load plugin ---- */
    if ( inlist("`c(os)'", "MacOSX") | strpos("`c(machine_type)'", "Mac") ) local c_os_ macosx
    else local c_os_: di lower("`c(os)'")

    capture program grf_plugin, plugin using("grf_plugin_`c_os_'.plugin")

    if "`yhatinput'" != "" & "`whatinput'" != "" & "`zhatinput'" != "" {
        display as text "Nuisance: Using user-supplied estimates"
        tempvar Y_hat W_hat Z_hat
        quietly gen double `Y_hat' = `yhatinput' if `touse'
        quietly gen double `W_hat' = `whatinput' if `touse'
        quietly gen double `Z_hat' = `zhatinput' if `touse'
    }
    else if "`yhatinput'" != "" | "`whatinput'" != "" | "`zhatinput'" != "" {
        display as error "all three of yhatinput(), whatinput(), zhatinput() must be specified together"
        exit 198
    }
    else {

    display as text "Nuisance estimation (3 regression forests)..."

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

    /* Build extra_vars for nuisance calls (must be before nuisance plugin calls) */
    local extra_vars ""
    if "`cluster_var'" != "" {
        local extra_vars `extra_vars' `cluster_var'
    }
    if "`weight_var'" != "" {
        local extra_vars `extra_vars' `weight_var'
    }

    /* ---- Step 1: Y.hat ---- */
    tempvar Y_hat
    quietly gen double `Y_hat' = .

    display as text "  Fitting Y ~ X..."
    plugin call grf_plugin `indepvars' `depvar' `extra_vars' `Y_hat' ///
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
        "1"                                              ///
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

    /* ---- Step 2: W.hat ---- */
    tempvar W_hat
    quietly gen double `W_hat' = .

    display as text "  Fitting W ~ X..."
    plugin call grf_plugin `indepvars' `treatment' `extra_vars' `W_hat' ///
        if `touse',                                         ///
        "regression"                                        ///
        "`nuisancetrees'"                                   ///
        "`seed'"                                            ///
        "`mtry'"                                            ///
        "`minnodesize'"                                     ///
        "`samplefrac'"                                      ///
        "`do_honesty'"                                      ///
        "`honestyfrac'"                                     ///
        "`do_honesty_prune'"                                ///
        "`alpha'"                                           ///
        "`imbalancepenalty'"                                 ///
        "1"                                                 ///
        "`numthreads'"                                      ///
        "0"                                                 ///
        "0"                                                 ///
        "`nindep'"                                          ///
        "1"                                                 ///
        "0"                                                 ///
        "0"                                                 ///
        "1"                                                 ///
        "`allow_missing_x'"                                 ///
        "`_nuis_cluster_idx'"                               ///
        "`_nuis_weight_idx'"

    /* ---- Step 3: Z.hat ---- */
    tempvar Z_hat
    quietly gen double `Z_hat' = .

    display as text "  Fitting Z ~ X..."
    plugin call grf_plugin `indepvars' `instrument' `extra_vars' `Z_hat' ///
        if `touse',                                          ///
        "regression"                                         ///
        "`nuisancetrees'"                                    ///
        "`seed'"                                             ///
        "`mtry'"                                             ///
        "`minnodesize'"                                      ///
        "`samplefrac'"                                       ///
        "`do_honesty'"                                       ///
        "`honestyfrac'"                                      ///
        "`do_honesty_prune'"                                 ///
        "`alpha'"                                            ///
        "`imbalancepenalty'"                                  ///
        "1"                                                  ///
        "`numthreads'"                                       ///
        "0"                                                  ///
        "0"                                                  ///
        "`nindep'"                                           ///
        "1"                                                  ///
        "0"                                                  ///
        "0"                                                  ///
        "1"                                                  ///
        "`allow_missing_x'"                                  ///
        "`_nuis_cluster_idx'"                                ///
        "`_nuis_weight_idx'"

    } /* end else: nuisance forest pipeline */

    /* ---- Step 4: Center Y, W, Z ---- */
    tempvar Y_centered W_centered Z_centered
    quietly gen double `Y_centered' = `depvar' - `Y_hat' if `touse'
    quietly gen double `W_centered' = `treatment' - `W_hat' if `touse'
    quietly gen double `Z_centered' = `instrument' - `Z_hat' if `touse'

    display as text "  Centering complete."

    /* ---- Save nuisance estimates if requested ---- */
    if "`yhatgenerate'" != "" {
        if "`replace'" != "" {
            capture drop `yhatgenerate'
        }
        confirm new variable `yhatgenerate'
        quietly gen double `yhatgenerate' = `Y_hat'
        label variable `yhatgenerate' "Y.hat = E[Y|X]"
    }
    if "`whatgenerate'" != "" {
        if "`replace'" != "" {
            capture drop `whatgenerate'
        }
        confirm new variable `whatgenerate'
        quietly gen double `whatgenerate' = `W_hat'
        label variable `whatgenerate' "W.hat = E[W|X]"
    }
    if "`zhatgenerate'" != "" {
        if "`replace'" != "" {
            capture drop `zhatgenerate'
        }
        confirm new variable `zhatgenerate'
        quietly gen double `zhatgenerate' = `Z_hat'
        label variable `zhatgenerate' "Z.hat = E[Z|X]"
    }

    display as text ""

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
    /* Instrumental forest data vars: X1..Xp Y.c W.c Z.c => nindep + 3 columns */
    local _data_col_count = `nindep' + 3
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

    /* ---- Call plugin for instrumental forest ----
     *
     * Variable order: X1..Xp Y.c W.c Z.c [cluster] [weight] out1 [out2]
     * argv: forest_type num_trees seed mtry min_node_size sample_fraction
     *       honesty honesty_fraction honesty_prune alpha imbalance_penalty
     *       ci_group_size num_threads estimate_variance compute_oob
     *       n_x n_y n_w n_z n_output
     *       allow_missing_x cluster_col_idx weight_col_idx
     *       reduced_form_weight stabilize_splits
     */
    plugin call grf_plugin `indepvars' `Y_centered' `W_centered' ///
        `Z_centered' `extra_vars' `output_vars'                   ///
        if `touse',                                               ///
        "instrumental"                                            ///
        "`ntrees'"                                                ///
        "`seed'"                                                  ///
        "`mtry'"                                                  ///
        "`minnodesize'"                                           ///
        "`samplefrac'"                                            ///
        "`do_honesty'"                                            ///
        "`honestyfrac'"                                           ///
        "`do_honesty_prune'"                                      ///
        "`alpha'"                                                 ///
        "`imbalancepenalty'"                                       ///
        "`cigroupsize'"                                           ///
        "`numthreads'"                                            ///
        "`do_est_var'"                                            ///
        "0"                                                       ///
        "`nindep'"                                                ///
        "1"                                                       ///
        "1"                                                       ///
        "1"                                                       ///
        "`n_output'"                                              ///
        "`allow_missing_x'"                                       ///
        "`cluster_col_idx'"                                       ///
        "`weight_col_idx'"                                        ///
        "`reducedformweight'"                                     ///
        "`do_stabilize'"

    /* ---- Store results ---- */
    ereturn clear
    ereturn scalar N                  = `n_use'
    ereturn scalar n_trees            = `ntrees'
    ereturn scalar seed               = `seed'
    ereturn scalar mtry               = `mtry'
    ereturn scalar min_node           = `minnodesize'
    ereturn scalar alpha              = `alpha'
    ereturn scalar honesty            = `do_honesty'
    ereturn scalar honesty_prune = `do_honesty_prune'
    ereturn scalar sample_fraction    = `samplefrac'
    ereturn scalar honesty_fraction   = `honestyfrac'
    ereturn scalar imbalance_penalty  = `imbalancepenalty'
    ereturn scalar ci_group_size      = `cigroupsize'
    ereturn scalar reduced_form_wt    = `reducedformweight'
    ereturn scalar stabilize_splits   = `do_stabilize'
    ereturn local  cmd                  "grf_instrumental_forest"
    ereturn scalar allow_missing_x = `allow_missing_x'
    ereturn local  forest_type          "instrumental"
    ereturn local  depvar               "`depvar'"
    ereturn local  treatment            "`treatment'"
    ereturn local  instrument           "`instrument'"
    ereturn local  indepvars            "`indepvars'"
    ereturn local  predict_var          "`generate'"
    if `do_est_var' {
        ereturn local variance_var "`vargenerate'"
    }
    if "`cluster_var'" != "" {
        ereturn local cluster_var "`cluster_var'"
    }
    if "`weight_var'" != "" {
        ereturn local weight_var "`weight_var'"
    }
    if "`yhatgenerate'" != "" {
        ereturn local yhat_var "`yhatgenerate'"
    }
    if "`whatgenerate'" != "" {
        ereturn local what_var "`whatgenerate'"
    }
    if "`zhatgenerate'" != "" {
        ereturn local zhat_var "`zhatgenerate'"
    }
    if "`yhatinput'" != "" {
        ereturn local yhat_input "`yhatinput'"
        ereturn local what_input "`whatinput'"
        ereturn local zhat_input "`zhatinput'"
    }

    /* ---- Summary stats ---- */
    quietly summarize `generate' if `touse'
    local n_pred = r(N)
    local pred_mean = r(mean)
    local pred_sd = r(sd)

    display as text ""
    display as text "LATE estimates written to: " as result "`generate'"
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
