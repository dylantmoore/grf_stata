*! grf_tune.ado -- Cross-validation tuning for GRF forests
*! Version 0.2.0
*! Random search over hyperparameters with OOB prediction error evaluation
*! Supports all GRF forest types

program define grf_tune, rclass
    version 14.0

    syntax varlist(min=2 numeric) [if] [in],  ///
        FORESTtype(string)                     ///
        [                                      ///
            NUMReps(integer 50)                ///
            TUNETrees(integer 200)             ///
            SEED(integer 42)                   ///
            NUMThreads(integer 0)              ///
            noHONesty                          ///
            noHONestyprune                     ///
            noSTABilizesplits                  ///
            noMIA                              ///
            CLuster(varname numeric)           ///
            WEIghts(varname numeric)           ///
            Xvars(varlist numeric)             ///
            NTReat(integer 0)                  ///
            NDEp(integer 0)                    ///
            NCLasses(integer 0)                ///
            NUMFailures(integer 0)             ///
            PREDtype(integer 1)                ///
            REDucedformweight(real 0.0)        ///
            HORizon(real 0)                    ///
            TARget(integer 1)                  ///
        ]

    /* ---- Validate forest_type ---- */
    local valid_types "regression causal quantile probability instrumental"
    local valid_types "`valid_types' survival ll_regression lm_forest"
    local valid_types "`valid_types' multi_arm_causal multi_regression"
    local valid_types "`valid_types' causal_survival boosted_regression"

    local type_ok 0
    foreach vt of local valid_types {
        if "`foresttype'" == "`vt'" {
            local type_ok 1
        }
    }
    if !`type_ok' {
        display as error "forest_type() must be one of: `valid_types'"
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

    /* ---- Parse MIA ---- */
    local allow_missing_x 1
    if "`mia'" == "nomia" {
        local allow_missing_x 0
    }

    /* ---- Parse cluster ---- */
    local cluster_var ""
    local cluster_col_idx 0
    if "`cluster'" != "" {
        local cluster_var "`cluster'"
        confirm numeric variable `cluster_var'
    }

    /* ---- Parse weights ---- */
    local weight_var ""
    local weight_col_idx 0
    if "`weights'" != "" {
        local weight_var "`weights'"
        confirm numeric variable `weight_var'
    }

    /* ---- Parse varlist depending on forest type ---- */
    if "`foresttype'" == "causal" {
        /* depvar treatvar indepvars */
        gettoken depvar rest : varlist
        gettoken treatvar indepvars : rest
        local nindep : word count `indepvars'
        if `nindep' < 1 {
            display as error "need at least 1 predictor variable for causal forest"
            exit 198
        }
    }
    else if "`foresttype'" == "instrumental" {
        /* depvar treatvar instrument indepvars */
        gettoken depvar rest : varlist
        gettoken treatvar rest : rest
        gettoken instrument indepvars : rest
        local nindep : word count `indepvars'
        if `nindep' < 1 {
            display as error "need at least 1 predictor variable for instrumental forest"
            exit 198
        }
    }
    else if "`foresttype'" == "survival" {
        /* timevar statusvar indepvars */
        gettoken timevar rest : varlist
        gettoken statusvar indepvars : rest
        local nindep : word count `indepvars'
        if `nindep' < 1 {
            display as error "need at least 1 predictor variable for survival forest"
            exit 198
        }
    }
    else if "`foresttype'" == "lm_forest" {
        /* depvar regressor1 [regressor2 ...] -- xvars required separately */
        if "`xvars'" == "" {
            display as error "xvars() required for lm_forest tuning"
            exit 198
        }
        gettoken depvar regvars : varlist
        local n_regressors : word count `regvars'
        if `n_regressors' < 1 {
            display as error "need at least 1 regressor variable for lm_forest"
            exit 198
        }
        local indepvars `xvars'
        local nindep : word count `indepvars'
        if `nindep' < 1 {
            display as error "xvars() must contain at least 1 covariate"
            exit 198
        }
    }
    else if "`foresttype'" == "multi_arm_causal" {
        /* depvar treat1 treat2 ... indepvars -- ntreat required */
        if `ntreat' < 1 {
            display as error "ntreat() must be at least 1 for multi_arm_causal"
            exit 198
        }
        gettoken depvar rest : varlist
        local treatvars ""
        forvalues j = 1/`ntreat' {
            gettoken tv rest : rest
            local treatvars `treatvars' `tv'
        }
        local indepvars `rest'
        local nindep : word count `indepvars'
        if `nindep' < 1 {
            display as error "need at least 1 predictor variable for multi_arm_causal"
            exit 198
        }
    }
    else if "`foresttype'" == "multi_regression" {
        /* dep1 dep2 ... indepvars -- ndep required */
        if `ndep' < 2 {
            display as error "ndep() must be at least 2 for multi_regression"
            exit 198
        }
        local depvars ""
        local rest `varlist'
        forvalues j = 1/`ndep' {
            gettoken dv rest : rest
            local depvars `depvars' `dv'
        }
        local indepvars `rest'
        local nindep : word count `indepvars'
        if `nindep' < 1 {
            display as error "need at least 1 predictor variable for multi_regression"
            exit 198
        }
    }
    else if "`foresttype'" == "causal_survival" {
        /* timevar statusvar treatvar indepvars */
        gettoken timevar rest : varlist
        gettoken statusvar rest : rest
        gettoken treatvar indepvars : rest
        local nindep : word count `indepvars'
        if `nindep' < 1 {
            display as error "need at least 1 predictor variable for causal_survival"
            exit 198
        }
    }
    else {
        /* regression, quantile, probability, ll_regression, boosted_regression */
        /* depvar indepvars */
        gettoken depvar indepvars : varlist
        local nindep : word count `indepvars'
        if `nindep' < 1 {
            display as error "need at least 1 predictor variable"
            exit 198
        }
    }

    /* ---- Build extra vars for cluster/weight ---- */
    local extra_vars ""
    if "`cluster_var'" != "" {
        local extra_vars `extra_vars' `cluster_var'
    }
    if "`weight_var'" != "" {
        local extra_vars `extra_vars' `weight_var'
    }

    /* ---- Mark sample ---- */
    marksample touse
    if "`foresttype'" == "lm_forest" {
        markout `touse' `xvars'
    }
    quietly count if `touse'
    local n_use = r(N)

    if `n_use' < 10 {
        display as error "need at least 10 non-missing observations for tuning"
        exit 2000
    }

    /* ---- Validate parameters ---- */
    if `numreps' < 1 {
        display as error "num_reps must be at least 1"
        exit 198
    }
    if `tunetrees' < 10 {
        display as error "num_tune_trees must be at least 10"
        exit 198
    }

    /* ---- Type-specific validation ---- */
    if "`foresttype'" == "probability" {
        /* Auto-detect nclasses if not specified */
        if `nclasses' == 0 {
            quietly summarize `depvar' if `touse'
            local nclasses = r(max) + 1
        }
        if `nclasses' < 2 {
            display as error "need at least 2 classes for probability forest"
            exit 198
        }
    }

    if "`foresttype'" == "survival" | "`foresttype'" == "causal_survival" {
        quietly summarize `timevar' if `touse'
        if r(min) <= 0 {
            display as error "survival time must be strictly positive"
            exit 198
        }
    }

    if "`foresttype'" == "causal_survival" & `horizon' == 0 {
        quietly summarize `timevar' if `statusvar' == 1 & `touse', detail
        local horizon = r(p50)
        display as text "(horizon not specified; using median failure time: " ///
            as result %9.3f `horizon' as text ")"
    }

    /* ---- Display header ---- */
    display as text ""
    display as text "GRF Hyperparameter Tuning (Random Search)"
    display as text "{hline 55}"
    display as text "Forest type:           " as result "`foresttype'"
    if "`foresttype'" == "causal" {
        display as text "Outcome variable:      " as result "`depvar'"
        display as text "Treatment variable:    " as result "`treatvar'"
    }
    else if "`foresttype'" == "instrumental" {
        display as text "Outcome variable:      " as result "`depvar'"
        display as text "Treatment variable:    " as result "`treatvar'"
        display as text "Instrument variable:   " as result "`instrument'"
    }
    else if "`foresttype'" == "survival" {
        display as text "Time variable:         " as result "`timevar'"
        display as text "Status variable:       " as result "`statusvar'"
    }
    else if "`foresttype'" == "lm_forest" {
        display as text "Outcome variable:      " as result "`depvar'"
        display as text "Regressor variables:   " as result "`regvars'"
        display as text "Splitting covariates:  " as result "`indepvars'"
    }
    else if "`foresttype'" == "multi_arm_causal" {
        display as text "Outcome variable:      " as result "`depvar'"
        display as text "Treatment arms:        " as result "`treatvars'"
    }
    else if "`foresttype'" == "multi_regression" {
        display as text "Outcome variables:     " as result "`depvars'"
    }
    else if "`foresttype'" == "causal_survival" {
        display as text "Time variable:         " as result "`timevar'"
        display as text "Status variable:       " as result "`statusvar'"
        display as text "Treatment variable:    " as result "`treatvar'"
        display as text "Horizon:               " as result %9.3f `horizon'
    }
    else {
        display as text "Dependent variable:    " as result "`depvar'"
    }
    display as text "Predictors:            " as result "`indepvars'"
    display as text "Observations:          " as result `n_use'
    display as text "Random search reps:    " as result `numreps'
    display as text "Trees per candidate:   " as result `tunetrees'
    display as text "Seed:                  " as result `seed'
    display as text "{hline 55}"
    display as text ""

    /* ---- Load plugin ---- */
    if ( inlist("`c(os)'", "MacOSX") | strpos("`c(machine_type)'", "Mac") ) local c_os_ macosx
    else local c_os_: di lower("`c(os)'")

    capture program grf_plugin, plugin using("grf_plugin_`c_os_'.plugin")

    /* ---- Set seed for reproducibility ---- */
    set seed `seed'

    /* ---- Compute bounds for min_node_size (log-uniform) ---- */
    local log_mn_lo = ln(5)
    local log_mn_hi = ln(max(6, `n_use' / 4))

    /* ---- Initialize best tracking ---- */
    local best_mse       = .
    local best_mtry       = 0
    local best_minnodesize = 5
    local best_samplefrac = 0.5
    local best_honestyfrac = 0.5
    local best_alpha      = 0.05
    local best_imbpen     = 0.0

    /* ==================================================================
     * Nuisance estimation (once, with defaults) for forests that need it
     * Following R's grf approach: estimate nuisance once, then only tune
     * the main forest hyperparameters.
     * ================================================================== */

    if "`foresttype'" == "causal" {
        display as text "Fitting nuisance models (Y ~ X and W ~ X) ..."
        tempvar yhat_nuis what_nuis
        quietly gen double `yhat_nuis' = .
        quietly gen double `what_nuis' = .

        /* Compute nuisance column indices for cluster/weight */
        local _nuis_cluster_idx 0
        local _nuis_weight_idx 0
        if "`cluster_var'" != "" | "`weight_var'" != "" {
            local _nuis_data_col_count = `nindep' + 1
            local _nuis_offset 0
            if "`cluster_var'" != "" {
                local _nuis_cluster_idx = `_nuis_data_col_count' + 1
                local _nuis_offset 1
            }
            if "`weight_var'" != "" {
                local _nuis_weight_idx = `_nuis_data_col_count' + `_nuis_offset' + 1
            }
        }

        /* Y ~ X regression forest */
        plugin call grf_plugin `indepvars' `depvar' `extra_vars' `yhat_nuis' ///
            if `touse',                                          ///
            "regression"                                         ///
            "`tunetrees'"                                        ///
            "`seed'"                                             ///
            "0"                                                  ///
            "5"                                                  ///
            "0.5"                                                ///
            "`do_honesty'"                                       ///
            "0.5"                                                ///
            "`do_honesty_prune'"                                 ///
            "0.05"                                               ///
            "0.0"                                                ///
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

        /* W ~ X regression forest */
        plugin call grf_plugin `indepvars' `treatvar' `extra_vars' `what_nuis' ///
            if `touse',                                            ///
            "regression"                                           ///
            "`tunetrees'"                                          ///
            "`seed'"                                               ///
            "0"                                                    ///
            "5"                                                    ///
            "0.5"                                                  ///
            "`do_honesty'"                                         ///
            "0.5"                                                  ///
            "`do_honesty_prune'"                                   ///
            "0.05"                                                 ///
            "0.0"                                                  ///
            "1"                                                    ///
            "`numthreads'"                                         ///
            "0"                                                    ///
            "0"                                                    ///
            "`nindep'"                                             ///
            "1"                                                    ///
            "0"                                                    ///
            "0"                                                    ///
            "1"                                                    ///
            "`allow_missing_x'"                                    ///
            "`_nuis_cluster_idx'"                                  ///
            "`_nuis_weight_idx'"

        /* Center Y and W */
        tempvar y_centered w_centered
        quietly gen double `y_centered' = `depvar' - `yhat_nuis' if `touse'
        quietly gen double `w_centered' = `treatvar' - `what_nuis' if `touse'

        display as text "Nuisance models fitted. Beginning random search ..."
        display as text ""
    }
    else if "`foresttype'" == "instrumental" {
        display as text "Fitting nuisance models (Y ~ X, W ~ X, Z ~ X) ..."
        tempvar yhat_nuis what_nuis zhat_nuis
        quietly gen double `yhat_nuis' = .
        quietly gen double `what_nuis' = .
        quietly gen double `zhat_nuis' = .

        /* Compute nuisance column indices for cluster/weight */
        local _nuis_cluster_idx 0
        local _nuis_weight_idx 0
        if "`cluster_var'" != "" | "`weight_var'" != "" {
            local _nuis_data_col_count = `nindep' + 1
            local _nuis_offset 0
            if "`cluster_var'" != "" {
                local _nuis_cluster_idx = `_nuis_data_col_count' + 1
                local _nuis_offset 1
            }
            if "`weight_var'" != "" {
                local _nuis_weight_idx = `_nuis_data_col_count' + `_nuis_offset' + 1
            }
        }

        /* Y ~ X */
        plugin call grf_plugin `indepvars' `depvar' `extra_vars' `yhat_nuis' ///
            if `touse',                                          ///
            "regression"                                         ///
            "`tunetrees'"                                        ///
            "`seed'"                                             ///
            "0"                                                  ///
            "5"                                                  ///
            "0.5"                                                ///
            "`do_honesty'"                                       ///
            "0.5"                                                ///
            "`do_honesty_prune'"                                 ///
            "0.05"                                               ///
            "0.0"                                                ///
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

        /* W ~ X */
        plugin call grf_plugin `indepvars' `treatvar' `extra_vars' `what_nuis' ///
            if `touse',                                            ///
            "regression"                                           ///
            "`tunetrees'"                                          ///
            "`seed'"                                               ///
            "0"                                                    ///
            "5"                                                    ///
            "0.5"                                                  ///
            "`do_honesty'"                                         ///
            "0.5"                                                  ///
            "`do_honesty_prune'"                                   ///
            "0.05"                                                 ///
            "0.0"                                                  ///
            "1"                                                    ///
            "`numthreads'"                                         ///
            "0"                                                    ///
            "0"                                                    ///
            "`nindep'"                                             ///
            "1"                                                    ///
            "0"                                                    ///
            "0"                                                    ///
            "1"                                                    ///
            "`allow_missing_x'"                                    ///
            "`_nuis_cluster_idx'"                                  ///
            "`_nuis_weight_idx'"

        /* Z ~ X */
        plugin call grf_plugin `indepvars' `instrument' `extra_vars' `zhat_nuis' ///
            if `touse',                                              ///
            "regression"                                             ///
            "`tunetrees'"                                            ///
            "`seed'"                                                 ///
            "0"                                                      ///
            "5"                                                      ///
            "0.5"                                                    ///
            "`do_honesty'"                                           ///
            "0.5"                                                    ///
            "`do_honesty_prune'"                                     ///
            "0.05"                                                   ///
            "0.0"                                                    ///
            "1"                                                      ///
            "`numthreads'"                                           ///
            "0"                                                      ///
            "0"                                                      ///
            "`nindep'"                                               ///
            "1"                                                      ///
            "0"                                                      ///
            "0"                                                      ///
            "1"                                                      ///
            "`allow_missing_x'"                                      ///
            "`_nuis_cluster_idx'"                                    ///
            "`_nuis_weight_idx'"

        tempvar y_centered w_centered z_centered
        quietly gen double `y_centered' = `depvar' - `yhat_nuis' if `touse'
        quietly gen double `w_centered' = `treatvar' - `what_nuis' if `touse'
        quietly gen double `z_centered' = `instrument' - `zhat_nuis' if `touse'

        display as text "Nuisance models fitted. Beginning random search ..."
        display as text ""
    }
    else if "`foresttype'" == "lm_forest" {
        display as text "Fitting nuisance models (Y ~ X, W_k ~ X) ..."

        /* Compute nuisance column indices for cluster/weight */
        local _nuis_cluster_idx 0
        local _nuis_weight_idx 0
        if "`cluster_var'" != "" | "`weight_var'" != "" {
            local _nuis_data_col_count = `nindep' + 1
            local _nuis_offset 0
            if "`cluster_var'" != "" {
                local _nuis_cluster_idx = `_nuis_data_col_count' + 1
                local _nuis_offset 1
            }
            if "`weight_var'" != "" {
                local _nuis_weight_idx = `_nuis_data_col_count' + `_nuis_offset' + 1
            }
        }

        /* Y ~ X */
        tempvar yhat_nuis
        quietly gen double `yhat_nuis' = .
        plugin call grf_plugin `indepvars' `depvar' `extra_vars' `yhat_nuis' ///
            if `touse',                                          ///
            "regression"                                         ///
            "`tunetrees'"                                        ///
            "`seed'"                                             ///
            "0"                                                  ///
            "5"                                                  ///
            "0.5"                                                ///
            "`do_honesty'"                                       ///
            "0.5"                                                ///
            "`do_honesty_prune'"                                 ///
            "0.05"                                               ///
            "0.0"                                                ///
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

        /* W_k ~ X for each regressor */
        local w_centered_vars ""
        forvalues j = 1/`n_regressors' {
            local wv : word `j' of `regvars'
            tempvar what_`j'
            quietly gen double `what_`j'' = .
            plugin call grf_plugin `indepvars' `wv' `extra_vars' `what_`j'' ///
                if `touse',                                     ///
                "regression"                                    ///
                "`tunetrees'"                                   ///
                "`seed'"                                        ///
                "0"                                             ///
                "5"                                             ///
                "0.5"                                           ///
                "`do_honesty'"                                  ///
                "0.5"                                           ///
                "`do_honesty_prune'"                            ///
                "0.05"                                          ///
                "0.0"                                           ///
                "1"                                             ///
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

        tempvar y_centered
        quietly gen double `y_centered' = `depvar' - `yhat_nuis' if `touse'

        display as text "Nuisance models fitted. Beginning random search ..."
        display as text ""
    }
    else if "`foresttype'" == "multi_arm_causal" {
        display as text "Fitting nuisance models (Y ~ X, W_k ~ X) ..."

        /* Compute nuisance column indices for cluster/weight */
        local _nuis_cluster_idx 0
        local _nuis_weight_idx 0
        if "`cluster_var'" != "" | "`weight_var'" != "" {
            local _nuis_data_col_count = `nindep' + 1
            local _nuis_offset 0
            if "`cluster_var'" != "" {
                local _nuis_cluster_idx = `_nuis_data_col_count' + 1
                local _nuis_offset 1
            }
            if "`weight_var'" != "" {
                local _nuis_weight_idx = `_nuis_data_col_count' + `_nuis_offset' + 1
            }
        }

        /* Y ~ X */
        tempvar yhat_nuis
        quietly gen double `yhat_nuis' = .
        plugin call grf_plugin `indepvars' `depvar' `extra_vars' `yhat_nuis' ///
            if `touse',                                          ///
            "regression"                                         ///
            "`tunetrees'"                                        ///
            "`seed'"                                             ///
            "0"                                                  ///
            "5"                                                  ///
            "0.5"                                                ///
            "`do_honesty'"                                       ///
            "0.5"                                                ///
            "`do_honesty_prune'"                                 ///
            "0.05"                                               ///
            "0.0"                                                ///
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

        /* W_k ~ X for each treatment arm */
        local w_centered_vars ""
        forvalues j = 1/`ntreat' {
            local tv : word `j' of `treatvars'
            tempvar what_`j'
            quietly gen double `what_`j'' = .
            plugin call grf_plugin `indepvars' `tv' `extra_vars' `what_`j'' ///
                if `touse',                                     ///
                "regression"                                    ///
                "`tunetrees'"                                   ///
                "`seed'"                                        ///
                "0"                                             ///
                "5"                                             ///
                "0.5"                                           ///
                "`do_honesty'"                                  ///
                "0.5"                                           ///
                "`do_honesty_prune'"                            ///
                "0.05"                                          ///
                "0.0"                                           ///
                "1"                                             ///
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
            quietly gen double `wc_`j'' = `tv' - `what_`j'' if `touse'
            local w_centered_vars `w_centered_vars' `wc_`j''
        }

        tempvar y_centered
        quietly gen double `y_centered' = `depvar' - `yhat_nuis' if `touse'

        display as text "Nuisance models fitted. Beginning random search ..."
        display as text ""
    }
    else if "`foresttype'" == "causal_survival" {
        display as text "Fitting nuisance model (W ~ X) ..."

        /* Compute nuisance column indices for cluster/weight */
        local _nuis_cluster_idx 0
        local _nuis_weight_idx 0
        if "`cluster_var'" != "" | "`weight_var'" != "" {
            local _nuis_data_col_count = `nindep' + 1
            local _nuis_offset 0
            if "`cluster_var'" != "" {
                local _nuis_cluster_idx = `_nuis_data_col_count' + 1
                local _nuis_offset 1
            }
            if "`weight_var'" != "" {
                local _nuis_weight_idx = `_nuis_data_col_count' + `_nuis_offset' + 1
            }
        }

        tempvar what_nuis
        quietly gen double `what_nuis' = .
        plugin call grf_plugin `indepvars' `treatvar' `extra_vars' `what_nuis' ///
            if `touse',                                            ///
            "regression"                                           ///
            "`tunetrees'"                                          ///
            "`seed'"                                               ///
            "0"                                                    ///
            "5"                                                    ///
            "0.5"                                                  ///
            "`do_honesty'"                                         ///
            "0.5"                                                  ///
            "`do_honesty_prune'"                                   ///
            "0.05"                                                 ///
            "0.0"                                                  ///
            "1"                                                    ///
            "`numthreads'"                                         ///
            "0"                                                    ///
            "0"                                                    ///
            "`nindep'"                                             ///
            "1"                                                    ///
            "0"                                                    ///
            "0"                                                    ///
            "1"                                                    ///
            "`allow_missing_x'"                                    ///
            "`_nuis_cluster_idx'"                                  ///
            "`_nuis_weight_idx'"

        /* Simplified nuisance: IPCW approximation */
        tempvar w_centered cs_numer cs_denom
        quietly gen double `w_centered' = `treatvar' - `what_nuis' if `touse'
        quietly gen double `cs_numer' = `w_centered' * `statusvar' * ///
            min(`timevar', `horizon') if `touse'
        quietly gen double `cs_denom' = 1 if `touse'

        display as text "Nuisance model fitted. Beginning random search ..."
        display as text ""
    }
    else {
        display as text "Beginning random search ..."
        display as text ""
    }

    /* ---- Compute main tuning column indices for cluster/weight ---- */
    local _main_cluster_idx 0
    local _main_weight_idx 0
    if "`cluster_var'" != "" | "`weight_var'" != "" {
        if "`foresttype'" == "regression" {
            local _data_col_count = `nindep' + 1
        }
        else if "`foresttype'" == "causal" {
            local _data_col_count = `nindep' + 2
        }
        else if "`foresttype'" == "quantile" {
            local _data_col_count = `nindep' + 1
        }
        else if "`foresttype'" == "probability" {
            local _data_col_count = `nindep' + 1
        }
        else if "`foresttype'" == "instrumental" {
            local _data_col_count = `nindep' + 3
        }
        else if "`foresttype'" == "survival" {
            local _data_col_count = `nindep' + 2
        }
        else if "`foresttype'" == "ll_regression" {
            local _data_col_count = `nindep' + 1
        }
        else if "`foresttype'" == "lm_forest" {
            local _data_col_count = `nindep' + 1 + `n_regressors'
        }
        else if "`foresttype'" == "multi_arm_causal" {
            local _data_col_count = `nindep' + 1 + `ntreat'
        }
        else if "`foresttype'" == "multi_regression" {
            local _data_col_count = `nindep' + `ndep'
        }
        else if "`foresttype'" == "causal_survival" {
            local _data_col_count = `nindep' + 5
        }
        else if "`foresttype'" == "boosted_regression" {
            local _data_col_count = `nindep' + 1
        }

        local _main_offset 0
        if "`cluster_var'" != "" {
            local _main_cluster_idx = `_data_col_count' + 1
            local _main_offset 1
        }
        if "`weight_var'" != "" {
            local _main_weight_idx = `_data_col_count' + `_main_offset' + 1
        }
    }

    /* ---- Random search loop ---- */
    forvalues rep = 1/`numreps' {

        /* Generate random candidate parameters */
        local cand_samplefrac  = 0.05 + runiform() * (0.50 - 0.05)
        local cand_mtry        = max(1, ceil(runiform() * `nindep'))
        local cand_minnodesize = max(1, round(exp(`log_mn_lo' + runiform() * (`log_mn_hi' - `log_mn_lo'))))
        local cand_honestyfrac = 0.50 + runiform() * (0.80 - 0.50)
        local cand_alpha       = runiform() * 0.25
        local cand_imbpen      = runiform() * 3.0

        /* Create temp output variable */
        tempvar tune_pred
        quietly gen double `tune_pred' = .

        /* ---- Fit candidate forest and get OOB predictions ---- */

        if "`foresttype'" == "regression" {
            capture plugin call grf_plugin `indepvars' `depvar' `extra_vars' `tune_pred' ///
                if `touse',                                                  ///
                "regression"                                                 ///
                "`tunetrees'"                                                ///
                "`seed'"                                                     ///
                "`cand_mtry'"                                                ///
                "`cand_minnodesize'"                                         ///
                "`cand_samplefrac'"                                          ///
                "`do_honesty'"                                               ///
                "`cand_honestyfrac'"                                         ///
                "`do_honesty_prune'"                                         ///
                "`cand_alpha'"                                               ///
                "`cand_imbpen'"                                              ///
                "1"                                                          ///
                "`numthreads'"                                               ///
                "0"                                                          ///
                "0"                                                          ///
                "`nindep'"                                                   ///
                "1"                                                          ///
                "0"                                                          ///
                "0"                                                          ///
                "1"                                                          ///
                "`allow_missing_x'"                                          ///
                "`_main_cluster_idx'"                                        ///
                "`_main_weight_idx'"

            if _rc {
                capture drop `tune_pred'
                continue
            }

            /* OOB MSE: mean((Y - Yhat)^2) */
            tempvar sq_err
            quietly gen double `sq_err' = (`depvar' - `tune_pred')^2 if `touse' & !missing(`tune_pred')
            quietly summarize `sq_err' if `touse' & !missing(`tune_pred')
            local cand_mse = r(mean)
            capture drop `sq_err'
        }
        else if "`foresttype'" == "causal" {
            capture plugin call grf_plugin `indepvars' `y_centered' `w_centered' `extra_vars' `tune_pred' ///
                if `touse',                                                                   ///
                "causal"                                                                      ///
                "`tunetrees'"                                                                 ///
                "`seed'"                                                                      ///
                "`cand_mtry'"                                                                 ///
                "`cand_minnodesize'"                                                          ///
                "`cand_samplefrac'"                                                           ///
                "`do_honesty'"                                                                ///
                "`cand_honestyfrac'"                                                          ///
                "`do_honesty_prune'"                                                          ///
                "`cand_alpha'"                                                                ///
                "`cand_imbpen'"                                                               ///
                "1"                                                                           ///
                "`numthreads'"                                                                ///
                "0"                                                                           ///
                "0"                                                                           ///
                "`nindep'"                                                                    ///
                "1"                                                                           ///
                "1"                                                                           ///
                "0"                                                                           ///
                "1"                                                                           ///
                "`allow_missing_x'"                                                           ///
                "`_main_cluster_idx'"                                                         ///
                "`_main_weight_idx'"                                                          ///
                "`do_stabilize'"

            if _rc {
                capture drop `tune_pred'
                continue
            }

            /* OOB MSE for causal: mean((Y.centered - tau_hat * W.centered)^2) */
            tempvar sq_err
            quietly gen double `sq_err' = (`y_centered' - `tune_pred' * `w_centered')^2 ///
                if `touse' & !missing(`tune_pred')
            quietly summarize `sq_err' if `touse' & !missing(`tune_pred')
            local cand_mse = r(mean)
            capture drop `sq_err'
        }
        else if "`foresttype'" == "quantile" {
            /* Tune on median prediction MSE */
            capture plugin call grf_plugin `indepvars' `depvar' `extra_vars' `tune_pred' ///
                if `touse',                                                  ///
                "quantile"                                                   ///
                "`tunetrees'"                                                ///
                "`seed'"                                                     ///
                "`cand_mtry'"                                                ///
                "`cand_minnodesize'"                                         ///
                "`cand_samplefrac'"                                          ///
                "`do_honesty'"                                               ///
                "`cand_honestyfrac'"                                         ///
                "`do_honesty_prune'"                                         ///
                "`cand_alpha'"                                               ///
                "`cand_imbpen'"                                              ///
                "1"                                                          ///
                "`numthreads'"                                               ///
                "0"                                                          ///
                "0"                                                          ///
                "`nindep'"                                                   ///
                "1"                                                          ///
                "0"                                                          ///
                "0"                                                          ///
                "1"                                                          ///
                "`allow_missing_x'"                                          ///
                "`_main_cluster_idx'"                                        ///
                "`_main_weight_idx'"                                         ///
                "0.5"

            if _rc {
                capture drop `tune_pred'
                continue
            }

            /* OOB MSE: mean((Y - median_hat)^2) */
            tempvar sq_err
            quietly gen double `sq_err' = (`depvar' - `tune_pred')^2 if `touse' & !missing(`tune_pred')
            quietly summarize `sq_err' if `touse' & !missing(`tune_pred')
            local cand_mse = r(mean)
            capture drop `sq_err'
        }
        else if "`foresttype'" == "probability" {
            /* Tune on predicted probability for the true class (log-loss proxy) */
            /* Use class 0 probability only and compute MSE against indicator */
            tempvar tune_pred_c0 tune_pred_c1
            quietly gen double `tune_pred_c0' = .
            quietly gen double `tune_pred_c1' = .

            /* For tuning, use 2-class output to get P(Y=0) and P(Y=1) */
            capture plugin call grf_plugin `indepvars' `depvar' `extra_vars' `tune_pred_c0' `tune_pred_c1' ///
                if `touse',                                                                    ///
                "probability"                                                                  ///
                "`tunetrees'"                                                                  ///
                "`seed'"                                                                       ///
                "`cand_mtry'"                                                                  ///
                "`cand_minnodesize'"                                                           ///
                "`cand_samplefrac'"                                                            ///
                "`do_honesty'"                                                                 ///
                "`cand_honestyfrac'"                                                           ///
                "`do_honesty_prune'"                                                           ///
                "`cand_alpha'"                                                                 ///
                "`cand_imbpen'"                                                                ///
                "1"                                                                            ///
                "`numthreads'"                                                                 ///
                "0"                                                                            ///
                "0"                                                                            ///
                "`nindep'"                                                                     ///
                "1"                                                                            ///
                "0"                                                                            ///
                "0"                                                                            ///
                "`nclasses'"                                                                   ///
                "`allow_missing_x'"                                                            ///
                "`_main_cluster_idx'"                                                          ///
                "`_main_weight_idx'"                                                           ///
                "`nclasses'"

            if _rc {
                capture drop `tune_pred' `tune_pred_c0' `tune_pred_c1'
                continue
            }

            /* Brier score: mean((1{Y=1} - P(Y=1))^2) for binary; generalize later */
            tempvar sq_err
            quietly gen double `sq_err' = ((`depvar' == 1) - `tune_pred_c1')^2 ///
                if `touse' & !missing(`tune_pred_c1')
            quietly summarize `sq_err' if `touse' & !missing(`tune_pred_c1')
            local cand_mse = r(mean)
            capture drop `sq_err' `tune_pred_c0' `tune_pred_c1'
        }
        else if "`foresttype'" == "instrumental" {
            /* Instrumental forest: tune on pseudo-outcome residual MSE */
            capture plugin call grf_plugin `indepvars' `y_centered' `w_centered' ///
                `z_centered' `extra_vars' `tune_pred'                             ///
                if `touse',                                                       ///
                "instrumental"                                                    ///
                "`tunetrees'"                                                     ///
                "`seed'"                                                          ///
                "`cand_mtry'"                                                     ///
                "`cand_minnodesize'"                                              ///
                "`cand_samplefrac'"                                               ///
                "`do_honesty'"                                                    ///
                "`cand_honestyfrac'"                                              ///
                "`do_honesty_prune'"                                              ///
                "`cand_alpha'"                                                    ///
                "`cand_imbpen'"                                                   ///
                "1"                                                               ///
                "`numthreads'"                                                    ///
                "0"                                                               ///
                "0"                                                               ///
                "`nindep'"                                                        ///
                "1"                                                               ///
                "1"                                                               ///
                "1"                                                               ///
                "1"                                                               ///
                "`allow_missing_x'"                                               ///
                "`_main_cluster_idx'"                                             ///
                "`_main_weight_idx'"                                              ///
                "`reducedformweight'"                                             ///
                "`do_stabilize'"

            if _rc {
                capture drop `tune_pred'
                continue
            }

            /* MSE: mean((Y.c - tau_hat * W.c)^2) */
            tempvar sq_err
            quietly gen double `sq_err' = (`y_centered' - `tune_pred' * `w_centered')^2 ///
                if `touse' & !missing(`tune_pred')
            quietly summarize `sq_err' if `touse' & !missing(`tune_pred')
            local cand_mse = r(mean)
            capture drop `sq_err'
        }
        else if "`foresttype'" == "survival" {
            /* Survival forest: tune on first survival curve column MSE vs event indicator */
            capture plugin call grf_plugin `indepvars' `timevar' `statusvar' `extra_vars' `tune_pred' ///
                if `touse',                                                               ///
                "survival"                                                                ///
                "`tunetrees'"                                                             ///
                "`seed'"                                                                  ///
                "`cand_mtry'"                                                             ///
                "`cand_minnodesize'"                                                      ///
                "`cand_samplefrac'"                                                       ///
                "`do_honesty'"                                                            ///
                "`cand_honestyfrac'"                                                      ///
                "`do_honesty_prune'"                                                      ///
                "`cand_alpha'"                                                            ///
                "`cand_imbpen'"                                                           ///
                "1"                                                                       ///
                "`numthreads'"                                                            ///
                "0"                                                                       ///
                "0"                                                                       ///
                "`nindep'"                                                                ///
                "1"                                                                       ///
                "0"                                                                       ///
                "0"                                                                       ///
                "1"                                                                       ///
                "`allow_missing_x'"                                                       ///
                "`_main_cluster_idx'"                                                     ///
                "`_main_weight_idx'"                                                      ///
                "`numfailures'"                                                           ///
                "`=1 - `predtype''"

            if _rc {
                capture drop `tune_pred'
                continue
            }

            /* Brier-like score: mean((1-D_i - S_hat)^2) where S_hat is the survival prob */
            tempvar sq_err
            quietly gen double `sq_err' = ((1 - `statusvar') - `tune_pred')^2 ///
                if `touse' & !missing(`tune_pred')
            quietly summarize `sq_err' if `touse' & !missing(`tune_pred')
            local cand_mse = r(mean)
            capture drop `sq_err'
        }
        else if "`foresttype'" == "ll_regression" {
            /* Local linear regression forest: same as regression but different plugin call */
            capture plugin call grf_plugin `indepvars' `depvar' `extra_vars' `tune_pred' ///
                if `touse',                                                  ///
                "ll_regression"                                              ///
                "`tunetrees'"                                                ///
                "`seed'"                                                     ///
                "`cand_mtry'"                                                ///
                "`cand_minnodesize'"                                         ///
                "`cand_samplefrac'"                                          ///
                "`do_honesty'"                                               ///
                "`cand_honestyfrac'"                                         ///
                "`do_honesty_prune'"                                         ///
                "`cand_alpha'"                                               ///
                "`cand_imbpen'"                                              ///
                "1"                                                          ///
                "`numthreads'"                                               ///
                "0"                                                          ///
                "0"                                                          ///
                "`nindep'"                                                   ///
                "1"                                                          ///
                "0"                                                          ///
                "0"                                                          ///
                "1"                                                          ///
                "`allow_missing_x'"                                          ///
                "`_main_cluster_idx'"                                        ///
                "`_main_weight_idx'"                                         ///
                "0"                                                          ///
                "0.1"                                                        ///
                "0"                                                          ///
                "0"

            if _rc {
                capture drop `tune_pred'
                continue
            }

            /* OOB MSE */
            tempvar sq_err
            quietly gen double `sq_err' = (`depvar' - `tune_pred')^2 if `touse' & !missing(`tune_pred')
            quietly summarize `sq_err' if `touse' & !missing(`tune_pred')
            local cand_mse = r(mean)
            capture drop `sq_err'
        }
        else if "`foresttype'" == "lm_forest" {
            /* Linear model forest: tune on centered outcome residual */
            /* Output: one coefficient per regressor */
            local tune_output_vars ""
            forvalues j = 1/`n_regressors' {
                tempvar tune_coef_`j'
                quietly gen double `tune_coef_`j'' = .
                local tune_output_vars `tune_output_vars' `tune_coef_`j''
            }

            capture plugin call grf_plugin `indepvars' `y_centered' `w_centered_vars' `extra_vars' `tune_output_vars' ///
                if `touse',                                                                               ///
                "lm_forest"                                                                              ///
                "`tunetrees'"                                                                            ///
                "`seed'"                                                                                 ///
                "`cand_mtry'"                                                                            ///
                "`cand_minnodesize'"                                                                     ///
                "`cand_samplefrac'"                                                                      ///
                "`do_honesty'"                                                                           ///
                "`cand_honestyfrac'"                                                                     ///
                "`do_honesty_prune'"                                                                     ///
                "`cand_alpha'"                                                                           ///
                "`cand_imbpen'"                                                                          ///
                "1"                                                                                      ///
                "`numthreads'"                                                                           ///
                "0"                                                                                      ///
                "0"                                                                                      ///
                "`nindep'"                                                                               ///
                "1"                                                                                      ///
                "`n_regressors'"                                                                         ///
                "0"                                                                                      ///
                "`n_regressors'"                                                                         ///
                "`allow_missing_x'"                                                                      ///
                "`_main_cluster_idx'"                                                                    ///
                "`_main_weight_idx'"                                                                     ///
                "`do_stabilize'"

            if _rc {
                capture drop `tune_pred'
                forvalues j = 1/`n_regressors' {
                    capture drop `tune_coef_`j''
                }
                continue
            }

            /* MSE: mean((Y.c - sum_k(coef_k * W_k.c))^2) */
            tempvar fitted_vals
            quietly gen double `fitted_vals' = 0 if `touse'
            forvalues j = 1/`n_regressors' {
                local wcv : word `j' of `w_centered_vars'
                quietly replace `fitted_vals' = `fitted_vals' + `tune_coef_`j'' * `wcv' ///
                    if `touse' & !missing(`tune_coef_`j'')
            }
            tempvar sq_err
            quietly gen double `sq_err' = (`y_centered' - `fitted_vals')^2 ///
                if `touse' & !missing(`fitted_vals')
            quietly summarize `sq_err' if `touse' & !missing(`fitted_vals')
            local cand_mse = r(mean)
            capture drop `sq_err' `fitted_vals'
            forvalues j = 1/`n_regressors' {
                capture drop `tune_coef_`j''
            }
        }
        else if "`foresttype'" == "multi_arm_causal" {
            /* Multi-arm causal: output one CATE per treatment arm */
            local tune_output_vars ""
            forvalues j = 1/`ntreat' {
                tempvar tune_tau_`j'
                quietly gen double `tune_tau_`j'' = .
                local tune_output_vars `tune_output_vars' `tune_tau_`j''
            }

            capture plugin call grf_plugin `indepvars' `y_centered' `w_centered_vars' `extra_vars' `tune_output_vars' ///
                if `touse',                                                                               ///
                "multi_arm_causal"                                                                       ///
                "`tunetrees'"                                                                            ///
                "`seed'"                                                                                 ///
                "`cand_mtry'"                                                                            ///
                "`cand_minnodesize'"                                                                     ///
                "`cand_samplefrac'"                                                                      ///
                "`do_honesty'"                                                                           ///
                "`cand_honestyfrac'"                                                                     ///
                "`do_honesty_prune'"                                                                     ///
                "`cand_alpha'"                                                                           ///
                "`cand_imbpen'"                                                                          ///
                "1"                                                                                      ///
                "`numthreads'"                                                                           ///
                "0"                                                                                      ///
                "0"                                                                                      ///
                "`nindep'"                                                                               ///
                "1"                                                                                      ///
                "`ntreat'"                                                                               ///
                "0"                                                                                      ///
                "`ntreat'"                                                                               ///
                "`allow_missing_x'"                                                                      ///
                "`_main_cluster_idx'"                                                                    ///
                "`_main_weight_idx'"                                                                     ///
                "`do_stabilize'"                                                                         ///
                "`ntreat'"

            if _rc {
                capture drop `tune_pred'
                forvalues j = 1/`ntreat' {
                    capture drop `tune_tau_`j''
                }
                continue
            }

            /* MSE: mean((Y.c - sum_k(tau_k * W_k.c))^2) */
            tempvar fitted_vals
            quietly gen double `fitted_vals' = 0 if `touse'
            forvalues j = 1/`ntreat' {
                local wcv : word `j' of `w_centered_vars'
                quietly replace `fitted_vals' = `fitted_vals' + `tune_tau_`j'' * `wcv' ///
                    if `touse' & !missing(`tune_tau_`j'')
            }
            tempvar sq_err
            quietly gen double `sq_err' = (`y_centered' - `fitted_vals')^2 ///
                if `touse' & !missing(`fitted_vals')
            quietly summarize `sq_err' if `touse' & !missing(`fitted_vals')
            local cand_mse = r(mean)
            capture drop `sq_err' `fitted_vals'
            forvalues j = 1/`ntreat' {
                capture drop `tune_tau_`j''
            }
        }
        else if "`foresttype'" == "multi_regression" {
            /* Multi-output regression: one prediction per outcome */
            local tune_output_vars ""
            forvalues j = 1/`ndep' {
                tempvar tune_yhat_`j'
                quietly gen double `tune_yhat_`j'' = .
                local tune_output_vars `tune_output_vars' `tune_yhat_`j''
            }

            capture plugin call grf_plugin `indepvars' `depvars' `extra_vars' `tune_output_vars' ///
                if `touse',                                                          ///
                "multi_regression"                                                   ///
                "`tunetrees'"                                                        ///
                "`seed'"                                                             ///
                "`cand_mtry'"                                                        ///
                "`cand_minnodesize'"                                                 ///
                "`cand_samplefrac'"                                                  ///
                "`do_honesty'"                                                       ///
                "`cand_honestyfrac'"                                                 ///
                "`do_honesty_prune'"                                                 ///
                "`cand_alpha'"                                                       ///
                "`cand_imbpen'"                                                      ///
                "1"                                                                  ///
                "`numthreads'"                                                       ///
                "0"                                                                  ///
                "0"                                                                  ///
                "`nindep'"                                                           ///
                "`ndep'"                                                             ///
                "0"                                                                  ///
                "0"                                                                  ///
                "`ndep'"                                                             ///
                "`allow_missing_x'"                                                  ///
                "`_main_cluster_idx'"                                                ///
                "`_main_weight_idx'"                                                 ///
                "`ndep'"

            if _rc {
                capture drop `tune_pred'
                forvalues j = 1/`ndep' {
                    capture drop `tune_yhat_`j''
                }
                continue
            }

            /* Average MSE across all outcomes */
            tempvar total_sq_err
            quietly gen double `total_sq_err' = 0 if `touse'
            local any_missing 0
            forvalues j = 1/`ndep' {
                local dv : word `j' of `depvars'
                quietly replace `total_sq_err' = `total_sq_err' + (`dv' - `tune_yhat_`j'')^2 ///
                    if `touse' & !missing(`tune_yhat_`j'')
            }
            quietly summarize `total_sq_err' if `touse'
            local cand_mse = r(mean) / `ndep'
            capture drop `total_sq_err'
            forvalues j = 1/`ndep' {
                capture drop `tune_yhat_`j''
            }
        }
        else if "`foresttype'" == "causal_survival" {
            /* Causal survival: tune the main forest */
            capture plugin call grf_plugin `indepvars' `timevar' `w_centered' ///
                `statusvar' `cs_numer' `cs_denom' `extra_vars' `tune_pred'     ///
                if `touse',                                                    ///
                "causal_survival"                                              ///
                "`tunetrees'"                                                  ///
                "`seed'"                                                       ///
                "`cand_mtry'"                                                  ///
                "`cand_minnodesize'"                                           ///
                "`cand_samplefrac'"                                            ///
                "`do_honesty'"                                                 ///
                "`cand_honestyfrac'"                                           ///
                "`do_honesty_prune'"                                           ///
                "`cand_alpha'"                                                 ///
                "`cand_imbpen'"                                                ///
                "1"                                                            ///
                "`numthreads'"                                                 ///
                "0"                                                            ///
                "0"                                                            ///
                "`nindep'"                                                     ///
                "1"                                                            ///
                "1"                                                            ///
                "0"                                                            ///
                "1"                                                            ///
                "`allow_missing_x'"                                            ///
                "`_main_cluster_idx'"                                          ///
                "`_main_weight_idx'"                                           ///
                "`do_stabilize'"                                               ///
                "`=`nindep'+3'"                                                ///
                "`=`nindep'+4'"                                                ///
                "`=`nindep'+2'"                                                ///
                "`target'"

            if _rc {
                capture drop `tune_pred'
                continue
            }

            /* MSE of CATE predictions (variance of tau_hat is a proxy for tuning) */
            /* Use pseudo-outcome MSE: mean((numer/denom - tau_hat)^2) */
            tempvar sq_err pseudo_out
            quietly gen double `pseudo_out' = `cs_numer' / max(`cs_denom', 0.001) ///
                if `touse' & !missing(`tune_pred')
            quietly gen double `sq_err' = (`pseudo_out' - `tune_pred')^2 ///
                if `touse' & !missing(`tune_pred') & !missing(`pseudo_out')
            quietly summarize `sq_err' if `touse' & !missing(`sq_err')
            local cand_mse = r(mean)
            capture drop `sq_err' `pseudo_out'
        }
        else if "`foresttype'" == "boosted_regression" {
            /* Boosted regression: use single-step (booststeps=1) for tuning */
            capture plugin call grf_plugin `indepvars' `depvar' `extra_vars' `tune_pred' ///
                if `touse',                                                  ///
                "boosted_regression"                                         ///
                "`tunetrees'"                                                ///
                "`seed'"                                                     ///
                "`cand_mtry'"                                                ///
                "`cand_minnodesize'"                                         ///
                "`cand_samplefrac'"                                          ///
                "`do_honesty'"                                               ///
                "`cand_honestyfrac'"                                         ///
                "`do_honesty_prune'"                                         ///
                "`cand_alpha'"                                               ///
                "`cand_imbpen'"                                              ///
                "1"                                                          ///
                "`numthreads'"                                               ///
                "0"                                                          ///
                "0"                                                          ///
                "`nindep'"                                                   ///
                "1"                                                          ///
                "0"                                                          ///
                "0"                                                          ///
                "1"                                                          ///
                "`allow_missing_x'"                                          ///
                "`_main_cluster_idx'"                                        ///
                "`_main_weight_idx'"                                         ///
                "1"                                                          ///
                "0.97"                                                       ///
                "5"                                                          ///
                "10"                                                         ///
                "`do_stabilize'"

            if _rc {
                capture drop `tune_pred'
                continue
            }

            /* OOB MSE */
            tempvar sq_err
            quietly gen double `sq_err' = (`depvar' - `tune_pred')^2 if `touse' & !missing(`tune_pred')
            quietly summarize `sq_err' if `touse' & !missing(`tune_pred')
            local cand_mse = r(mean)
            capture drop `sq_err'
        }

        /* Drop temp prediction variable */
        capture drop `tune_pred'

        /* ---- Check if this is the best so far ---- */
        if !missing(`cand_mse') & (`cand_mse' < `best_mse' | `best_mse' == .) {
            local best_mse         = `cand_mse'
            local best_mtry        = `cand_mtry'
            local best_minnodesize = `cand_minnodesize'
            local best_samplefrac  = `cand_samplefrac'
            local best_honestyfrac = `cand_honestyfrac'
            local best_alpha       = `cand_alpha'
            local best_imbpen      = `cand_imbpen'
        }

        /* Progress indicator every 10 reps */
        if mod(`rep', 10) == 0 {
            display as text "  ... completed `rep'/`numreps' candidates (current best MSE = " ///
                as result %9.6f `best_mse' as text ")"
        }
    }

    /* ---- Check that we got at least one valid result ---- */
    if `best_mse' == . {
        display as error "all `numreps' candidate forests failed; cannot tune"
        exit 498
    }

    /* ---- Display results ---- */
    display as text ""
    display as text "Tuning Complete"
    display as text "{hline 55}"
    display as text "Best parameters from `numreps' random candidates:"
    display as text "{hline 55}"
    display as text %~30s "Parameter" %~20s "Value"
    display as text "{hline 55}"
    display as text %30s "mtry" _col(35) as result %12.0f `best_mtry'
    display as text %30s "min_node_size" _col(35) as result %12.0f `best_minnodesize'
    display as text %30s "sample_fraction" _col(35) as result %12.4f `best_samplefrac'
    display as text %30s "honesty_fraction" _col(35) as result %12.4f `best_honestyfrac'
    display as text %30s "alpha" _col(35) as result %12.4f `best_alpha'
    display as text %30s "imbalance_penalty" _col(35) as result %12.4f `best_imbpen'
    display as text "{hline 55}"
    display as text %30s "OOB MSE" _col(35) as result %12.6f `best_mse'
    display as text "{hline 55}"
    display as text ""

    /* ---- Construct usage hint ---- */
    local cmd_name "grf_`foresttype'_forest"
    if "`foresttype'" == "lm_forest" {
        local cmd_name "grf_lm_forest"
    }
    else if "`foresttype'" == "multi_arm_causal" {
        local cmd_name "grf_multi_arm_causal_forest"
    }
    display as text "Use these parameters in your forest call, e.g.:"
    display as text "  `cmd_name' ..., mtry(`best_mtry') minnodesize(`best_minnodesize') ///"
    display as text "    samplefrac(" %5.4f `best_samplefrac' ") honestyfrac(" %5.4f `best_honestyfrac' ") ///"
    display as text "    alpha(" %5.4f `best_alpha' ") imbalancepenalty(" %5.4f `best_imbpen' ")"
    display as text ""

    /* ---- Store results ---- */
    return scalar best_mtry             = `best_mtry'
    return scalar best_min_node_size    = `best_minnodesize'
    return scalar best_sample_fraction  = `best_samplefrac'
    return scalar best_honesty_fraction = `best_honestyfrac'
    return scalar best_alpha            = `best_alpha'
    return scalar best_imbalance_penalty = `best_imbpen'
    return scalar best_mse              = `best_mse'
    return scalar n_reps                = `numreps'
    return scalar N                     = `n_use'
    return local  forest_type             "`foresttype'"
end
