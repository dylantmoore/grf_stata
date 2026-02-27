*! grf_tune.ado -- Cross-validation tuning for GRF forests
*! Version 0.1.0
*! Random search over hyperparameters with OOB MSE evaluation

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
        ]

    /* ---- Validate forest_type ---- */
    if "`foresttype'" != "regression" & "`foresttype'" != "causal" & "`foresttype'" != "quantile" {
        display as error "forest_type() must be regression, causal, or quantile"
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

    /* ---- Parse stabilize_splits (causal only) ---- */
    local do_stabilize 1
    if "`stabilizesplits'" == "nostabilizesplits" {
        local do_stabilize 0
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
    else {
        /* depvar indepvars */
        gettoken depvar indepvars : varlist
        local nindep : word count `indepvars'
        if `nindep' < 1 {
            display as error "need at least 1 predictor variable"
            exit 198
        }
    }

    /* ---- Mark sample ---- */
    marksample touse
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

    /* ---- Display header ---- */
    display as text ""
    display as text "GRF Hyperparameter Tuning (Random Search)"
    display as text "{hline 55}"
    display as text "Forest type:           " as result "`foresttype'"
    if "`foresttype'" == "causal" {
        display as text "Outcome variable:      " as result "`depvar'"
        display as text "Treatment variable:    " as result "`treatvar'"
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

    cap program drop grf_plugin
    program grf_plugin, plugin using("grf_plugin_`c_os_'.plugin")

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

    /* ---- For causal forest: fit nuisance models once ---- */
    if "`foresttype'" == "causal" {
        display as text "Fitting nuisance models (Y ~ X and W ~ X) ..."
        tempvar yhat_nuis what_nuis
        quietly gen double `yhat_nuis' = .
        quietly gen double `what_nuis' = .

        /* Y ~ X regression forest (small, for tuning) */
        plugin call grf_plugin `indepvars' `depvar' `yhat_nuis' ///
            if `touse',                                          ///
            "regression"                                         ///
            "`tunetrees'"                                     ///
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
            "1"

        /* W ~ X regression forest */
        plugin call grf_plugin `indepvars' `treatvar' `what_nuis' ///
            if `touse',                                            ///
            "regression"                                           ///
            "`tunetrees'"                                       ///
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
            "1"

        /* Center Y and W */
        tempvar y_centered w_centered
        quietly gen double `y_centered' = `depvar' - `yhat_nuis' if `touse'
        quietly gen double `w_centered' = `treatvar' - `what_nuis' if `touse'

        display as text "Nuisance models fitted. Beginning random search ..."
        display as text ""
    }
    else {
        display as text "Beginning random search ..."
        display as text ""
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
            capture plugin call grf_plugin `indepvars' `depvar' `tune_pred' ///
                if `touse',                                                  ///
                "regression"                                                 ///
                "`tunetrees'"                                             ///
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
                "1"

            if _rc {
                capture drop `tune_pred'
                continue
            }

            /* Compute OOB MSE: mean((Y - Yhat)^2) */
            tempvar sq_err
            quietly gen double `sq_err' = (`depvar' - `tune_pred')^2 if `touse' & !missing(`tune_pred')
            quietly summarize `sq_err' if `touse' & !missing(`tune_pred')
            local cand_mse = r(mean)
            capture drop `sq_err'
        }
        else if "`foresttype'" == "causal" {
            capture plugin call grf_plugin `indepvars' `y_centered' `w_centered' `tune_pred' ///
                if `touse',                                                                   ///
                "causal"                                                                      ///
                "`tunetrees'"                                                              ///
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
                "`do_stabilize'"

            if _rc {
                capture drop `tune_pred'
                continue
            }

            /* Compute OOB MSE for causal: mean((Y.centered - tau_hat * W.centered)^2) */
            tempvar sq_err
            quietly gen double `sq_err' = (`y_centered' - `tune_pred' * `w_centered')^2 ///
                if `touse' & !missing(`tune_pred')
            quietly summarize `sq_err' if `touse' & !missing(`tune_pred')
            local cand_mse = r(mean)
            capture drop `sq_err'
        }
        else if "`foresttype'" == "quantile" {
            /* For quantile forest, tune on median prediction MSE */
            capture plugin call grf_plugin `indepvars' `depvar' `tune_pred' ///
                if `touse',                                                  ///
                "quantile"                                                   ///
                "`tunetrees'"                                             ///
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
                "0.5"

            if _rc {
                capture drop `tune_pred'
                continue
            }

            /* Compute OOB MSE: mean((Y - median_hat)^2) */
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
    display as text "Use these parameters in your forest call, e.g.:"
    display as text "  grf_`foresttype'_forest ..., mtry(`best_mtry') minnodesize(`best_minnodesize') ///"
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
