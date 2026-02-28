*! grf_multi_regression_forest.ado -- Multi-Output Regression Forest via grf C++ library
*! Version 0.1.0
*! Implements grf::multi_regression_forest()
*!
*! Jointly predicts multiple outcome variables using a single forest that
*! accounts for cross-outcome structure.

program define grf_multi_regression_forest, eclass
    version 14.0

    /* ---- Syntax ----
     * varlist: depvar1 depvar2 ... indepvar1 indepvar2 ...
     * ndep() specifies how many leading variables are outcomes.
     */
    syntax varlist(min=3 numeric) [if] [in],  ///
        GENerate(name)                         ///
        NDEp(integer)                          ///
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
            noMIA                              ///
            CLuster(varname numeric)           ///
            WEIghts(varname numeric)           ///
            EQUALizeclusterweights             ///
        ]

    /* ---- Validate ndep ---- */
    if `ndep' < 2 {
        display as error "ndep() must be at least 2 (use grf_regression_forest for single outcome)"
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

    /* ---- Parse estimate_variance ---- */
    /* NOTE: grf's multi_regression_forest does not support variance estimation.
     * If requested, warn and ignore. */
    local do_est_var 0
    if "`estimatevariance'" != "" {
        display as text "(warning: estimatevariance not supported for multi_regression_forest; ignored)"
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

    /* ---- Parse varlist: Y1 Y2 ... X1 X2 ... ---- */
    local depvars ""
    local rest `varlist'
    forvalues j = 1/`ndep' {
        gettoken dv rest : rest
        local depvars `depvars' `dv'
    }
    local indepvars `rest'
    local nindep : word count `indepvars'

    if `nindep' < 1 {
        display as error "need at least 1 predictor variable"
        exit 198
    }

    /* ---- Handle replace ---- */
    if "`replace'" != "" {
        forvalues j = 1/`ndep' {
            capture drop `generate'_y`j'
            if `do_est_var' {
                capture drop `generate'_y`j'_var
            }
        }
    }

    /* ---- Mark sample ---- */
    if `allow_missing_x' {
        marksample touse, novarlist
        foreach dv of local depvars {
            markout `touse' `dv'
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
    /* One prediction per outcome variable */
    local n_output `ndep'
    if `do_est_var' {
        local n_output = `ndep' * 2
    }

    forvalues j = 1/`ndep' {
        confirm new variable `generate'_y`j'
        quietly gen double `generate'_y`j' = .
        if `do_est_var' {
            confirm new variable `generate'_y`j'_var
            quietly gen double `generate'_y`j'_var = .
        }
    }

    /* ---- Build output varlist ----
     * Plugin writes contiguously: all predictions first, then all variances.
     * So output_vars must be: pred_y1 pred_y2 ... var_y1 var_y2 ...
     */
    local output_vars ""
    forvalues j = 1/`ndep' {
        local output_vars `output_vars' `generate'_y`j'
    }
    if `do_est_var' {
        forvalues j = 1/`ndep' {
            local output_vars `output_vars' `generate'_y`j'_var
        }
    }

    /* ---- Display header ---- */
    display as text ""
    display as text "Generalized Random Forest: Multi-Output Regression Forest"
    display as text "{hline 60}"
    display as text "Outcome variables:     " as result "`depvars'"
    display as text "  Number of outcomes:  " as result `ndep'
    display as text "Predictors:            " as result "`indepvars'"
    display as text "Observations:          " as result `n_use'
    display as text "Trees:                 " as result `ntrees'
    display as text "Honesty:               " as result cond(`do_honesty', "yes", "no")
    if `do_est_var' {
        display as text "Variance estimation:   " as result "yes (ci_group_size=`cigroupsize')"
    }
    display as text "{hline 60}"
    display as text ""

    /* ---- Load plugin ---- */
    if ( inlist("`c(os)'", "MacOSX") | strpos("`c(machine_type)'", "Mac") ) local c_os_ macosx
    else local c_os_: di lower("`c(os)'")

    capture program grf_plugin, plugin using("grf_plugin_`c_os_'.plugin")

    /* ---- Build extra vars for cluster/weight ---- */
    local extra_vars ""
    local n_data_before = `nindep' + `ndep'
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
     * Variable order: X1..Xp Y1 Y2 ... [cluster] [weight] out_y1 [out_y1_var] out_y2 [out_y2_var] ...
     *   n_x = nindep (covariates)
     *   n_y = ndep (multiple outcome columns)
     *   n_w = 0
     *   n_z = 0
     *   n_output = ndep (or ndep*2 with variance)
     *
     * argv[20-22]: allow_missing_x, cluster_col_idx, weight_col_idx
     * argv[23]: num_outcomes
     */
    display as text "Fitting multi-output regression forest ..."
    plugin call grf_plugin `indepvars' `depvars' `extra_vars' `output_vars' ///
        if `touse',                                             ///
        "multi_regression"                                      ///
        "`ntrees'"                                              ///
        "`seed'"                                                ///
        "`mtry'"                                                ///
        "`minnodesize'"                                         ///
        "`samplefrac'"                                          ///
        "`do_honesty'"                                          ///
        "`honestyfrac'"                                         ///
        "`do_honesty_prune'"                                    ///
        "`alpha'"                                               ///
        "`imbalancepenalty'"                                     ///
        "`cigroupsize'"                                         ///
        "`numthreads'"                                          ///
        "`do_est_var'"                                          ///
        "0"                                                     ///
        "`nindep'"                                              ///
        "`ndep'"                                                ///
        "0"                                                     ///
        "0"                                                     ///
        "`n_output'"                                            ///
        "`allow_missing_x'"                                     ///
        "`cluster_col_idx'"                                     ///
        "`weight_col_idx'"                                      ///
        "`ndep'"

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
    ereturn scalar n_outcomes  = `ndep'
    ereturn local  cmd           "grf_multi_regression_forest"
    ereturn scalar allow_missing_x = `allow_missing_x'
    ereturn local  forest_type   "multi_regression"
    ereturn local  depvars       "`depvars'"
    ereturn local  indepvars     "`indepvars'"
    ereturn local  predict_stub  "`generate'"
    if "`cluster_var'" != "" {
        ereturn local cluster_var "`cluster_var'"
    }
    if "`weight_var'" != "" {
        ereturn local weight_var "`weight_var'"
    }

    /* ---- Summary stats ---- */
    display as text ""
    display as text "Multi-Output Regression Forest Results"
    display as text "{hline 60}"

    forvalues j = 1/`ndep' {
        local dv : word `j' of `depvars'
        quietly summarize `generate'_y`j' if `touse'
        local n_pred = r(N)
        local pred_mean = r(mean)
        local pred_sd = r(sd)

        display as text ""
        display as text "Outcome `j' (`dv'):"
        display as text "  Variable:     " as result "`generate'_y`j'"
        display as text "  Non-missing:  " as result `n_pred'
        display as text "  Mean:         " as result %9.4f `pred_mean'
        display as text "  SD:           " as result %9.4f `pred_sd'
        if `do_est_var' {
            quietly summarize `generate'_y`j'_var if `touse'
            display as text "  Mean variance: " as result %9.6f r(mean)
        }
    }

    display as text ""
    display as text "{hline 60}"
    display as text ""
end
