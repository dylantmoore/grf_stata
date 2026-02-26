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
    marksample touse
    quietly count if `touse'
    local n_use = r(N)

    if `n_use' < 2 {
        display as error "need at least 2 non-missing observations"
        exit 2000
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

    /* ---- Call plugin ----
     *
     * Variable order: X1..Xp Y1 Y2 ... out_y1 [out_y1_var] out_y2 [out_y2_var] ...
     *   n_x = nindep (covariates)
     *   n_y = ndep (multiple outcome columns)
     *   n_w = 0
     *   n_z = 0
     *   n_output = ndep (or ndep*2 with variance)
     *
     * argv[20]: num_outcomes
     */
    display as text "Fitting multi-output regression forest ..."
    plugin call grf_plugin `indepvars' `depvars' `output_vars' ///
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
    ereturn local  forest_type   "multi_regression"
    ereturn local  depvars       "`depvars'"
    ereturn local  indepvars     "`indepvars'"
    ereturn local  predict_stub  "`generate'"

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
