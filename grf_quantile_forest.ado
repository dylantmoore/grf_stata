*! grf_quantile_forest.ado -- Quantile Forest via grf C++ library
*! Version 0.1.0
*! Implements Athey, Tibshirani, Wager (2019) quantile_forest()

program define grf_quantile_forest, eclass
    version 14.0

    syntax varlist(min=2 numeric) [if] [in],  ///
        GENerate(name)                         ///
        [                                      ///
            QUANTiles(string)                  ///
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
            NUMThreads(integer 0)              ///
            REPlace                            ///
        ]

    /* ---- Default quantiles ---- */
    if "`quantiles'" == "" {
        local quantiles "0.1 0.5 0.9"
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

    /* ---- Parse quantile list and build output variable names ---- */
    local n_quantiles 0
    local output_vars ""
    local quantile_csv ""
    local quantile_display ""

    foreach q of local quantiles {
        /* Validate quantile is in (0,1) */
        if `q' <= 0 | `q' >= 1 {
            display as error "quantiles must be between 0 and 1 (exclusive)"
            exit 198
        }

        local n_quantiles = `n_quantiles' + 1

        /* Build variable name: stub_qNN where NN = quantile*100 */
        local qint = round(`q' * 100)
        local varname `generate'_q`qint'

        /* Handle replace */
        if "`replace'" != "" {
            capture drop `varname'
        }
        confirm new variable `varname'

        local output_vars `output_vars' `varname'

        /* Build comma-separated string for plugin */
        if `n_quantiles' == 1 {
            local quantile_csv "`q'"
        }
        else {
            local quantile_csv "`quantile_csv',`q'"
        }

        /* Build display string */
        if `n_quantiles' == 1 {
            local quantile_display "`q'"
        }
        else {
            local quantile_display "`quantile_display' `q'"
        }
    }

    if `n_quantiles' < 1 {
        display as error "must specify at least one quantile"
        exit 198
    }

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

    /* ---- Create output variables ---- */
    foreach v of local output_vars {
        quietly gen double `v' = .
    }

    /* ---- Display header ---- */
    display as text ""
    display as text "Generalized Random Forest: Quantile Forest"
    display as text "{hline 55}"
    display as text "Dependent variable:    " as result "`depvar'"
    display as text "Predictors:            " as result "`indepvars'"
    display as text "Observations:          " as result `n_use'
    display as text "Trees:                 " as result `ntrees'
    display as text "Quantiles:             " as result "`quantile_display'"
    display as text "Honesty:               " as result cond(`do_honesty', "yes", "no")
    display as text "{hline 55}"
    display as text ""

    /* ---- Load plugin ---- */
    if ( inlist("`c(os)'", "MacOSX") | strpos("`c(machine_type)'", "Mac") ) local c_os_ macosx
    else local c_os_: di lower("`c(os)'")

    cap program drop grf_plugin
    program grf_plugin, plugin using("grf_plugin_`c_os_'.plugin")

    /* ---- Call plugin ----
     *
     * Variable order: X1..Xp Y out1 out2 ... out_nq
     * argv: forest_type num_trees seed mtry min_node_size sample_fraction
     *       honesty honesty_fraction honesty_prune alpha imbalance_penalty
     *       ci_group_size num_threads estimate_variance compute_oob
     *       n_x n_y n_w n_z n_output quantiles_csv
     */
    plugin call grf_plugin `indepvars' `depvar' `output_vars' ///
        if `touse',                                            ///
        "quantile"                                             ///
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
        "1"                                                    ///
        "`numthreads'"                                         ///
        "0"                                                    ///
        "0"                                                    ///
        "`nindep'"                                             ///
        "1"                                                    ///
        "0"                                                    ///
        "0"                                                    ///
        "`n_quantiles'"                                        ///
        "`quantile_csv'"

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
    ereturn scalar ci_group_size      = 1
    ereturn scalar n_quantiles = `n_quantiles'
    ereturn local  cmd           "grf_quantile_forest"
    ereturn local  forest_type   "quantile"
    ereturn local  depvar        "`depvar'"
    ereturn local  indepvars     "`indepvars'"
    ereturn local  quantiles     "`quantile_display'"
    ereturn local  predict_vars  "`output_vars'"

    /* ---- Summary stats ---- */
    display as text ""
    display as text "Predictions written to:"
    local qi 0
    foreach q of local quantiles {
        local qi = `qi' + 1
        local qint = round(`q' * 100)
        local v `generate'_q`qint'
        quietly summarize `v' if `touse'
        display as text "  `v'" _col(30) as text "q(`q')" ///
            _col(42) as text "mean=" as result %9.4f r(mean) ///
            _col(60) as text "sd=" as result %9.4f r(sd)
    }
    display as text ""
end
