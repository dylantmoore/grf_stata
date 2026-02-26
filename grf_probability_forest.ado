*! grf_probability_forest.ado -- Probability Forest via grf C++ library
*! Version 0.1.0
*! Implements Athey, Tibshirani, Wager (2019) probability_forest()

program define grf_probability_forest, eclass
    version 14.0

    syntax varlist(min=2 numeric) [if] [in],  ///
        GENerate(name)                         ///
        [                                      ///
            NCLasses(integer 0)                ///
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

    /* ---- Validate depvar is integer-valued ---- */
    tempvar is_int
    quietly gen byte `is_int' = (`depvar' == floor(`depvar')) if `touse'
    quietly count if `is_int' == 0 & `touse'
    if r(N) > 0 {
        display as error "`depvar' must be integer-valued (0, 1, 2, ...)"
        exit 198
    }

    quietly summarize `depvar' if `touse'
    local y_min = r(min)
    local y_max = r(max)

    if `y_min' < 0 {
        display as error "`depvar' must be non-negative (0, 1, 2, ...)"
        exit 198
    }

    /* ---- Auto-detect or validate nclasses ---- */
    if `nclasses' == 0 {
        local nclasses = `y_max' + 1
    }
    else {
        if `nclasses' <= `y_max' {
            display as error "nclasses(`nclasses') must be > max(`depvar') = `y_max'"
            exit 198
        }
    }

    if `nclasses' < 2 {
        display as error "need at least 2 classes"
        exit 198
    }

    /* ---- Create output variables: stub_c0, stub_c1, ... ---- */
    local output_vars ""
    forvalues c = 0/`=`nclasses'-1' {
        local varname `generate'_c`c'
        if "`replace'" != "" {
            capture drop `varname'
        }
        confirm new variable `varname'
        quietly gen double `varname' = .
        local output_vars `output_vars' `varname'
    }

    /* ---- Display header ---- */
    display as text ""
    display as text "Generalized Random Forest: Probability Forest"
    display as text "{hline 55}"
    display as text "Dependent variable:    " as result "`depvar'"
    display as text "Predictors:            " as result "`indepvars'"
    display as text "Observations:          " as result `n_use'
    display as text "Classes:               " as result `nclasses' as text " (0 to `=`nclasses'-1')"
    display as text "Trees:                 " as result `ntrees'
    display as text "Honesty:               " as result cond(`do_honesty', "yes", "no")
    display as text "{hline 55}"
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

    /* ---- Call plugin ----
     *
     * Variable order: X1..Xp Y out1 out2 ... out_nclasses
     * argv: forest_type num_trees seed mtry min_node_size sample_fraction
     *       honesty honesty_fraction honesty_prune alpha imbalance_penalty
     *       ci_group_size num_threads estimate_variance compute_oob
     *       n_x n_y n_w n_z n_output num_classes
     */
    plugin call grf_plugin `indepvars' `depvar' `output_vars' ///
        if `touse',                                            ///
        "probability"                                          ///
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
        "`nclasses'"                                           ///
        "`nclasses'"

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
    ereturn scalar n_classes   = `nclasses'
    ereturn local  cmd           "grf_probability_forest"
    ereturn local  forest_type   "probability"
    ereturn local  depvar        "`depvar'"
    ereturn local  indepvars     "`indepvars'"
    ereturn local  predict_vars  "`output_vars'"

    /* ---- Summary stats ---- */
    display as text ""
    display as text "Class probability estimates written to:"
    forvalues c = 0/`=`nclasses'-1' {
        local v `generate'_c`c'
        quietly summarize `v' if `touse'
        display as text "  `v'" _col(30) as text "P(Y=`c')" ///
            _col(42) as text "mean=" as result %7.4f r(mean) ///
            _col(58) as text "sd=" as result %7.4f r(sd)
    }
    display as text ""
end
