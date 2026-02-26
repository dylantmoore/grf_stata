*! grf_survival_forest.ado -- Survival Forest via grf C++ library
*! Version 0.1.0
*! Implements Cui et al. (2023) survival_forest()

program define grf_survival_forest, eclass
    version 14.0

    syntax varlist(min=3 numeric) [if] [in],  ///
        GENerate(name)                         ///
        [                                      ///
            NTrees(integer 2000)               ///
            SEED(integer 42)                   ///
            MTRY(integer 0)                    ///
            MINNodesize(integer 15)            ///
            SAMPLEfrac(real 0.5)               ///
            noHONesty                          ///
            HONestyfrac(real 0.5)              ///
            noHONestyprune                     ///
            ALPha(real 0.05)                   ///
            IMBalancepenalty(real 0.0)          ///
            CIGroupsize(integer 1)             ///
            NUMThreads(integer 0)              ///
            NOUTput(integer 20)                ///
            NUMFailures(integer 0)             ///
            PREDtype(integer 1)                ///
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

    /* ---- Validate noutput ---- */
    if `noutput' < 1 {
        display as error "noutput() must be at least 1"
        exit 198
    }

    /* ---- Validate predtype ---- */
    if `predtype' < 0 | `predtype' > 1 {
        display as error "predtype() must be 0 (Nelson-Aalen) or 1 (Kaplan-Meier)"
        exit 198
    }

    /* ---- Handle replace ---- */
    if "`replace'" != "" {
        forvalues j = 1/`noutput' {
            capture drop `generate'_s`j'
        }
    }

    /* ---- Parse varlist: timevar statusvar indepvars ---- */
    gettoken timevar rest : varlist
    gettoken statusvar indepvars : rest
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

    /* ---- Validate survival time (must be positive) ---- */
    quietly summarize `timevar' if `touse'
    if r(min) <= 0 {
        display as error "survival time variable `timevar' must be strictly positive"
        exit 198
    }

    /* ---- Validate status variable (0/1) ---- */
    quietly summarize `statusvar' if `touse'
    if r(min) < 0 | r(max) > 1 {
        display as error "status variable `statusvar' must be 0 (censored) or 1 (event)"
        exit 198
    }
    local n_events = 0
    quietly count if `statusvar' == 1 & `touse'
    local n_events = r(N)
    local n_censored = `n_use' - `n_events'

    /* ---- Create output variables ---- */
    forvalues j = 1/`noutput' {
        confirm new variable `generate'_s`j'
        quietly gen double `generate'_s`j' = .
    }

    /* ---- Build output varlist ---- */
    local output_vars ""
    forvalues j = 1/`noutput' {
        local output_vars `output_vars' `generate'_s`j'
    }

    /* ---- Display header ---- */
    display as text ""
    display as text "Generalized Random Forest: Survival Forest"
    display as text "{hline 55}"
    display as text "Time variable:         " as result "`timevar'"
    display as text "Status variable:       " as result "`statusvar'"
    display as text "Predictors:            " as result "`indepvars'"
    display as text "Observations:          " as result `n_use'
    display as text "  Events:              " as result `n_events'
    display as text "  Censored:            " as result `n_censored'
    display as text "Trees:                 " as result `ntrees'
    display as text "Honesty:               " as result cond(`do_honesty', "yes", "no")
    display as text "Output columns:        " as result `noutput'
    local predlabel = cond(`predtype' == 0, "Nelson-Aalen", "Kaplan-Meier")
    display as text "Prediction type:       " as result "`predlabel'"
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
     * Variable order: X1..Xp time status out1..outN
     *   n_x = nindep (covariates)
     *   n_y = 1 (survival time)
     *   n_w = 0 (censor column handled via set_censor_index, not as treatment)
     *   n_z = 0
     *   n_output = noutput (survival curve columns)
     *
     * argv[0..19]: standard forest params
     * argv[20]: num_failures (0 = auto-detect)
     * argv[21]: prediction_type (C++: 0 = Kaplan-Meier, 1 = Nelson-Aalen)
     *           User-facing: predtype(0) = Nelson-Aalen, predtype(1) = Kaplan-Meier
     */
    /* Remap user-facing predtype to C++ convention (inverted) */
    local cpp_predtype = 1 - `predtype'
    display as text "Fitting survival forest ..."
    plugin call grf_plugin `indepvars' `timevar' `statusvar' `output_vars' ///
        if `touse',                                                         ///
        "survival"                                                          ///
        "`ntrees'"                                                          ///
        "`seed'"                                                            ///
        "`mtry'"                                                            ///
        "`minnodesize'"                                                     ///
        "`samplefrac'"                                                      ///
        "`do_honesty'"                                                      ///
        "`honestyfrac'"                                                     ///
        "`do_honesty_prune'"                                                ///
        "`alpha'"                                                           ///
        "`imbalancepenalty'"                                                 ///
        "`cigroupsize'"                                                     ///
        "`numthreads'"                                                      ///
        "0"                                                                 ///
        "0"                                                                 ///
        "`nindep'"                                                          ///
        "1"                                                                 ///
        "0"                                                                 ///
        "0"                                                                 ///
        "`noutput'"                                                         ///
        "`numfailures'"                                                     ///
        "`cpp_predtype'"

    /* ---- Store results ---- */
    ereturn clear
    ereturn scalar N           = `n_use'
    ereturn scalar n_events    = `n_events'
    ereturn scalar n_censored  = `n_censored'
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
    ereturn scalar n_output    = `noutput'
    ereturn scalar pred_type   = `predtype'
    ereturn local  cmd           "grf_survival_forest"
    ereturn local  forest_type   "survival"
    ereturn local  timevar       "`timevar'"
    ereturn local  statusvar     "`statusvar'"
    ereturn local  indepvars     "`indepvars'"
    ereturn local  predict_stub  "`generate'"

    /* ---- Summary stats ---- */
    /* Report stats from first output column as a representative summary */
    quietly summarize `generate'_s1 if `touse'
    local n_pred = r(N)
    local pred_mean = r(mean)
    local pred_sd = r(sd)

    display as text ""
    display as text "Survival Forest Results"
    display as text "{hline 55}"
    display as text "Predictions written to: " as result "`generate'_s1 ... `generate'_s`noutput'"
    display as text "  (`noutput' columns at evaluated failure times)"
    display as text ""
    display as text "First column (`generate'_s1):"
    display as text "  Non-missing:  " as result `n_pred'
    display as text "  Mean:         " as result %9.4f `pred_mean'
    display as text "  SD:           " as result %9.4f `pred_sd'
    display as text "{hline 55}"
    display as text ""
end
