*! grf_survival_forest.ado -- Survival Forest via grf C++ library
*! Version 0.1.0
*! Implements Cui et al. (2023) survival_forest()

program define grf_survival_forest, eclass
    version 14.0

    syntax varlist(min=3 numeric) [if] [in],  ///
        GENerate(name)                         ///
        [                                      ///
            NTrees(integer 1000)               ///
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
            noMIA                              ///
            noFASTlogrank                      ///
            CLuster(varname numeric)           ///
            WEIghts(varname numeric)           ///
            EQUALizeclusterweights             ///
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

    /* ---- Parse MIA ---- */
    local allow_missing_x 1
    if "`mia'" == "nomia" {
        local allow_missing_x 0
    }

    /* ---- Parse fast_logrank ---- */
    /* R defaults fast.logrank = TRUE. noFASTlogrank disables it. */
    local do_fast_logrank 1
    if "`fastlogrank'" == "nofastlogrank" {
        local do_fast_logrank 0
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

    /* ---- Parse varlist: timevar statusvar indepvars ---- */
    gettoken timevar rest : varlist
    gettoken statusvar indepvars : rest
    local nindep : word count `indepvars'

    if `nindep' < 1 {
        display as error "need at least 1 predictor variable"
        exit 198
    }

    /* ---- Mark sample ---- */
    if `allow_missing_x' {
        marksample touse, novarlist
        markout `touse' `timevar'
        markout `touse' `statusvar'
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
    if ( inlist("`c(os)'", "MacOSX") | strpos("`c(machine_type)'", "Mac") ) local c_os_ macosx
    else local c_os_: di lower("`c(os)'")

    capture program grf_plugin, plugin using("grf_plugin_`c_os_'.plugin")

    /* ---- Build plugin varlist with optional cluster/weight columns ---- */
    local extra_vars ""
    if "`cluster_var'" != "" {
        local extra_vars `extra_vars' `cluster_var'
    }
    if "`weight_var'" != "" {
        local extra_vars `extra_vars' `weight_var'
    }

    /* Compute cluster/weight column indices in the plugin varlist */
    /* Survival forest data vars: X1..Xp time status => nindep + 2 columns */
    local _data_col_count = `nindep' + 2
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

    /* ---- Call plugin ----
     *
     * Variable order: X1..Xp time status [cluster] [weight] out1..outN
     *   n_x = nindep (covariates)
     *   n_y = 1 (survival time)
     *   n_w = 0 (censor column handled via set_censor_index, not as treatment)
     *   n_z = 0
     *   n_output = noutput (survival curve columns)
     *
     * argv[0..19]: standard forest params
     * argv[20..22]: allow_missing_x cluster_col_idx weight_col_idx
     * argv[23]: num_failures (0 = auto-detect)
     * argv[24]: prediction_type (C++: 0 = Kaplan-Meier, 1 = Nelson-Aalen)
     *           User-facing: predtype(0) = Nelson-Aalen, predtype(1) = Kaplan-Meier
     */
    /* Remap user-facing predtype to C++ convention (inverted) */
    local cpp_predtype = 1 - `predtype'
    display as text "Fitting survival forest ..."
    plugin call grf_plugin `indepvars' `timevar' `statusvar' `extra_vars' `output_vars' ///
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
        "`allow_missing_x'"                                                 ///
        "`cluster_col_idx'"                                                 ///
        "`weight_col_idx'"                                                  ///
        "`numfailures'"                                                     ///
        "`cpp_predtype'"                                                    ///
        "`do_fast_logrank'"

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
    ereturn scalar allow_missing_x = `allow_missing_x'
    ereturn local  forest_type   "survival"
    ereturn local  timevar       "`timevar'"
    ereturn local  statusvar     "`statusvar'"
    ereturn local  indepvars     "`indepvars'"
    ereturn local  predict_stub  "`generate'"
    if "`cluster_var'" != "" {
        ereturn local cluster_var "`cluster_var'"
    }
    if "`weight_var'" != "" {
        ereturn local weight_var "`weight_var'"
    }

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
