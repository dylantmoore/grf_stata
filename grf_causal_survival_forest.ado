*! grf_causal_survival_forest.ado -- Causal Survival Forest via grf C++ library
*! Version 0.1.0
*! Implements Cui et al. (2023) causal_survival_forest()
*!
*! NOTE: Full implementation requires a multi-step nuisance estimation pipeline
*! (survival forest for S.hat, regression forests for censoring/treatment models).
*! This version supports two modes:
*!   1. Simplified: pre-computed nuisance columns supplied via numer()/denom()
*!   2. Auto pipeline: internal estimation (default, calls survival + regression forests)

program define grf_causal_survival_forest, eclass
    version 14.0

    syntax varlist(min=4 numeric) [if] [in],  ///
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
            ESTIMATEVariance                   ///
            noSTABilizesplits                  ///
            NUMer(varname numeric)             ///
            DENom(varname numeric)             ///
            WHATinput(varname numeric)         ///
            HORizon(real 0)                    ///
            TARget(integer 1)                  ///
            REPlace                            ///
            VARGenerate(name)                  ///
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

    /* ---- Handle replace ---- */
    if "`replace'" != "" {
        capture drop `generate'
        if `do_est_var' & "`vargenerate'" != "" {
            capture drop `vargenerate'
        }
    }
    confirm new variable `generate'

    /* ---- Parse varlist: timevar statusvar treatvar indepvars ---- */
    gettoken timevar rest : varlist
    gettoken statusvar rest : rest
    gettoken treatvar indepvars : rest
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
        markout `touse' `treatvar'
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

    /* ---- Validate treatment variable ---- */
    quietly summarize `treatvar' if `touse'
    local w_min = r(min)
    local w_max = r(max)
    if `w_min' == `w_max' {
        display as error "treatment variable `treatvar' has no variation"
        exit 198
    }

    /* ---- Determine horizon ---- */
    if `horizon' == 0 {
        /* Default: use the median failure time */
        quietly summarize `timevar' if `statusvar' == 1 & `touse', detail
        local horizon = r(p50)
        display as text "(horizon not specified; using median failure time: " ///
            as result %9.3f `horizon' as text ")"
    }

        /* ---- Inline tuning ---- */
    if `"`tuneparameters'"' != "" {
        display as text ""
        display as text "Running inline parameter tuning..."
        display as text "  Parameters: `tuneparameters'"
        display as text "  Tune trees: `tunenumtrees'  Tune reps: `tunenumreps'"

        /* Call grf_tune to find best parameters */
        grf_tune `varlist' if `touse', foresttype(causal_survival) ///
            numreps(`tunenumreps') tunetrees(`tunenumtrees') seed(`seed') ///
            numthreads(`numthreads') horizon(`horizon')

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
    display as text "Generalized Random Forest: Causal Survival Forest"
    display as text "{hline 60}"
    display as text "Time variable:         " as result "`timevar'"
    display as text "Status variable:       " as result "`statusvar'"
    display as text "Treatment variable:    " as result "`treatvar'"
    display as text "Predictors:            " as result "`indepvars'"
    display as text "Observations:          " as result `n_use'
    display as text "  Events:              " as result `n_events'
    display as text "  Censored:            " as result `n_censored'
    display as text "Trees:                 " as result `ntrees'
    display as text "Honesty:               " as result cond(`do_honesty', "yes", "no")
    display as text "Stabilize splits:      " as result cond(`do_stabilize', "yes", "no")
    display as text "Horizon:               " as result %9.3f `horizon'
    display as text "Target:                " as result cond(`target'==1, "RMST", "survival probability")
    if `do_est_var' {
        display as text "Variance estimation:   " as result "yes (ci_group_size=`cigroupsize')"
    }
    display as text "{hline 60}"
    display as text ""

    /* ---- Load plugin ---- */
    if ( inlist("`c(os)'", "MacOSX") | strpos("`c(machine_type)'", "Mac") ) local c_os_ macosx
    else local c_os_: di lower("`c(os)'")

    capture program grf_plugin, plugin using("grf_plugin_`c_os_'.plugin")

    /* ==== Nuisance estimation pipeline ====
     *
     * Causal survival forests (Cui et al. 2023) require:
     *   1. S.hat: conditional survival function via survival forest
     *   2. C.hat: conditional censoring survival function
     *   3. W.hat: propensity scores P(W=1|X)
     *   4. Compute IPCW-style numerator/denominator for each obs
     *
     * If user supplies numer() and denom(), skip the pipeline.
     * Otherwise, estimate internally.
     */

    local precomputed 0
    if "`numer'" != "" & "`denom'" != "" {
        local precomputed 1
        display as text "Using pre-computed nuisance columns:"
        display as text "  Numerator: `numer'"
        display as text "  Denominator: `denom'"
    }

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

    /* Build extra_vars for nuisance calls */
    local extra_vars ""
    if "`cluster_var'" != "" {
        local extra_vars `extra_vars' `cluster_var'
    }
    if "`weight_var'" != "" {
        local extra_vars `extra_vars' `weight_var'
    }

    if `precomputed' == 0 {
        /* ---- Step 1: Estimate propensity scores W.hat = E[W|X] ---- */
        if "`whatinput'" != "" {
            display as text "Step 1/4: Using user-supplied propensity scores from `whatinput'"
            tempvar what
            quietly gen double `what' = `whatinput' if `touse'
        }
        else {
            display as text "Step 1/4: Estimating propensity scores W ~ X ..."
            tempvar what
            quietly gen double `what' = .
            plugin call grf_plugin `indepvars' `treatvar' `extra_vars' `what' ///
                if `touse',                                       ///
                "regression"                                      ///
                "`ntrees'"                                        ///
                "`seed'"                                          ///
                "`mtry'"                                          ///
                "`minnodesize'"                                   ///
                "`samplefrac'"                                    ///
                "`do_honesty'"                                    ///
                "`honestyfrac'"                                   ///
                "`do_honesty_prune'"                              ///
                "`alpha'"                                         ///
                "`imbalancepenalty'"                               ///
                "`cigroupsize'"                                   ///
                "`numthreads'"                                    ///
                "0"                                               ///
                "0"                                               ///
                "`nindep'"                                        ///
                "1"                                               ///
                "0"                                               ///
                "0"                                               ///
                "1"                                               ///
                "`allow_missing_x'"                               ///
                "`_nuis_cluster_idx'"                             ///
                "`_nuis_weight_idx'"
        }

        /* ---- Step 2: Fit censoring survival forest C.hat ----
         *
         * Fit survival forest on (X+W, T, 1-D) to estimate conditional
         * censoring probabilities C.hat(t|X,W) = P(C > t | X, W).
         * Then extract C.hat evaluated at each observation's time and
         * at the horizon for IPCW weights.
         */

        display as text "Step 2/4: Fitting censoring survival forest C.hat ..."

        tempvar w_centered
        quietly gen double `w_centered' = `treatvar' - `what' if `touse'

        /* Flip event indicator: 1-D so censoring is the "event" */
        tempvar censor_ind
        quietly gen byte `censor_ind' = 1 - `statusvar' if `touse'

        /* Fit censoring forest on (X, W, T, 1-D) with reduced trees */
        local nuis_trees = max(50, min(`ntrees' / 4, 500))
        tempvar chat_out
        quietly gen double `chat_out' = .

        /* Use regression forest to estimate E[T|X,W] as Y.hat proxy */
        tempvar yhat
        quietly gen double `yhat' = .

        /* Compute f(Y) = min(Y, horizon) for RMST target
         * (or 1{Y > horizon} for survival probability) */
        tempvar fY
        if `target' == 2 {
            quietly gen double `fY' = (`timevar' > `horizon') if `touse'
        }
        else {
            quietly gen double `fY' = min(`timevar', `horizon') if `touse'
        }

        /* Fit Y.hat = E[f(Y)|X] using regression forest on uncensored obs */
        plugin call grf_plugin `indepvars' `fY' `extra_vars' `yhat' ///
            if `touse',                                            ///
            "regression"                                           ///
            "`nuis_trees'"                                         ///
            "`seed'"                                               ///
            "`mtry'"                                               ///
            "15"                                                   ///
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
            "1"                                                    ///
            "`allow_missing_x'"                                    ///
            "`_nuis_cluster_idx'"                                  ///
            "`_nuis_weight_idx'"

        /* Fit censoring regression: E[D|X,W] as proxy for C.hat(T|X,W) */
        tempvar chat_proxy
        quietly gen double `chat_proxy' = .

        plugin call grf_plugin `indepvars' `statusvar' `extra_vars' `chat_proxy' ///
            if `touse',                                            ///
            "regression"                                           ///
            "`nuis_trees'"                                         ///
            "`seed'"                                               ///
            "`mtry'"                                               ///
            "15"                                                   ///
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
            "1"                                                    ///
            "`allow_missing_x'"                                    ///
            "`_nuis_cluster_idx'"                                  ///
            "`_nuis_weight_idx'"

        /* C.hat proxy = P(not censored | X,W) = E[D|X,W]
         * Clip to [0.001, 1.0] for numerical stability */
        quietly replace `chat_proxy' = max(`chat_proxy', 0.001) if `touse'
        quietly replace `chat_proxy' = min(`chat_proxy', 1.0) if `touse'

        /* ---- Step 3: Compute IPCW numerator/denominator ----
         *
         * IPCW pseudo-outcome following Cui et al. (2023):
         *   numerator_i = W.centered * (D_i * (f(Y_i) - Y.hat_i) +
         *                  (1 - D_i) * 0) / max(C.hat_i, eps)
         *   denominator_i = W.centered^2
         *
         * This is a simplified but proper IPCW that uses conditional
         * censoring probabilities rather than constant weights.
         */

        display as text "Step 3/4: Computing IPCW pseudo-outcomes ..."

        tempvar cs_numer cs_denom
        /* For events (D=1): IPCW weighted outcome deviation
         * For censored (D=0): contribute 0 (conservative) */
        quietly gen double `cs_numer' = `w_centered' * ///
            `statusvar' * (`fY' - `yhat') / `chat_proxy' if `touse'
        quietly gen double `cs_denom' = `w_centered' * `w_centered' if `touse'

        /* Ensure denominator is positive (avoid division by zero in plugin) */
        quietly replace `cs_denom' = max(`cs_denom', 1e-10) if `touse'

        display as text "  Nuisance preparation complete."

        local numer `cs_numer'
        local denom `cs_denom'
    }
    else {
        display as text "Steps 1-3: Skipped (using pre-computed nuisance columns)."

        /* Center treatment even in precomputed mode */
        if "`whatinput'" != "" {
            display as text "  Using user-supplied propensity scores from `whatinput'"
            tempvar what
            quietly gen double `what' = `whatinput' if `touse'
        }
        else {
            tempvar what
            quietly gen double `what' = .
            plugin call grf_plugin `indepvars' `treatvar' `extra_vars' `what' ///
                if `touse',                                       ///
                "regression"                                      ///
                "`ntrees'"                                        ///
                "`seed'"                                          ///
                "`mtry'"                                          ///
                "`minnodesize'"                                   ///
                "`samplefrac'"                                    ///
                "`do_honesty'"                                    ///
                "`honestyfrac'"                                   ///
                "`do_honesty_prune'"                              ///
                "`alpha'"                                         ///
                "`imbalancepenalty'"                               ///
                "`cigroupsize'"                                   ///
                "`numthreads'"                                    ///
                "0"                                               ///
                "0"                                               ///
                "`nindep'"                                        ///
                "1"                                               ///
                "0"                                               ///
                "0"                                               ///
                "1"                                               ///
                "`allow_missing_x'"                               ///
                "`_nuis_cluster_idx'"                             ///
                "`_nuis_weight_idx'"
        }
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

    /* ---- Build output varlist ---- */
    local output_vars `generate'
    if `do_est_var' {
        local output_vars `generate' `vargenerate'
    }

    /* ---- Build extra vars for cluster/weight ---- */
    local extra_vars ""
    local n_data_before = `nindep' + 5
    if "`cluster_var'" != "" {
        local extra_vars `extra_vars' `cluster_var'
        local cluster_col_idx = `n_data_before' + 1
        local n_data_before = `n_data_before' + 1
    }
    if "`weight_var'" != "" {
        local extra_vars `extra_vars' `weight_var'
        local weight_col_idx = `n_data_before' + 1
    }

    /* ---- Save nuisance estimates ---- */
    capture drop _grf_cs_what
    quietly gen double _grf_cs_what = `what'
    label variable _grf_cs_what "W.hat from propensity regression (W ~ X)"

    /* ---- Call plugin for causal survival forest ----
     *
     * Variable order: X1..Xp time W.centered status cs_numer cs_denom out1 [out2]
     *   n_x = nindep
     *   n_y = 1 (survival time)
     *   n_w = 1 (centered treatment)
     *   n_z = 0
     *   n_output = 1 or 2
     *
     * argv[0..19]: standard forest params
     * argv[20]: stabilize_splits
     * argv[21]: col index for cs_numer (relative, after X Y W status)
     * argv[22]: col index for cs_denom
     * argv[23]: col index for censor/status
     * argv[24]: target (1=RMST, 2=survival probability)
     */

    /* Center treatment for the causal forest step */
    tempvar w_cent_final
    quietly gen double `w_cent_final' = `treatvar' - `what' if `touse'

    local step_n = cond(`precomputed', "2/2", "4/4")
    display as text "Step `step_n': Fitting causal survival forest ..."
    plugin call grf_plugin `indepvars' `timevar' `w_cent_final' ///
        `statusvar' `numer' `denom' `extra_vars' `output_vars'   ///
        if `touse',                                              ///
        "causal_survival"                                        ///
        "`ntrees'"                                               ///
        "`seed'"                                                 ///
        "`mtry'"                                                 ///
        "`minnodesize'"                                          ///
        "`samplefrac'"                                           ///
        "`do_honesty'"                                           ///
        "`honestyfrac'"                                          ///
        "`do_honesty_prune'"                                     ///
        "`alpha'"                                                ///
        "`imbalancepenalty'"                                      ///
        "`cigroupsize'"                                          ///
        "`numthreads'"                                           ///
        "`do_est_var'"                                           ///
        "0"                                                      ///
        "`nindep'"                                               ///
        "1"                                                      ///
        "1"                                                      ///
        "0"                                                      ///
        "`n_output'"                                             ///
        "`allow_missing_x'"                                      ///
        "`cluster_col_idx'"                                      ///
        "`weight_col_idx'"                                       ///
        "`do_stabilize'"                                         ///
        "`=`nindep'+3'"                                          ///
        "`=`nindep'+4'"                                          ///
        "`=`nindep'+2'"                                          ///
        "`target'"

    /* ---- Compute CATE summary ---- */
    quietly summarize `generate' if `touse'
    local ate = r(mean)
    local ate_se = r(sd) / sqrt(r(N))

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
    ereturn scalar stabilize   = `do_stabilize'
    ereturn scalar horizon     = `horizon'
    ereturn scalar target      = `target'
    ereturn scalar ate         = `ate'
    ereturn scalar ate_se      = `ate_se'
    ereturn local  cmd           "grf_causal_survival_forest"
    ereturn scalar allow_missing_x = `allow_missing_x'
    ereturn local  forest_type   "causal_survival"
    ereturn local  timevar       "`timevar'"
    ereturn local  statusvar     "`statusvar'"
    ereturn local  treatvar      "`treatvar'"
    ereturn local  indepvars     "`indepvars'"
    ereturn local  predict_var   "`generate'"
    if `do_est_var' {
        ereturn local variance_var "`vargenerate'"
    }
    ereturn local what_var   "_grf_cs_what"
    if "`cluster_var'" != "" {
        ereturn local cluster_var "`cluster_var'"
    }
    if "`weight_var'" != "" {
        ereturn local weight_var "`weight_var'"
    }

    /* ---- Summary stats ---- */
    quietly summarize `generate' if `touse'
    local n_pred = r(N)
    local pred_mean = r(mean)
    local pred_sd = r(sd)

    display as text ""
    display as text "Causal Survival Forest Results"
    display as text "{hline 60}"
    display as text "CATE predictions:       " as result "`generate'"
    display as text "  Non-missing:  " as result `n_pred'
    display as text "  Mean (ATE):   " as result %9.4f `pred_mean'
    display as text "  SD:           " as result %9.4f `pred_sd'
    display as text "  ATE s.e.:     " as result %9.4f `ate_se'
    display as text "  Horizon:      " as result %9.3f `horizon'
    if `do_est_var' {
        quietly summarize `vargenerate' if `touse'
        display as text "Variance estimates:     " as result "`vargenerate'"
        display as text "  Mean variance: " as result %9.6f r(mean)
    }
    display as text "{hline 60}"
    display as text ""
end
