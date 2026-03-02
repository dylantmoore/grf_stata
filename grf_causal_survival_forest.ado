*! grf_causal_survival_forest.ado -- Causal Survival Forest via grf C++ library
*! Version 0.3.0
*! Implements Cui et al. (2023) causal_survival_forest()

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
            IMBalancepenalty(real 0.0)         ///
            CIGroupsize(integer 1)             ///
            NUMThreads(integer 0)              ///
            ESTIMATEVariance                   ///
            noSTABilizesplits                  ///
            NUMer(varname numeric)             ///
            DENom(varname numeric)             ///
            WHATinput(varname numeric)         ///
            YHATinput(varname numeric)         ///
            SHATinput(varname numeric)         ///
            CHATinput(varname numeric)         ///
            WHATGenerate(name)                 ///
            YHATGenerate(name)                 ///
            SHATGenerate(name)                 ///
            CHATGenerate(name)                 ///
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

    /* ---- Validate nuisance-input combinations ---- */
    if ("`numer'" != "" & "`denom'" == "") | ("`numer'" == "" & "`denom'" != "") {
        display as error "numer() and denom() must be supplied together"
        exit 198
    }

    local n_partial_inputs 0
    foreach _opt in whatinput yhatinput shatinput chatinput {
        if "``_opt''" != "" {
            local n_partial_inputs = `n_partial_inputs' + 1
        }
    }

    if "`numer'" != "" & `n_partial_inputs' > 0 {
        display as error "numer()/denom() is mutually exclusive with whatinput()/yhatinput()/shatinput()/chatinput()"
        exit 198
    }

    if "`numer'" == "" & `n_partial_inputs' > 0 & `n_partial_inputs' < 4 {
        display as error "full nuisance input mode requires all of whatinput(), yhatinput(), shatinput(), and chatinput()"
        exit 198
    }

    if !inlist(`target', 1, 2) {
        display as error "target() must be 1 (RMST) or 2 (survival probability)"
        exit 198
    }

    local nuisance_mode "auto"
    if "`numer'" != "" {
        local nuisance_mode "moment_input"
    }
    else if `n_partial_inputs' == 4 {
        local nuisance_mode "full_input"
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
        foreach _g in whatgenerate yhatgenerate shatgenerate chatgenerate {
            if "``_g''" != "" {
                capture drop ``_g''
            }
        }
    }
    confirm new variable `generate'
    foreach _g in whatgenerate yhatgenerate shatgenerate chatgenerate {
        if "``_g''" != "" {
            confirm new variable ``_g''
        }
    }

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

    if "`numer'" != "" {
        markout `touse' `numer' `denom'
    }
    if "`whatinput'" != "" {
        markout `touse' `whatinput'
    }
    if "`yhatinput'" != "" {
        markout `touse' `yhatinput'
    }
    if "`shatinput'" != "" {
        markout `touse' `shatinput'
    }
    if "`chatinput'" != "" {
        markout `touse' `chatinput'
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
        tempvar _eq_clsize _eq_wt
        quietly bysort `cluster': gen long `_eq_clsize' = _N if `touse'
        quietly gen double `_eq_wt' = 1.0 / `_eq_clsize' if `touse'
        if "`weight_var'" != "" {
            quietly replace `_eq_wt' = `_eq_wt' * `weight_var' if `touse'
        }
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
    quietly count if !inlist(`treatvar', 0, 1) & `touse'
    local binary_treat = (r(N) == 0)

    /* ---- Determine horizon ---- */
    if `horizon' == 0 {
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

        grf_tune `varlist' if `touse', foresttype(causal_survival) ///
            numreps(`tunenumreps') tunetrees(`tunenumtrees') seed(`seed') ///
            numthreads(`numthreads') horizon(`horizon')

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
    display as text "Nuisance mode:         " as result "`nuisance_mode'"
    if `do_est_var' {
        display as text "Variance estimation:   " as result "yes (ci_group_size=`cigroupsize')"
    }
    display as text "{hline 60}"
    display as text ""

    /* ---- Load plugin ---- */
    if ( inlist("`c(os)'", "MacOSX") | strpos("`c(machine_type)'", "Mac") ) local c_os_ macosx
    else local c_os_: di lower("`c(os)'")

    capture program grf_plugin, plugin using("grf_plugin_`c_os_'.plugin")

    /* Nuisance regression indices: X + Y (+ optional cluster/weight) */
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

    local extra_vars ""
    if "`cluster_var'" != "" {
        local extra_vars `extra_vars' `cluster_var'
    }
    if "`weight_var'" != "" {
        local extra_vars `extra_vars' `weight_var'
    }

    local what_source ""
    local yhat_source ""
    local shat_source ""
    local chat_source ""
    local numer_source "`numer'"
    local denom_source "`denom'"

    /* Shared transformed outcome f(Y) */
    tempvar fY
    if `target' == 2 {
        quietly gen double `fY' = (`timevar' > `horizon') if `touse'
    }
    else {
        quietly gen double `fY' = min(`timevar', `horizon') if `touse'
    }

    if "`nuisance_mode'" == "full_input" {
        display as text "Step 1/4: Using full user-supplied nuisance inputs"

        tempvar what yhat shat chat_proxy w_centered cs_numer cs_denom
        quietly gen double `what' = `whatinput' if `touse'
        quietly gen double `yhat' = `yhatinput' if `touse'
        quietly gen double `shat' = `shatinput' if `touse'
        quietly gen double `chat_proxy' = `chatinput' if `touse'

        if `binary_treat' {
            quietly replace `what' = min(max(`what', 1e-6), 1 - 1e-6) if `touse'
        }
        quietly replace `shat' = min(max(`shat', 0.001), 1.0) if `touse'
        quietly replace `chat_proxy' = min(max(`chat_proxy', 0.001), 1.0) if `touse'

        quietly gen double `w_centered' = `treatvar' - `what' if `touse'
        quietly gen double `cs_numer' = .
        if `target' == 2 {
            quietly replace `cs_numer' = `w_centered' * ///
                (`statusvar' * (`fY' - `yhat') + (1 - `statusvar') * (`shat' - `yhat')) / `chat_proxy' if `touse'
        }
        else {
            quietly replace `cs_numer' = `w_centered' * ///
                `statusvar' * (`fY' - `yhat') / `chat_proxy' if `touse'
        }
        quietly gen double `cs_denom' = max(`w_centered' * `w_centered', 1e-10) if `touse'

        local what_source `what'
        local yhat_source `yhat'
        local shat_source `shat'
        local chat_source `chat_proxy'
        local numer_source `cs_numer'
        local denom_source `cs_denom'

        display as text "  Full-input nuisance moments prepared."
    }
    else if "`nuisance_mode'" == "auto" {
        local nuis_trees = max(50, min(`ntrees' / 4, 500))

        display as text "Step 1/4: Estimating propensity scores W.hat = E[W|X]"
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
            "`imbalancepenalty'"                              ///
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

        if `binary_treat' {
            quietly replace `what' = min(max(`what', 1e-6), 1 - 1e-6) if `touse'
        }
        local what_source `what'

        display as text "Step 2/4: Estimating Y.hat, S.hat, and C.hat"

        tempvar yhat
        quietly gen double `yhat' = .
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
            "`imbalancepenalty'"                                   ///
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
        local yhat_source `yhat'

        tempvar surv_ind shat
        quietly gen double `surv_ind' = (`timevar' > `horizon') if `touse'
        quietly gen double `shat' = .
        plugin call grf_plugin `indepvars' `surv_ind' `extra_vars' `shat' ///
            if `touse',                                               ///
            "regression"                                              ///
            "`nuis_trees'"                                            ///
            "`seed'"                                                  ///
            "`mtry'"                                                  ///
            "15"                                                      ///
            "`samplefrac'"                                            ///
            "`do_honesty'"                                            ///
            "`honestyfrac'"                                           ///
            "`do_honesty_prune'"                                      ///
            "`alpha'"                                                 ///
            "`imbalancepenalty'"                                      ///
            "1"                                                       ///
            "`numthreads'"                                            ///
            "0"                                                       ///
            "0"                                                       ///
            "`nindep'"                                                ///
            "1"                                                       ///
            "0"                                                       ///
            "0"                                                       ///
            "1"                                                       ///
            "`allow_missing_x'"                                       ///
            "`_nuis_cluster_idx'"                                     ///
            "`_nuis_weight_idx'"
        quietly replace `shat' = min(max(`shat', 0.001), 1.0) if `touse'
        local shat_source `shat'

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
            "`imbalancepenalty'"                                   ///
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
        quietly replace `chat_proxy' = min(max(`chat_proxy', 0.001), 1.0) if `touse'
        local chat_source `chat_proxy'

        display as text "Step 3/4: Computing nuisance moments"
        tempvar w_centered cs_numer cs_denom
        quietly gen double `w_centered' = `treatvar' - `what' if `touse'
        quietly gen double `cs_numer' = .
        if `target' == 2 {
            quietly replace `cs_numer' = `w_centered' * ///
                (`statusvar' * (`fY' - `yhat') + (1 - `statusvar') * (`shat' - `yhat')) / `chat_proxy' if `touse'
        }
        else {
            quietly replace `cs_numer' = `w_centered' * ///
                `statusvar' * (`fY' - `yhat') / `chat_proxy' if `touse'
        }
        quietly gen double `cs_denom' = max(`w_centered' * `w_centered', 1e-10) if `touse'
        local numer_source `cs_numer'
        local denom_source `cs_denom'
    }
    else {
        display as text "Step 1/2: Using moment input override numer()/denom()"
        display as text "  Numerator:   `numer'"
        display as text "  Denominator: `denom'"

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
            "`imbalancepenalty'"                              ///
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
        if `binary_treat' {
            quietly replace `what' = min(max(`what', 1e-6), 1 - 1e-6) if `touse'
        }
        local what_source `what'
    }

    /* ---- Persist canonical nuisance artifacts ---- */
    capture drop _grf_cs_what
    quietly gen double _grf_cs_what = .
    quietly replace _grf_cs_what = `what_source' if `touse'
    label variable _grf_cs_what "W.hat nuisance used by causal-survival wrapper"

    capture drop _grf_cs_yhat
    quietly gen double _grf_cs_yhat = .
    if "`yhat_source'" != "" {
        quietly replace _grf_cs_yhat = `yhat_source' if `touse'
    }
    label variable _grf_cs_yhat "Y.hat nuisance used by causal-survival wrapper"

    capture drop _grf_cs_shat
    quietly gen double _grf_cs_shat = .
    if "`shat_source'" != "" {
        quietly replace _grf_cs_shat = `shat_source' if `touse'
    }
    label variable _grf_cs_shat "S.hat nuisance used by causal-survival wrapper"

    capture drop _grf_cs_chat
    quietly gen double _grf_cs_chat = .
    if "`chat_source'" != "" {
        quietly replace _grf_cs_chat = `chat_source' if `touse'
    }
    label variable _grf_cs_chat "C.hat nuisance used by causal-survival wrapper"

    capture drop _grf_cs_numer
    quietly gen double _grf_cs_numer = .
    quietly replace _grf_cs_numer = `numer_source' if `touse'
    label variable _grf_cs_numer "Causal-survival nuisance numerator"

    capture drop _grf_cs_denom
    quietly gen double _grf_cs_denom = .
    quietly replace _grf_cs_denom = `denom_source' if `touse'
    label variable _grf_cs_denom "Causal-survival nuisance denominator"

    /* ---- Optional user-facing nuisance outputs ---- */
    if "`whatgenerate'" != "" {
        quietly gen double `whatgenerate' = _grf_cs_what if `touse'
        label variable `whatgenerate' "Stored W.hat nuisance from causal-survival fit"
    }
    if "`yhatgenerate'" != "" {
        quietly gen double `yhatgenerate' = _grf_cs_yhat if `touse'
        label variable `yhatgenerate' "Stored Y.hat nuisance from causal-survival fit"
    }
    if "`shatgenerate'" != "" {
        quietly gen double `shatgenerate' = _grf_cs_shat if `touse'
        label variable `shatgenerate' "Stored S.hat nuisance from causal-survival fit"
    }
    if "`chatgenerate'" != "" {
        quietly gen double `chatgenerate' = _grf_cs_chat if `touse'
        label variable `chatgenerate' "Stored C.hat nuisance from causal-survival fit"
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

    /* ---- Call plugin for causal survival forest ---- */
    tempvar w_cent_final
    quietly gen double `w_cent_final' = `treatvar' - _grf_cs_what if `touse'

    local step_n = cond("`nuisance_mode'" == "moment_input", "2/2", "4/4")
    display as text "Step `step_n': Fitting causal survival forest ..."
    plugin call grf_plugin `indepvars' `timevar' `w_cent_final' ///
        `statusvar' _grf_cs_numer _grf_cs_denom `extra_vars' `output_vars' ///
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
        "`imbalancepenalty'"                                     ///
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
    capture scalar __grf_model_counter = __grf_model_counter + 1
    if _rc {
        scalar __grf_model_counter = 1
    }
    ereturn scalar model_id = __grf_model_counter
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
    ereturn scalar allow_missing_x = `allow_missing_x'
    ereturn local cmd          "grf_causal_survival_forest"
    ereturn local forest_type  "causal_survival"
    ereturn local timevar      "`timevar'"
    ereturn local statusvar    "`statusvar'"
    ereturn local treatvar     "`treatvar'"
    ereturn local indepvars    "`indepvars'"
    ereturn local predict_var  "`generate'"
    if `do_est_var' {
        ereturn local variance_var "`vargenerate'"
    }

    ereturn local what_var  "_grf_cs_what"
    ereturn local yhat_var  "_grf_cs_yhat"
    ereturn local shat_var  "_grf_cs_shat"
    ereturn local chat_var  "_grf_cs_chat"
    ereturn local numer_var "_grf_cs_numer"
    ereturn local denom_var "_grf_cs_denom"

    if "`whatgenerate'" != "" {
        ereturn local what_generate "`whatgenerate'"
    }
    if "`yhatgenerate'" != "" {
        ereturn local yhat_generate "`yhatgenerate'"
    }
    if "`shatgenerate'" != "" {
        ereturn local shat_generate "`shatgenerate'"
    }
    if "`chatgenerate'" != "" {
        ereturn local chat_generate "`chatgenerate'"
    }

    ereturn local nuisance_mode "`nuisance_mode'"

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
