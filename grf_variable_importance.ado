*! grf_variable_importance.ado -- Variable importance from grf
*! Version 0.2.0
*! Split-frequency weighted variable importance via regression forest

program define grf_variable_importance, rclass
    version 14.0

    syntax varlist(min=2 numeric) [if] [in],  ///
        [                                      ///
            NTrees(integer 2000)               ///
            SEED(integer 42)                   ///
            MAXDepth(integer 4)                ///
            DECAYexponent(real 2.0)            ///
            MTRY(integer 0)                    ///
            MINNodesize(integer 5)             ///
            SAMPLEfrac(real 0.5)               ///
            noHONesty                          ///
            HONestyfrac(real 0.5)              ///
            noHONestyprune                     ///
            ALPha(real 0.05)                   ///
            IMBalancepenalty(real 0.0)          ///
            NUMThreads(integer 0)              ///
            noMIA                              ///
            CLuster(varname numeric)           ///
            WEIghts(varname numeric)           ///
            EQUALizeclusterweights             ///
        ]

    /* ---- Parse varlist ---- */
    gettoken depvar indepvars : varlist
    local nindep : word count `indepvars'

    if `nindep' < 1 {
        display as error "need at least 1 predictor variable"
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

    /* ---- Mark sample ---- */
    if `allow_missing_x' {
        marksample touse, novarlist
        markout `touse' `depvar'
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
        tempvar _eq_clsize _eq_wt
        quietly bysort `cluster': gen long `_eq_clsize' = _N if `touse'
        quietly gen double `_eq_wt' = 1.0 / `_eq_clsize' if `touse'
        if "`weight_var'" != "" {
            quietly replace `_eq_wt' = `_eq_wt' * `weight_var' if `touse'
        }
        local weight_var `_eq_wt'
    }

    /* ---- Create dummy output variable (plugin requires it) ---- */
    tempvar outvar
    quietly gen double `outvar' = .

    /* ---- Display header ---- */
    display as text ""
    display as text "GRF Variable Importance"
    display as text "{hline 55}"
    display as text "Dependent variable:    " as result "`depvar'"
    display as text "Predictors:            " as result "`indepvars'"
    display as text "Observations:          " as result `n_use'
    display as text "Trees:                 " as result `ntrees'
    display as text "Max depth:             " as result `maxdepth'
    display as text "Decay exponent:        " as result `decayexponent'
    display as text "Honesty:               " as result cond(`do_honesty', "yes", "no")
    if "`cluster_var'" != "" {
        display as text "Cluster variable:      " as result "`cluster_var'"
    }
    if "`weight_var'" != "" {
        display as text "Weight variable:       " as result "`weight_var'"
    }
    display as text "{hline 55}"
    display as text ""

    /* ---- Load plugin ---- */
    if ( inlist("`c(os)'", "MacOSX") | strpos("`c(machine_type)'", "Mac") ) local c_os_ macosx
    else local c_os_: di lower("`c(os)'")

    capture program grf_plugin, plugin using("grf_plugin_`c_os_'.plugin")

    /* ---- Build extra vars for cluster/weight ---- */
    local extra_vars ""
    local n_data_before = `nindep' + 1
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
     * Variable order: X1..Xp Y [cluster] [weight] out
     * argv: forest_type num_trees seed mtry min_node_size sample_fraction
     *       honesty honesty_fraction honesty_prune alpha imbalance_penalty
     *       ci_group_size num_threads estimate_variance compute_oob
     *       n_x n_y n_w n_z n_output
     *       allow_missing_x cluster_col_idx weight_col_idx
     *       max_depth decay_exponent
     */
    plugin call grf_plugin `indepvars' `depvar' `extra_vars' `outvar' ///
        if `touse',                                        ///
        "variable_importance"                              ///
        "`ntrees'"                                         ///
        "`seed'"                                           ///
        "`mtry'"                                           ///
        "`minnodesize'"                                    ///
        "`samplefrac'"                                     ///
        "`do_honesty'"                                     ///
        "`honestyfrac'"                                    ///
        "`do_honesty_prune'"                               ///
        "`alpha'"                                          ///
        "`imbalancepenalty'"                                ///
        "1"                                                ///
        "`numthreads'"                                     ///
        "0"                                                ///
        "0"                                                ///
        "`nindep'"                                         ///
        "1"                                                ///
        "0"                                                ///
        "0"                                                ///
        "1"                                                ///
        "`allow_missing_x'"                                ///
        "`cluster_col_idx'"                                ///
        "`weight_col_idx'"                                 ///
        "`maxdepth'"                                       ///
        "`decayexponent'"

    /* ---- Read scalars and build results ---- */
    local vi_n = scalar(_grf_vi_n)

    if missing(`vi_n') | `vi_n' < 1 {
        display as error "plugin did not return variable importance scores"
        exit 498
    }

    /* Build importance matrix */
    tempname imp_mat
    matrix `imp_mat' = J(1, `nindep', 0)

    local colnames ""
    forvalues v = 1/`nindep' {
        local vi_val = scalar(_grf_vi_`v')
        matrix `imp_mat'[1, `v'] = `vi_val'
        local vname : word `v' of `indepvars'
        local colnames "`colnames' `vname'"
    }
    matrix colnames `imp_mat' = `colnames'
    matrix rownames `imp_mat' = "importance"

    /* ---- Display results table ---- */
    display as text ""
    display as text %~20s "Variable" %~15s "Importance"
    display as text "{hline 38}"

    forvalues v = 1/`nindep' {
        local vname : word `v' of `indepvars'
        local vi_val = `imp_mat'[1, `v']
        display as text %20s "`vname'" _col(25) as result %10.6f `vi_val'
    }

    display as text "{hline 38}"
    display as text ""

    /* ---- Clean up scalars ---- */
    forvalues v = 1/`nindep' {
        scalar drop _grf_vi_`v'
    }
    scalar drop _grf_vi_n

    /* ---- Store results ---- */
    return matrix importance = `imp_mat'
    return scalar N = `n_use'
    return scalar n_trees = `ntrees'
    return scalar seed = `seed'
    return scalar max_depth = `maxdepth'
    return scalar decay_exponent = `decayexponent'
    return scalar mtry = `mtry'
    return scalar min_node_size = `minnodesize'
    return scalar sample_fraction = `samplefrac'
    return scalar honesty = `do_honesty'
    return scalar alpha = `alpha'
end
