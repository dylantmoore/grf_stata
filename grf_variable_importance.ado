*! grf_variable_importance.ado -- Variable importance from grf
*! Version 0.1.0
*! Split-frequency weighted variable importance via regression forest

program define grf_variable_importance, rclass
    version 14.0

    syntax varlist(min=2 numeric) [if] [in],  ///
        [                                      ///
            NTrees(integer 2000)               ///
            SEED(integer 42)                   ///
            MAXDepth(integer 4)                ///
        ]

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
     * Variable order: X1..Xp Y out
     * argv: forest_type num_trees seed mtry min_node_size sample_fraction
     *       honesty honesty_fraction honesty_prune alpha imbalance_penalty
     *       ci_group_size num_threads estimate_variance compute_oob
     *       n_x n_y n_w n_z n_output
     *       max_depth (forest-specific)
     */
    plugin call grf_plugin `indepvars' `depvar' `outvar' ///
        if `touse',                                       ///
        "variable_importance"                              ///
        "`ntrees'"                                         ///
        "`seed'"                                           ///
        "0"                                                ///
        "5"                                                ///
        "0.5"                                              ///
        "1"                                                ///
        "0.5"                                              ///
        "1"                                                ///
        "0.05"                                             ///
        "0.0"                                              ///
        "1"                                                ///
        "0"                                                ///
        "0"                                                ///
        "0"                                                ///
        "`nindep'"                                         ///
        "1"                                                ///
        "0"                                                ///
        "0"                                                ///
        "1"                                                ///
        "`maxdepth'"

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
end
