*! grf_predict.ado -- Predict on new data using previously estimated GRF
*! Version 0.1.0
*!
*! Usage:
*!   1. Estimate a forest: grf_regression_forest y x1 x2 x3, gen(yhat)
*!   2. Append test data:  append using test_data
*!   3. Predict:           grf_predict, gen(yhat_new)
*!
*! The command reads e() results to determine:
*!   - forest_type, indepvars, depvar, treatvar, etc.
*!   - n_trees, seed, mtry, min_node, alpha, honesty, etc.
*!   - e(N) as the number of training observations
*!
*! Training obs are the first e(N) observations in the current dataset.
*! Test obs are the remaining observations (must have non-missing X values).
*!
*! Strategy: Instead of modifying the user's data columns, we create tempvar
*! copies of Y/W/Z with missing values filled in for test obs.  The plugin
*! reads from tempvar copies, writes to the output variable, and tempvars
*! are automatically cleaned up when the program exits.

program define grf_predict, rclass
    version 14.0

    syntax , GENerate(name) [REPlace NUMThreads(integer 0)]

    /* ================================================================
     * Step 1: Read e() results from prior estimation
     * ================================================================ */
    local forest_type "`e(forest_type)'"
    if "`forest_type'" == "" {
        display as error "no prior GRF estimation found"
        display as error "run grf_regression_forest, grf_causal_forest, etc. first"
        exit 301
    }

    local est_cmd     "`e(cmd)'"
    local n_train     = e(N)
    local n_trees     = e(n_trees)
    local seed        = e(seed)
    local mtry        = e(mtry)
    local min_node    = e(min_node)
    local alpha       = e(alpha)
    local do_honesty  = e(honesty)
    local indepvars   "`e(indepvars)'"

    /* Read tuning parameters from e() -- use estimation-time values,
     * not arbitrary defaults, so the predict forest matches the original */
    local samplefrac        = e(sample_fraction)
    local honestyfrac       = e(honesty_fraction)
    local imbalancepenalty  = e(imbalance_penalty)
    local cigroupsize       = e(ci_group_size)

    /* Fallback defaults if e() doesn't store these (older estimations) */
    if missing(`samplefrac')       local samplefrac       = 0.5
    if missing(`honestyfrac')      local honestyfrac      = 0.5
    if missing(`imbalancepenalty') local imbalancepenalty  = 0.0
    if missing(`cigroupsize')      local cigroupsize      = 1

    /* Read allow_missing_x from e() -- inherit from estimation */
    local allow_missing_x = e(allow_missing_x)
    if missing(`allow_missing_x') {
        local allow_missing_x 1
    }

    /* Forest-type-specific e() results */
    local depvar      "`e(depvar)'"

    if "`forest_type'" == "causal" {
        local treatvar     "`e(treatvar)'"
        local do_stabilize = e(stabilize)
    }

    if "`forest_type'" == "instrumental" {
        local treatvar        "`e(treatment)'"
        local instrvar        "`e(instrument)'"
        local do_stabilize    = e(stabilize_splits)
        local reduced_form_wt = e(reduced_form_wt)
        if missing(`reduced_form_wt') local reduced_form_wt 0
    }

    if "`forest_type'" == "survival" {
        local timevar     "`e(timevar)'"
        local statusvar   "`e(statusvar)'"
        local n_output_sv = e(n_output)
        local pred_type   = e(pred_type)
        if missing(`pred_type') local pred_type 0
    }

    if "`forest_type'" == "quantile" {
        local quantiles   "`e(quantiles)'"
        local n_quantiles = e(n_quantiles)
    }

    if "`forest_type'" == "probability" {
        local n_classes = e(n_classes)
    }

    if "`forest_type'" == "causal_survival" {
        local treatvar      "`e(treatvar)'"
        local timevar       "`e(timevar)'"
        local statusvar     "`e(statusvar)'"
        local cs_horizon    = e(horizon)
        local cs_target     = e(target)
        local do_stabilize  = e(stabilize)
    }

    if "`forest_type'" == "multi_causal" {
        local depvar        "`e(depvar)'"
        local n_treat       = e(n_treat)
        local treatvars     "`e(treatvars)'"
        local do_stabilize  = e(stabilize)
    }

    if "`forest_type'" == "multi_regression" {
        local n_outcomes    = e(n_outcomes)
        local depvars       "`e(depvars)'"
    }

    if "`forest_type'" == "ll_regression" {
        local enable_ll_split    = e(enable_ll_split)
        local ll_lambda          = e(ll_lambda)
        local ll_weight_penalty  = e(ll_weight_penalty)
        local ll_split_cutoff    = e(ll_split_cutoff)
    }

    if "`forest_type'" == "lm_forest" {
        local regvars        "`e(regvars)'"
        local n_regressors   = e(n_regressors)
        local do_stabilize   = e(stabilize)
    }

    /* boosted_regression: no extra e() results needed */

    /* ================================================================
     * Step 2: Validate current dataset
     * ================================================================ */

    /* Validate that we have more rows than the training set */
    quietly count
    local n_total = r(N)

    if `n_total' <= `n_train' {
        display as error "current dataset has `n_total' observations,"
        display as error "but the training set had `n_train'"
        display as error "append test data before calling grf_predict"
        exit 198
    }

    local n_test = `n_total' - `n_train'

    /* Validate that X variables are non-missing in test obs */
    foreach v of local indepvars {
        quietly count if missing(`v') & _n > `n_train'
        if r(N) > 0 {
            display as error "predictor `v' has missing values in test observations"
            exit 416
        }
    }

    /* ================================================================
     * Step 3: Handle replace and create output variable(s)
     * ================================================================ */

    if "`forest_type'" == "regression"      | ///
       "`forest_type'" == "causal"           | ///
       "`forest_type'" == "instrumental"     | ///
       "`forest_type'" == "causal_survival"  | ///
       "`forest_type'" == "ll_regression" {

        /* Single output variable */
        if "`replace'" != "" {
            capture drop `generate'
        }
        confirm new variable `generate'
        quietly gen double `generate' = .
    }
    else if "`forest_type'" == "quantile" {
        /* Multiple output: stub_qNN for each quantile */
        foreach q of local quantiles {
            local qint = round(`q' * 100)
            local varname `generate'_q`qint'
            if "`replace'" != "" {
                capture drop `varname'
            }
            confirm new variable `varname'
            quietly gen double `varname' = .
        }
    }
    else if "`forest_type'" == "probability" {
        /* Multiple output: stub_cN for each class */
        forvalues c = 0/`=`n_classes'-1' {
            local varname `generate'_c`c'
            if "`replace'" != "" {
                capture drop `varname'
            }
            confirm new variable `varname'
            quietly gen double `varname' = .
        }
    }
    else if "`forest_type'" == "survival" {
        /* Multiple output: stub_sN for each time point */
        forvalues j = 1/`n_output_sv' {
            local varname `generate'_s`j'
            if "`replace'" != "" {
                capture drop `varname'
            }
            confirm new variable `varname'
            quietly gen double `varname' = .
        }
    }
    else if "`forest_type'" == "multi_causal" {
        /* Per-arm output: stub_t1, stub_t2, ... */
        forvalues j = 1/`n_treat' {
            local varname `generate'_t`j'
            if "`replace'" != "" {
                capture drop `varname'
            }
            confirm new variable `varname'
            quietly gen double `varname' = .
        }
    }
    else if "`forest_type'" == "multi_regression" {
        /* Per-outcome output: stub_y1, stub_y2, ... */
        forvalues j = 1/`n_outcomes' {
            local varname `generate'_y`j'
            if "`replace'" != "" {
                capture drop `varname'
            }
            confirm new variable `varname'
            quietly gen double `varname' = .
        }
    }
    else if "`forest_type'" == "lm_forest" {
        /* Per-coefficient output: stub_1, stub_2, ... */
        forvalues j = 1/`n_regressors' {
            local varname `generate'_`j'
            if "`replace'" != "" {
                capture drop `varname'
            }
            confirm new variable `varname'
            quietly gen double `varname' = .
        }
    }
    else if "`forest_type'" == "boosted_regression" {
        /* Single output variable -- but predict not supported */
        display as error "grf_predict does not support boosted_regression forests"
        display as error "the C++ plugin only supports OOB predictions for boosted regression;"
        display as error "predict mode (on new data) is not available for this forest type"
        exit 198
    }
    else {
        display as error "forest type `forest_type' is not yet supported for predict"
        exit 198
    }

    /* ================================================================
     * Step 4: Display header
     * ================================================================ */
    display as text ""
    display as text "GRF Predict: `forest_type' forest"
    display as text "{hline 55}"
    display as text "Prior estimation:      " as result "`est_cmd'"
    display as text "Training observations: " as result `n_train'
    display as text "Test observations:     " as result `n_test'
    display as text "Predictors:            " as result "`indepvars'"
    display as text "Trees:                 " as result `n_trees'
    display as text "{hline 55}"
    display as text ""

    /* ================================================================
     * Step 5: Load plugin
     * ================================================================ */
    if ( inlist("`c(os)'", "MacOSX") | strpos("`c(machine_type)'", "Mac") ) local c_os_ macosx
    else local c_os_: di lower("`c(os)'")

    capture program grf_plugin, plugin using("grf_plugin_`c_os_'.plugin")

    /* ================================================================
     * Step 6: Dispatch by forest type
     * ================================================================
     *
     * For predict mode, we pass n_train as argv[14].  The C++ plugin
     * uses the first n_train rows as training data and the remaining
     * rows as test data.
     *
     * To avoid modifying the user's data, we create tempvar copies of
     * Y/W/Z columns with missing values in test rows filled with 0.
     * The plugin reads from these tempvar copies.  Tempvars are
     * automatically dropped when the program exits.
     *
     * argv layout (common 23 args):
     *   [0]  forest_type
     *   [1]  num_trees
     *   [2]  seed
     *   [3]  mtry
     *   [4]  min_node_size
     *   [5]  sample_fraction
     *   [6]  honesty
     *   [7]  honesty_fraction
     *   [8]  honesty_prune_leaves
     *   [9]  alpha
     *   [10] imbalance_penalty
     *   [11] ci_group_size
     *   [12] num_threads
     *   [13] estimate_variance (0 for predict)
     *   [14] n_train (>0 triggers predict mode)
     *   [15] n_x
     *   [16] n_y
     *   [17] n_w
     *   [18] n_z
     *   [19] n_output
     *   [20] allow_missing_x (1=MIA enabled)
     *   [21] cluster_col_idx (0=no clustering)
     *   [22] weight_col_idx (0=no weights)
     */

    local nindep : word count `indepvars'

    /* ---- honesty_prune: read from e(), default to 1 ---- */
    local do_honesty_prune = e(honesty_prune)
    if missing(`do_honesty_prune') local do_honesty_prune 1

    /* ----------------------------------------------------------------
     * REGRESSION FOREST predict
     * ----------------------------------------------------------------
     * Pass: X1..Xp Y_copy output_var
     * Y_copy is a tempvar with missing values in test rows set to 0.
     * The plugin trains on the first n_train rows and predicts on the rest.
     * ---------------------------------------------------------------- */
    if "`forest_type'" == "regression" {

        display as text "Predicting with regression forest ..."

        /* Create tempvar copy of depvar with missing filled for test obs */
        tempvar y_safe
        quietly gen double `y_safe' = `depvar'
        quietly replace `y_safe' = 0 if _n > `n_train' & missing(`y_safe')

        /* Call plugin: X1..Xp Y_safe output_var */
        plugin call grf_plugin `indepvars' `y_safe' `generate', ///
            "regression"                                         ///
            "`n_trees'"                                          ///
            "`seed'"                                             ///
            "`mtry'"                                             ///
            "`min_node'"                                         ///
            "`samplefrac'"                                       ///
            "`do_honesty'"                                       ///
            "`honestyfrac'"                                      ///
            "`do_honesty_prune'"                                 ///
            "`alpha'"                                            ///
            "`imbalancepenalty'"                                  ///
            "`cigroupsize'"                                      ///
            "`numthreads'"                                       ///
            "0"                                                  ///
            "`n_train'"                                          ///
            "`nindep'"                                           ///
            "1"                                                  ///
            "0"                                                  ///
            "0"                                                  ///
            "1"                                                  ///
            "`allow_missing_x'"                                  ///
            "0"                                                  ///
            "0"

        /* Clear predictions for training obs (they got OOB predictions,
         * but the user only asked for test predictions) */
        quietly replace `generate' = . if _n <= `n_train'
    }

    /* ----------------------------------------------------------------
     * CAUSAL FOREST predict
     * ----------------------------------------------------------------
     *
     * The causal forest requires nuisance estimation:
     *   Step 1: Regression forest Y ~ X (OOB on training) => Y.hat for train obs
     *   Step 2: Regression forest W ~ X (OOB on training) => W.hat for train obs
     *   Step 3: Center training: Y.c = Y - Y.hat, W.c = W - W.hat
     *   Step 4: Causal forest on (X, Y.c, W.c) with n_train => CATE for test obs
     *
     * The nuisance models run in OOB mode on TRAINING data only.
     * We only need Y.hat/W.hat for centering training obs.
     * For test obs, Y.c/W.c are 0 placeholders (the plugin only
     * uses X columns from test data for prediction).
     * ---------------------------------------------------------------- */
    else if "`forest_type'" == "causal" {

        /* --- Step 1: Nuisance model Y ~ X (OOB on training only) --- */
        display as text "Step 1/3: Nuisance model Y ~ X (OOB on training) ..."

        tempvar yhat_pred
        quietly gen double `yhat_pred' = .

        plugin call grf_plugin `indepvars' `depvar' `yhat_pred' ///
            if _n <= `n_train',                                  ///
            "regression"                                          ///
            "`n_trees'"                                             ///
            "`seed'"                                              ///
            "`mtry'"                                              ///
            "`min_node'"                                          ///
            "`samplefrac'"                                        ///
            "`do_honesty'"                                        ///
            "`honestyfrac'"                                       ///
            "`do_honesty_prune'"                                  ///
            "`alpha'"                                             ///
            "`imbalancepenalty'"                                   ///
            "`cigroupsize'"                                       ///
            "`numthreads'"                                        ///
            "0"                                                   ///
            "0"                                                   ///
            "`nindep'"                                            ///
            "1"                                                   ///
            "0"                                                   ///
            "0"                                                   ///
            "1"                                                   ///
            "`allow_missing_x'"                                   ///
            "0"                                                   ///
            "0"

        /* --- Step 2: Nuisance model W ~ X (OOB on training only) --- */
        display as text "Step 2/3: Nuisance model W ~ X (OOB on training) ..."

        tempvar what_pred
        quietly gen double `what_pred' = .

        plugin call grf_plugin `indepvars' `treatvar' `what_pred' ///
            if _n <= `n_train',                                    ///
            "regression"                                           ///
            "`n_trees'"                                            ///
            "`seed'"                                               ///
            "`mtry'"                                               ///
            "`min_node'"                                           ///
            "`samplefrac'"                                         ///
            "`do_honesty'"                                         ///
            "`honestyfrac'"                                        ///
            "`do_honesty_prune'"                                   ///
            "`alpha'"                                              ///
            "`imbalancepenalty'"                                    ///
            "`cigroupsize'"                                        ///
            "`numthreads'"                                         ///
            "0"                                                    ///
            "0"                                                    ///
            "`nindep'"                                             ///
            "1"                                                    ///
            "0"                                                    ///
            "0"                                                    ///
            "1"                                                    ///
            "`allow_missing_x'"                                    ///
            "0"                                                    ///
            "0"

        /* --- Step 3: Center Y and W, then run causal forest --- */
        display as text "Step 3/3: Causal forest on centered data (with predict) ..."

        tempvar y_centered w_centered

        /* Center training obs using OOB Y.hat and W.hat */
        quietly gen double `y_centered' = `depvar' - `yhat_pred' if _n <= `n_train'
        quietly gen double `w_centered' = `treatvar' - `what_pred' if _n <= `n_train'

        /* For test obs: fill with 0 placeholder.
         * The plugin only reads Y.c/W.c from training rows for fitting;
         * test rows only need X columns for prediction. */
        quietly replace `y_centered' = 0 if _n > `n_train'
        quietly replace `w_centered' = 0 if _n > `n_train'

        /* Call causal forest with n_train */
        plugin call grf_plugin `indepvars' `y_centered' `w_centered' `generate', ///
            "causal"                                                              ///
            "`n_trees'"                                                           ///
            "`seed'"                                                              ///
            "`mtry'"                                                              ///
            "`min_node'"                                                          ///
            "`samplefrac'"                                                        ///
            "`do_honesty'"                                                        ///
            "`honestyfrac'"                                                       ///
            "`do_honesty_prune'"                                                  ///
            "`alpha'"                                                             ///
            "`imbalancepenalty'"                                                   ///
            "`cigroupsize'"                                                       ///
            "`numthreads'"                                                        ///
            "0"                                                                   ///
            "`n_train'"                                                           ///
            "`nindep'"                                                            ///
            "1"                                                                   ///
            "1"                                                                   ///
            "0"                                                                   ///
            "1"                                                                   ///
            "`allow_missing_x'"                                                   ///
            "0"                                                                   ///
            "0"                                                                   ///
            "`do_stabilize'"

        /* Clear predictions for training obs */
        quietly replace `generate' = . if _n <= `n_train'
    }

    /* ----------------------------------------------------------------
     * QUANTILE FOREST predict
     * ----------------------------------------------------------------
     * Pass: X1..Xp Y_copy out_q10 out_q25 ...
     * The quantile CSV tells the plugin which quantiles to estimate.
     * ---------------------------------------------------------------- */
    else if "`forest_type'" == "quantile" {

        display as text "Predicting with quantile forest ..."

        /* Build comma-separated quantile list for the plugin */
        local quantile_csv ""
        local first 1
        foreach q of local quantiles {
            if `first' {
                local quantile_csv "`q'"
                local first 0
            }
            else {
                local quantile_csv "`quantile_csv',`q'"
            }
        }

        /* Build output varlist */
        local output_vars ""
        foreach q of local quantiles {
            local qint = round(`q' * 100)
            local output_vars `output_vars' `generate'_q`qint'
        }

        /* Create tempvar copy of depvar with missing filled for test obs */
        tempvar y_safe
        quietly gen double `y_safe' = `depvar'
        quietly replace `y_safe' = 0 if _n > `n_train' & missing(`y_safe')

        /* Call plugin */
        plugin call grf_plugin `indepvars' `y_safe' `output_vars', ///
            "quantile"                                              ///
            "`n_trees'"                                             ///
            "`seed'"                                                ///
            "`mtry'"                                                ///
            "`min_node'"                                            ///
            "`samplefrac'"                                          ///
            "`do_honesty'"                                          ///
            "`honestyfrac'"                                         ///
            "`do_honesty_prune'"                                    ///
            "`alpha'"                                               ///
            "`imbalancepenalty'"                                     ///
            "`cigroupsize'"                                         ///
            "`numthreads'"                                          ///
            "0"                                                     ///
            "`n_train'"                                             ///
            "`nindep'"                                              ///
            "1"                                                     ///
            "0"                                                     ///
            "0"                                                     ///
            "`n_quantiles'"                                         ///
            "`allow_missing_x'"                                     ///
            "0"                                                     ///
            "0"                                                     ///
            "`quantile_csv'"

        /* Clear predictions for training obs */
        foreach q of local quantiles {
            local qint = round(`q' * 100)
            quietly replace `generate'_q`qint' = . if _n <= `n_train'
        }
    }

    /* ----------------------------------------------------------------
     * PROBABILITY FOREST predict
     * ----------------------------------------------------------------
     * Pass: X1..Xp Y_copy out_c0 out_c1 ...
     * ---------------------------------------------------------------- */
    else if "`forest_type'" == "probability" {

        display as text "Predicting with probability forest ..."

        /* Build output varlist */
        local output_vars ""
        forvalues c = 0/`=`n_classes'-1' {
            local output_vars `output_vars' `generate'_c`c'
        }

        /* Create tempvar copy of depvar with missing filled for test obs */
        tempvar y_safe
        quietly gen double `y_safe' = `depvar'
        quietly replace `y_safe' = 0 if _n > `n_train' & missing(`y_safe')

        /* Call plugin */
        plugin call grf_plugin `indepvars' `y_safe' `output_vars', ///
            "probability"                                           ///
            "`n_trees'"                                             ///
            "`seed'"                                                ///
            "`mtry'"                                                ///
            "`min_node'"                                            ///
            "`samplefrac'"                                          ///
            "`do_honesty'"                                          ///
            "`honestyfrac'"                                         ///
            "`do_honesty_prune'"                                    ///
            "`alpha'"                                               ///
            "`imbalancepenalty'"                                     ///
            "`cigroupsize'"                                         ///
            "`numthreads'"                                          ///
            "0"                                                     ///
            "`n_train'"                                             ///
            "`nindep'"                                              ///
            "1"                                                     ///
            "0"                                                     ///
            "0"                                                     ///
            "`n_classes'"                                            ///
            "`allow_missing_x'"                                     ///
            "0"                                                     ///
            "0"                                                     ///
            "`n_classes'"

        /* Clear predictions for training obs */
        forvalues c = 0/`=`n_classes'-1' {
            quietly replace `generate'_c`c' = . if _n <= `n_train'
        }
    }

    /* ----------------------------------------------------------------
     * INSTRUMENTAL FOREST predict
     * ----------------------------------------------------------------
     * Like causal but with 3 nuisance regressions: Y~X, W~X, Z~X
     * ---------------------------------------------------------------- */
    else if "`forest_type'" == "instrumental" {

        /* --- Step 1: Nuisance model Y ~ X (OOB on training only) --- */
        display as text "Step 1/4: Nuisance model Y ~ X (OOB on training) ..."

        tempvar yhat_pred
        quietly gen double `yhat_pred' = .

        plugin call grf_plugin `indepvars' `depvar' `yhat_pred' ///
            if _n <= `n_train',                                  ///
            "regression"                                          ///
            "`n_trees'"                                           ///
            "`seed'"                                              ///
            "`mtry'"                                              ///
            "`min_node'"                                          ///
            "`samplefrac'"                                        ///
            "`do_honesty'"                                        ///
            "`honestyfrac'"                                       ///
            "`do_honesty_prune'"                                  ///
            "`alpha'"                                             ///
            "`imbalancepenalty'"                                   ///
            "`cigroupsize'"                                       ///
            "`numthreads'"                                        ///
            "0"                                                   ///
            "0"                                                   ///
            "`nindep'"                                            ///
            "1"                                                   ///
            "0"                                                   ///
            "0"                                                   ///
            "1"                                                   ///
            "`allow_missing_x'"                                   ///
            "0"                                                   ///
            "0"

        /* --- Step 2: Nuisance model W ~ X (OOB on training only) --- */
        display as text "Step 2/4: Nuisance model W ~ X (OOB on training) ..."

        tempvar what_pred
        quietly gen double `what_pred' = .

        plugin call grf_plugin `indepvars' `treatvar' `what_pred' ///
            if _n <= `n_train',                                    ///
            "regression"                                           ///
            "`n_trees'"                                            ///
            "`seed'"                                               ///
            "`mtry'"                                               ///
            "`min_node'"                                           ///
            "`samplefrac'"                                         ///
            "`do_honesty'"                                         ///
            "`honestyfrac'"                                        ///
            "`do_honesty_prune'"                                   ///
            "`alpha'"                                              ///
            "`imbalancepenalty'"                                    ///
            "`cigroupsize'"                                        ///
            "`numthreads'"                                         ///
            "0"                                                    ///
            "0"                                                    ///
            "`nindep'"                                             ///
            "1"                                                    ///
            "0"                                                    ///
            "0"                                                    ///
            "1"                                                    ///
            "`allow_missing_x'"                                    ///
            "0"                                                    ///
            "0"

        /* --- Step 3: Nuisance model Z ~ X (OOB on training only) --- */
        display as text "Step 3/4: Nuisance model Z ~ X (OOB on training) ..."

        tempvar zhat_pred
        quietly gen double `zhat_pred' = .

        plugin call grf_plugin `indepvars' `instrvar' `zhat_pred' ///
            if _n <= `n_train',                                    ///
            "regression"                                           ///
            "`n_trees'"                                            ///
            "`seed'"                                               ///
            "`mtry'"                                               ///
            "`min_node'"                                           ///
            "`samplefrac'"                                         ///
            "`do_honesty'"                                         ///
            "`honestyfrac'"                                        ///
            "`do_honesty_prune'"                                   ///
            "`alpha'"                                              ///
            "`imbalancepenalty'"                                    ///
            "`cigroupsize'"                                        ///
            "`numthreads'"                                         ///
            "0"                                                    ///
            "0"                                                    ///
            "`nindep'"                                             ///
            "1"                                                    ///
            "0"                                                    ///
            "0"                                                    ///
            "1"                                                    ///
            "`allow_missing_x'"                                    ///
            "0"                                                    ///
            "0"

        /* --- Step 4: Center and run instrumental forest --- */
        display as text "Step 4/4: Instrumental forest on centered data (with predict) ..."

        tempvar y_centered w_centered z_centered

        /* Center training obs using OOB estimates */
        quietly gen double `y_centered' = `depvar'    - `yhat_pred' if _n <= `n_train'
        quietly gen double `w_centered' = `treatvar'  - `what_pred' if _n <= `n_train'
        quietly gen double `z_centered' = `instrvar'  - `zhat_pred' if _n <= `n_train'

        /* For test obs: fill with 0 placeholder */
        quietly replace `y_centered' = 0 if _n > `n_train'
        quietly replace `w_centered' = 0 if _n > `n_train'
        quietly replace `z_centered' = 0 if _n > `n_train'

        /* Call instrumental forest with n_train */
        plugin call grf_plugin `indepvars' `y_centered' `w_centered' ///
            `z_centered' `generate',                                  ///
            "instrumental"                                            ///
            "`n_trees'"                                               ///
            "`seed'"                                                  ///
            "`mtry'"                                                  ///
            "`min_node'"                                              ///
            "`samplefrac'"                                            ///
            "`do_honesty'"                                            ///
            "`honestyfrac'"                                           ///
            "`do_honesty_prune'"                                      ///
            "`alpha'"                                                 ///
            "`imbalancepenalty'"                                       ///
            "`cigroupsize'"                                           ///
            "`numthreads'"                                            ///
            "0"                                                       ///
            "`n_train'"                                               ///
            "`nindep'"                                                ///
            "1"                                                       ///
            "1"                                                       ///
            "1"                                                       ///
            "1"                                                       ///
            "`allow_missing_x'"                                       ///
            "0"                                                       ///
            "0"                                                       ///
            "`reduced_form_wt'"                                       ///
            "`do_stabilize'"

        /* Clear predictions for training obs */
        quietly replace `generate' = . if _n <= `n_train'
    }

    /* ----------------------------------------------------------------
     * SURVIVAL FOREST predict
     * ----------------------------------------------------------------
     * Pass: X1..Xp time_safe status_safe out_s1 ... out_sN
     * ---------------------------------------------------------------- */
    else if "`forest_type'" == "survival" {

        display as text "Predicting with survival forest ..."

        /* Build output varlist */
        local output_vars ""
        forvalues j = 1/`n_output_sv' {
            local output_vars `output_vars' `generate'_s`j'
        }

        /* Create tempvar copies of time and status, fill test obs with 0 */
        tempvar time_safe status_safe
        quietly gen double `time_safe' = `timevar'
        quietly replace `time_safe' = 0 if _n > `n_train' & missing(`time_safe')
        quietly gen double `status_safe' = `statusvar'
        quietly replace `status_safe' = 0 if _n > `n_train' & missing(`status_safe')

        /* Remap user-facing predtype to C++ convention (inverted) */
        local cpp_predtype = 1 - `pred_type'

        /* Call plugin: X1..Xp time_safe status_safe out1..outN */
        plugin call grf_plugin `indepvars' `time_safe' `status_safe' `output_vars', ///
            "survival"                                                               ///
            "`n_trees'"                                                              ///
            "`seed'"                                                                 ///
            "`mtry'"                                                                 ///
            "`min_node'"                                                             ///
            "`samplefrac'"                                                           ///
            "`do_honesty'"                                                           ///
            "`honestyfrac'"                                                          ///
            "`do_honesty_prune'"                                                     ///
            "`alpha'"                                                                ///
            "`imbalancepenalty'"                                                      ///
            "`cigroupsize'"                                                          ///
            "`numthreads'"                                                           ///
            "0"                                                                      ///
            "`n_train'"                                                              ///
            "`nindep'"                                                               ///
            "1"                                                                      ///
            "0"                                                                      ///
            "0"                                                                      ///
            "`n_output_sv'"                                                          ///
            "`allow_missing_x'"                                                      ///
            "0"                                                                      ///
            "0"                                                                      ///
            "0"                                                                      ///
            "`cpp_predtype'"

        /* Clear predictions for training obs */
        forvalues j = 1/`n_output_sv' {
            quietly replace `generate'_s`j' = . if _n <= `n_train'
        }
    }

    /* ----------------------------------------------------------------
     * CAUSAL SURVIVAL FOREST predict
     * ----------------------------------------------------------------
     * Nuisance pipeline:
     *   Step 1: W.hat via regression forest on training
     *   Step 2: W.centered = W - W.hat; compute simplified IPCW numer/denom
     *   Step 3: Call "causal_survival" plugin
     * Variable order: X1..Xp time w_centered status numer denom output
     * ---------------------------------------------------------------- */
    else if "`forest_type'" == "causal_survival" {

        /* --- Step 1: Propensity W.hat via regression on training --- */
        display as text "Step 1/3: Propensity model W ~ X (OOB on training) ..."

        tempvar what_pred
        quietly gen double `what_pred' = .

        plugin call grf_plugin `indepvars' `treatvar' `what_pred' ///
            if _n <= `n_train',                                    ///
            "regression"                                           ///
            "`n_trees'"                                            ///
            "`seed'"                                               ///
            "`mtry'"                                               ///
            "`min_node'"                                           ///
            "`samplefrac'"                                         ///
            "`do_honesty'"                                         ///
            "`honestyfrac'"                                        ///
            "`do_honesty_prune'"                                   ///
            "`alpha'"                                              ///
            "`imbalancepenalty'"                                    ///
            "`cigroupsize'"                                        ///
            "`numthreads'"                                         ///
            "0"                                                    ///
            "0"                                                    ///
            "`nindep'"                                             ///
            "1"                                                    ///
            "0"                                                    ///
            "0"                                                    ///
            "1"                                                    ///
            "`allow_missing_x'"                                    ///
            "0"                                                    ///
            "0"

        /* --- Step 2: Center W and compute simplified IPCW nuisance --- */
        display as text "Step 2/3: Centering W and computing nuisance estimates ..."

        tempvar w_cent cs_numer cs_denom time_safe status_safe

        /* Center treatment on training, fill test with 0 */
        quietly gen double `w_cent' = `treatvar' - `what_pred' if _n <= `n_train'
        quietly replace `w_cent' = 0 if _n > `n_train'

        /* Simplified IPCW nuisance (same as estimation command) */
        quietly gen double `cs_numer' = `w_cent' * `statusvar' * ///
            min(`timevar', `cs_horizon') if _n <= `n_train'
        quietly replace `cs_numer' = 0 if _n > `n_train'

        quietly gen double `cs_denom' = 1 if _n <= `n_train'
        quietly replace `cs_denom' = 0 if _n > `n_train'

        /* Tempvar copies of time and status, fill test with 0 */
        quietly gen double `time_safe' = `timevar'
        quietly replace `time_safe' = 0 if _n > `n_train' & missing(`time_safe')
        quietly gen double `status_safe' = `statusvar'
        quietly replace `status_safe' = 0 if _n > `n_train' & missing(`status_safe')

        /* --- Step 3: Causal survival forest with predict --- */
        display as text "Step 3/3: Causal survival forest (with predict) ..."

        plugin call grf_plugin `indepvars' `time_safe' `w_cent' ///
            `status_safe' `cs_numer' `cs_denom' `generate',      ///
            "causal_survival"                                     ///
            "`n_trees'"                                           ///
            "`seed'"                                              ///
            "`mtry'"                                              ///
            "`min_node'"                                          ///
            "`samplefrac'"                                        ///
            "`do_honesty'"                                        ///
            "`honestyfrac'"                                       ///
            "`do_honesty_prune'"                                  ///
            "`alpha'"                                             ///
            "`imbalancepenalty'"                                    ///
            "`cigroupsize'"                                       ///
            "`numthreads'"                                        ///
            "0"                                                   ///
            "`n_train'"                                           ///
            "`nindep'"                                            ///
            "1"                                                   ///
            "1"                                                   ///
            "0"                                                   ///
            "1"                                                   ///
            "`allow_missing_x'"                                   ///
            "0"                                                   ///
            "0"                                                   ///
            "`do_stabilize'"                                      ///
            "`=`nindep'+3'"                                       ///
            "`=`nindep'+4'"                                       ///
            "`=`nindep'+2'"                                       ///
            "`cs_target'"

        /* Clear predictions for training obs */
        quietly replace `generate' = . if _n <= `n_train'
    }

    /* ----------------------------------------------------------------
     * MULTI-ARM CAUSAL FOREST predict
     * ----------------------------------------------------------------
     * Nuisance pipeline (same pattern as estimation):
     *   Step 1: Y.hat via regression on training
     *   Step 2: For each arm k, W_k.hat via regression on training
     *   Step 3: Center Y.c = Y - Y.hat, W_k.c = W_k - W_k.hat
     *   Step 4: Call "multi_arm_causal" plugin
     * ---------------------------------------------------------------- */
    else if "`forest_type'" == "multi_causal" {

        local total_steps = `n_treat' + 2

        /* --- Step 1: Nuisance model Y ~ X (OOB on training) --- */
        display as text "Step 1/`total_steps': Nuisance model Y ~ X (OOB on training) ..."

        tempvar yhat_pred
        quietly gen double `yhat_pred' = .

        plugin call grf_plugin `indepvars' `depvar' `yhat_pred' ///
            if _n <= `n_train',                                  ///
            "regression"                                          ///
            "`n_trees'"                                           ///
            "`seed'"                                              ///
            "`mtry'"                                              ///
            "`min_node'"                                          ///
            "`samplefrac'"                                        ///
            "`do_honesty'"                                        ///
            "`honestyfrac'"                                       ///
            "`do_honesty_prune'"                                  ///
            "`alpha'"                                             ///
            "`imbalancepenalty'"                                   ///
            "`cigroupsize'"                                       ///
            "`numthreads'"                                        ///
            "0"                                                   ///
            "0"                                                   ///
            "`nindep'"                                            ///
            "1"                                                   ///
            "0"                                                   ///
            "0"                                                   ///
            "1"                                                   ///
            "`allow_missing_x'"                                   ///
            "0"                                                   ///
            "0"

        /* --- Step 2: For each treatment arm, W_k ~ X --- */
        local w_centered_vars ""
        forvalues j = 1/`n_treat' {
            local step = `j' + 1
            local tv : word `j' of `treatvars'
            display as text "Step `step'/`total_steps': Nuisance model `tv' ~ X (OOB on training) ..."

            tempvar what_`j'
            quietly gen double `what_`j'' = .

            plugin call grf_plugin `indepvars' `tv' `what_`j'' ///
                if _n <= `n_train',                              ///
                "regression"                                     ///
                "`n_trees'"                                      ///
                "`seed'"                                         ///
                "`mtry'"                                         ///
                "`min_node'"                                     ///
                "`samplefrac'"                                   ///
                "`do_honesty'"                                   ///
                "`honestyfrac'"                                  ///
                "`do_honesty_prune'"                             ///
                "`alpha'"                                        ///
                "`imbalancepenalty'"                               ///
                "`cigroupsize'"                                  ///
                "`numthreads'"                                   ///
                "0"                                              ///
                "0"                                              ///
                "`nindep'"                                       ///
                "1"                                              ///
                "0"                                              ///
                "0"                                              ///
                "1"                                              ///
                "`allow_missing_x'"                              ///
                "0"                                              ///
                "0"

            /* Center on training, fill test with 0 */
            tempvar wc_`j'
            quietly gen double `wc_`j'' = `tv' - `what_`j'' if _n <= `n_train'
            quietly replace `wc_`j'' = 0 if _n > `n_train'
            local w_centered_vars `w_centered_vars' `wc_`j''
        }

        /* --- Center Y --- */
        tempvar y_centered
        quietly gen double `y_centered' = `depvar' - `yhat_pred' if _n <= `n_train'
        quietly replace `y_centered' = 0 if _n > `n_train'

        /* --- Build output varlist --- */
        local output_vars ""
        forvalues j = 1/`n_treat' {
            local output_vars `output_vars' `generate'_t`j'
        }

        /* --- Step final: Multi-arm causal forest with predict --- */
        local final_step = `n_treat' + 2
        display as text "Step `final_step'/`final_step': Multi-arm causal forest (with predict) ..."

        plugin call grf_plugin `indepvars' `y_centered' `w_centered_vars' `output_vars', ///
            "multi_arm_causal"                                                            ///
            "`n_trees'"                                                                   ///
            "`seed'"                                                                      ///
            "`mtry'"                                                                      ///
            "`min_node'"                                                                  ///
            "`samplefrac'"                                                                ///
            "`do_honesty'"                                                                ///
            "`honestyfrac'"                                                               ///
            "`do_honesty_prune'"                                                          ///
            "`alpha'"                                                                     ///
            "`imbalancepenalty'"                                                            ///
            "`cigroupsize'"                                                               ///
            "`numthreads'"                                                                ///
            "0"                                                                           ///
            "`n_train'"                                                                   ///
            "`nindep'"                                                                    ///
            "1"                                                                           ///
            "`n_treat'"                                                                   ///
            "0"                                                                           ///
            "`n_treat'"                                                                   ///
            "`allow_missing_x'"                                                           ///
            "0"                                                                           ///
            "0"                                                                           ///
            "`do_stabilize'"                                                              ///
            "`n_treat'"

        /* Clear predictions for training obs */
        forvalues j = 1/`n_treat' {
            quietly replace `generate'_t`j' = . if _n <= `n_train'
        }
    }

    /* ----------------------------------------------------------------
     * MULTI-REGRESSION FOREST predict
     * ----------------------------------------------------------------
     * Pass: X1..Xp Y1_safe Y2_safe ... out_y1 out_y2 ...
     * ---------------------------------------------------------------- */
    else if "`forest_type'" == "multi_regression" {

        display as text "Predicting with multi-regression forest ..."

        /* Create tempvar copies of all depvars, fill test obs with 0 */
        local safe_depvars ""
        forvalues j = 1/`n_outcomes' {
            local dv : word `j' of `depvars'
            tempvar ysafe_`j'
            quietly gen double `ysafe_`j'' = `dv'
            quietly replace `ysafe_`j'' = 0 if _n > `n_train' & missing(`ysafe_`j'')
            local safe_depvars `safe_depvars' `ysafe_`j''
        }

        /* Build output varlist */
        local output_vars ""
        forvalues j = 1/`n_outcomes' {
            local output_vars `output_vars' `generate'_y`j'
        }

        /* Call plugin */
        plugin call grf_plugin `indepvars' `safe_depvars' `output_vars', ///
            "multi_regression"                                            ///
            "`n_trees'"                                                   ///
            "`seed'"                                                      ///
            "`mtry'"                                                      ///
            "`min_node'"                                                  ///
            "`samplefrac'"                                                ///
            "`do_honesty'"                                                ///
            "`honestyfrac'"                                               ///
            "`do_honesty_prune'"                                          ///
            "`alpha'"                                                     ///
            "`imbalancepenalty'"                                           ///
            "`cigroupsize'"                                               ///
            "`numthreads'"                                                ///
            "0"                                                           ///
            "`n_train'"                                                   ///
            "`nindep'"                                                    ///
            "`n_outcomes'"                                                ///
            "0"                                                           ///
            "0"                                                           ///
            "`n_outcomes'"                                                ///
            "`allow_missing_x'"                                           ///
            "0"                                                           ///
            "0"                                                           ///
            "`n_outcomes'"

        /* Clear predictions for training obs */
        forvalues j = 1/`n_outcomes' {
            quietly replace `generate'_y`j' = . if _n <= `n_train'
        }
    }

    /* ----------------------------------------------------------------
     * LOCAL LINEAR REGRESSION FOREST predict
     * ----------------------------------------------------------------
     * Pass: X1..Xp Y_safe output
     * ---------------------------------------------------------------- */
    else if "`forest_type'" == "ll_regression" {

        display as text "Predicting with local linear regression forest ..."

        /* Create tempvar copy of depvar with missing filled for test obs */
        tempvar y_safe
        quietly gen double `y_safe' = `depvar'
        quietly replace `y_safe' = 0 if _n > `n_train' & missing(`y_safe')

        /* Call plugin */
        plugin call grf_plugin `indepvars' `y_safe' `generate', ///
            "ll_regression"                                      ///
            "`n_trees'"                                          ///
            "`seed'"                                             ///
            "`mtry'"                                             ///
            "`min_node'"                                         ///
            "`samplefrac'"                                       ///
            "`do_honesty'"                                       ///
            "`honestyfrac'"                                      ///
            "`do_honesty_prune'"                                 ///
            "`alpha'"                                            ///
            "`imbalancepenalty'"                                  ///
            "`cigroupsize'"                                      ///
            "`numthreads'"                                       ///
            "0"                                                  ///
            "`n_train'"                                          ///
            "`nindep'"                                           ///
            "1"                                                  ///
            "0"                                                  ///
            "0"                                                  ///
            "1"                                                  ///
            "`allow_missing_x'"                                  ///
            "0"                                                  ///
            "0"                                                  ///
            "`enable_ll_split'"                                  ///
            "`ll_lambda'"                                        ///
            "`ll_weight_penalty'"                                ///
            "`ll_split_cutoff'"

        /* Clear predictions for training obs */
        quietly replace `generate' = . if _n <= `n_train'
    }

    /* ----------------------------------------------------------------
     * LINEAR MODEL FOREST predict
     * ----------------------------------------------------------------
     * Nuisance pipeline (same pattern as estimation):
     *   Step 1: Y.hat via regression on training using indepvars (xvars)
     *   Step 2: For each regressor k, W_k.hat via regression on training
     *   Step 3: Center Y.c, W_k.c on training; 0 for test
     *   Step 4: Call "lm_forest" plugin
     * ---------------------------------------------------------------- */
    else if "`forest_type'" == "lm_forest" {

        local total_steps = `n_regressors' + 2

        /* --- Step 1: Nuisance model Y ~ X (OOB on training) --- */
        display as text "Step 1/`total_steps': Nuisance model Y ~ X (OOB on training) ..."

        tempvar yhat_pred
        quietly gen double `yhat_pred' = .

        plugin call grf_plugin `indepvars' `depvar' `yhat_pred' ///
            if _n <= `n_train',                                  ///
            "regression"                                          ///
            "`n_trees'"                                           ///
            "`seed'"                                              ///
            "`mtry'"                                              ///
            "`min_node'"                                          ///
            "`samplefrac'"                                        ///
            "`do_honesty'"                                        ///
            "`honestyfrac'"                                       ///
            "`do_honesty_prune'"                                  ///
            "`alpha'"                                             ///
            "`imbalancepenalty'"                                   ///
            "`cigroupsize'"                                       ///
            "`numthreads'"                                        ///
            "0"                                                   ///
            "0"                                                   ///
            "`nindep'"                                            ///
            "1"                                                   ///
            "0"                                                   ///
            "0"                                                   ///
            "1"                                                   ///
            "`allow_missing_x'"                                   ///
            "0"                                                   ///
            "0"

        /* --- Step 2: For each regressor, W_k ~ X --- */
        local w_centered_vars ""
        forvalues j = 1/`n_regressors' {
            local step = `j' + 1
            local wv : word `j' of `regvars'
            display as text "Step `step'/`total_steps': Nuisance model `wv' ~ X (OOB on training) ..."

            tempvar what_`j'
            quietly gen double `what_`j'' = .

            plugin call grf_plugin `indepvars' `wv' `what_`j'' ///
                if _n <= `n_train',                              ///
                "regression"                                     ///
                "`n_trees'"                                      ///
                "`seed'"                                         ///
                "`mtry'"                                         ///
                "`min_node'"                                     ///
                "`samplefrac'"                                   ///
                "`do_honesty'"                                   ///
                "`honestyfrac'"                                  ///
                "`do_honesty_prune'"                             ///
                "`alpha'"                                        ///
                "`imbalancepenalty'"                               ///
                "`cigroupsize'"                                  ///
                "`numthreads'"                                   ///
                "0"                                              ///
                "0"                                              ///
                "`nindep'"                                       ///
                "1"                                              ///
                "0"                                              ///
                "0"                                              ///
                "1"                                              ///
                "`allow_missing_x'"                              ///
                "0"                                              ///
                "0"

            /* Center on training, fill test with 0 */
            tempvar wc_`j'
            quietly gen double `wc_`j'' = `wv' - `what_`j'' if _n <= `n_train'
            quietly replace `wc_`j'' = 0 if _n > `n_train'
            local w_centered_vars `w_centered_vars' `wc_`j''
        }

        /* --- Center Y --- */
        tempvar y_centered
        quietly gen double `y_centered' = `depvar' - `yhat_pred' if _n <= `n_train'
        quietly replace `y_centered' = 0 if _n > `n_train'

        /* --- Build output varlist --- */
        local output_vars ""
        forvalues j = 1/`n_regressors' {
            local output_vars `output_vars' `generate'_`j'
        }

        /* --- Step final: Linear model forest with predict --- */
        local final_step = `n_regressors' + 2
        display as text "Step `final_step'/`final_step': Linear model forest (with predict) ..."

        plugin call grf_plugin `indepvars' `y_centered' `w_centered_vars' `output_vars', ///
            "lm_forest"                                                                   ///
            "`n_trees'"                                                                   ///
            "`seed'"                                                                      ///
            "`mtry'"                                                                      ///
            "`min_node'"                                                                  ///
            "`samplefrac'"                                                                ///
            "`do_honesty'"                                                                ///
            "`honestyfrac'"                                                               ///
            "`do_honesty_prune'"                                                          ///
            "`alpha'"                                                                     ///
            "`imbalancepenalty'"                                                            ///
            "`cigroupsize'"                                                               ///
            "`numthreads'"                                                                ///
            "0"                                                                           ///
            "`n_train'"                                                                   ///
            "`nindep'"                                                                    ///
            "1"                                                                           ///
            "`n_regressors'"                                                              ///
            "0"                                                                           ///
            "`n_regressors'"                                                              ///
            "`allow_missing_x'"                                                           ///
            "0"                                                                           ///
            "0"                                                                           ///
            "`do_stabilize'"

        /* Clear predictions for training obs */
        forvalues j = 1/`n_regressors' {
            quietly replace `generate'_`j' = . if _n <= `n_train'
        }
    }

    /* ================================================================
     * Step 7: Summary statistics
     * ================================================================ */

    if "`forest_type'" == "regression" | "`forest_type'" == "ll_regression" {

        quietly summarize `generate' if _n > `n_train'
        local n_pred = r(N)
        local pred_mean = r(mean)
        local pred_sd = r(sd)

        display as text ""
        display as text "Predictions written to: " as result "`generate'"
        display as text "  Test obs predicted: " as result `n_pred'
        display as text "  Mean:               " as result %9.4f `pred_mean'
        display as text "  SD:                 " as result %9.4f `pred_sd'
        display as text ""

        return scalar N_test    = `n_pred'
        return scalar mean      = `pred_mean'
        return scalar sd        = `pred_sd'
        return local  predict_var "`generate'"
    }

    else if "`forest_type'" == "causal" | ///
            "`forest_type'" == "instrumental" {

        quietly summarize `generate' if _n > `n_train'
        local n_pred = r(N)
        local pred_mean = r(mean)
        local pred_sd = r(sd)

        local label = cond("`forest_type'" == "causal", "Causal", "Instrumental")
        local effect_label = cond("`forest_type'" == "causal", "CATE", "LATE")

        display as text ""
        display as text "`label' Forest Predictions"
        display as text "{hline 55}"
        display as text "`effect_label' predictions:     " as result "`generate'"
        display as text "  Test obs predicted: " as result `n_pred'
        display as text "  Mean `effect_label':          " as result %9.4f `pred_mean'
        display as text "  SD `effect_label':            " as result %9.4f `pred_sd'
        display as text "{hline 55}"
        display as text ""

        return scalar N_test    = `n_pred'
        return scalar mean_cate = `pred_mean'
        return scalar sd_cate   = `pred_sd'
        return local  predict_var "`generate'"
    }

    else if "`forest_type'" == "causal_survival" {

        quietly summarize `generate' if _n > `n_train'
        local n_pred = r(N)
        local pred_mean = r(mean)
        local pred_sd = r(sd)

        display as text ""
        display as text "Causal Survival Forest Predictions"
        display as text "{hline 55}"
        display as text "CATE predictions:       " as result "`generate'"
        display as text "  Test obs predicted: " as result `n_pred'
        display as text "  Mean CATE:          " as result %9.4f `pred_mean'
        display as text "  SD CATE:            " as result %9.4f `pred_sd'
        display as text "{hline 55}"
        display as text ""

        return scalar N_test    = `n_pred'
        return scalar mean_cate = `pred_mean'
        return scalar sd_cate   = `pred_sd'
        return local  predict_var "`generate'"
    }

    else if "`forest_type'" == "quantile" {

        display as text ""
        display as text "Quantile Forest Predictions"
        display as text "{hline 55}"
        foreach q of local quantiles {
            local qint = round(`q' * 100)
            quietly summarize `generate'_q`qint' if _n > `n_train'
            local n_pred = r(N)
            local pred_mean = r(mean)
            local pred_sd = r(sd)
            display as text "  q`qint' (`generate'_q`qint'): " ///
                as result "mean=" %9.4f `pred_mean' " sd=" %9.4f `pred_sd' ///
                " N=" `n_pred'
        }
        display as text "{hline 55}"
        display as text ""

        /* Return stats for first quantile */
        local q1 : word 1 of `quantiles'
        local qint1 = round(`q1' * 100)
        quietly summarize `generate'_q`qint1' if _n > `n_train'
        return scalar N_test = r(N)
        return local  predict_stub "`generate'"
    }

    else if "`forest_type'" == "probability" {

        display as text ""
        display as text "Probability Forest Predictions"
        display as text "{hline 55}"
        forvalues c = 0/`=`n_classes'-1' {
            quietly summarize `generate'_c`c' if _n > `n_train'
            local n_pred = r(N)
            local pred_mean = r(mean)
            local pred_sd = r(sd)
            display as text "  class `c' (`generate'_c`c'): " ///
                as result "mean=" %9.4f `pred_mean' " sd=" %9.4f `pred_sd' ///
                " N=" `n_pred'
        }
        display as text "{hline 55}"
        display as text ""

        quietly summarize `generate'_c0 if _n > `n_train'
        return scalar N_test = r(N)
        return local  predict_stub "`generate'"
    }

    else if "`forest_type'" == "survival" {

        display as text ""
        display as text "Survival Forest Predictions"
        display as text "{hline 55}"
        display as text "Predictions: " as result "`generate'_s1 ... `generate'_s`n_output_sv'"
        quietly summarize `generate'_s1 if _n > `n_train'
        local n_pred = r(N)
        display as text "  Test obs predicted: " as result `n_pred'
        display as text "  (`n_output_sv' columns at evaluated failure times)"
        quietly summarize `generate'_s1 if _n > `n_train'
        display as text "  First col mean:     " as result %9.4f r(mean)
        display as text "  First col SD:       " as result %9.4f r(sd)
        display as text "{hline 55}"
        display as text ""

        return scalar N_test = `n_pred'
        return local  predict_stub "`generate'"
    }

    else if "`forest_type'" == "multi_causal" {

        display as text ""
        display as text "Multi-Arm Causal Forest Predictions"
        display as text "{hline 55}"
        forvalues j = 1/`n_treat' {
            local tv : word `j' of `treatvars'
            quietly summarize `generate'_t`j' if _n > `n_train'
            local n_pred = r(N)
            local pred_mean = r(mean)
            local pred_sd = r(sd)
            display as text "  arm `j' (`tv'): " ///
                as result "mean=" %9.4f `pred_mean' " sd=" %9.4f `pred_sd' ///
                " N=" `n_pred'
        }
        display as text "{hline 55}"
        display as text ""

        quietly summarize `generate'_t1 if _n > `n_train'
        return scalar N_test = r(N)
        return local  predict_stub "`generate'"
    }

    else if "`forest_type'" == "multi_regression" {

        display as text ""
        display as text "Multi-Regression Forest Predictions"
        display as text "{hline 55}"
        forvalues j = 1/`n_outcomes' {
            local dv : word `j' of `depvars'
            quietly summarize `generate'_y`j' if _n > `n_train'
            local n_pred = r(N)
            local pred_mean = r(mean)
            local pred_sd = r(sd)
            display as text "  outcome `j' (`dv'): " ///
                as result "mean=" %9.4f `pred_mean' " sd=" %9.4f `pred_sd' ///
                " N=" `n_pred'
        }
        display as text "{hline 55}"
        display as text ""

        quietly summarize `generate'_y1 if _n > `n_train'
        return scalar N_test = r(N)
        return local  predict_stub "`generate'"
    }

    else if "`forest_type'" == "lm_forest" {

        display as text ""
        display as text "Linear Model Forest Predictions"
        display as text "{hline 55}"
        forvalues j = 1/`n_regressors' {
            local wv : word `j' of `regvars'
            quietly summarize `generate'_`j' if _n > `n_train'
            local n_pred = r(N)
            local pred_mean = r(mean)
            local pred_sd = r(sd)
            display as text "  coef `j' (`wv'): " ///
                as result "mean=" %9.4f `pred_mean' " sd=" %9.4f `pred_sd' ///
                " N=" `n_pred'
        }
        display as text "{hline 55}"
        display as text ""

        quietly summarize `generate'_1 if _n > `n_train'
        return scalar N_test = r(N)
        return local  predict_stub "`generate'"
    }

    return local  forest_type "`forest_type'"
    return scalar N_train   = `n_train'
end
