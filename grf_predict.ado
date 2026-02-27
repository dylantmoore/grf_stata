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
    }

    if "`forest_type'" == "survival" {
        local timevar     "`e(timevar)'"
        local statusvar   "`e(statusvar)'"
        local n_output_sv = e(n_output)
        local pred_type   = e(pred_type)
    }

    if "`forest_type'" == "quantile" {
        local quantiles   "`e(quantiles)'"
        local n_quantiles = e(n_quantiles)
    }

    if "`forest_type'" == "probability" {
        local n_classes = e(n_classes)
    }

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

    if "`forest_type'" == "regression" | ///
       "`forest_type'" == "causal"     | ///
       "`forest_type'" == "instrumental" {

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

    cap program drop grf_plugin
    program grf_plugin, plugin using("grf_plugin_`c_os_'.plugin")

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
     * argv layout (common 20 args):
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
            "1"

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
            "1"

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
            "1"

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
            "`do_stabilize'"

        /* Clear predictions for training obs */
        quietly replace `generate' = . if _n <= `n_train'
    }

    /* ----------------------------------------------------------------
     * QUANTILE FOREST predict (stub -- not yet supported)
     * ---------------------------------------------------------------- */
    else if "`forest_type'" == "quantile" {
        display as error "grf_predict does not yet support quantile forests"
        display as error "this feature is planned for a future version"

        /* Clean up output variables we already created */
        foreach q of local quantiles {
            local qint = round(`q' * 100)
            capture drop `generate'_q`qint'
        }
        exit 198
    }

    /* ----------------------------------------------------------------
     * PROBABILITY FOREST predict (stub -- not yet supported)
     * ---------------------------------------------------------------- */
    else if "`forest_type'" == "probability" {
        display as error "grf_predict does not yet support probability forests"
        display as error "this feature is planned for a future version"

        /* Clean up output variables we already created */
        forvalues c = 0/`=`n_classes'-1' {
            capture drop `generate'_c`c'
        }
        exit 198
    }

    /* ----------------------------------------------------------------
     * INSTRUMENTAL FOREST predict (stub -- not yet supported)
     * ---------------------------------------------------------------- */
    else if "`forest_type'" == "instrumental" {
        display as error "grf_predict does not yet support instrumental forests"
        display as error "this feature is planned for a future version"

        capture drop `generate'
        exit 198
    }

    /* ----------------------------------------------------------------
     * SURVIVAL FOREST predict (stub -- not yet supported)
     * ---------------------------------------------------------------- */
    else if "`forest_type'" == "survival" {
        display as error "grf_predict does not yet support survival forests"
        display as error "this feature is planned for a future version"

        /* Clean up output variables we already created */
        forvalues j = 1/`n_output_sv' {
            capture drop `generate'_s`j'
        }
        exit 198
    }

    /* ================================================================
     * Step 7: Summary statistics
     * ================================================================ */

    if "`forest_type'" == "regression" {

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

    else if "`forest_type'" == "causal" {

        quietly summarize `generate' if _n > `n_train'
        local n_pred = r(N)
        local pred_mean = r(mean)
        local pred_sd = r(sd)

        display as text ""
        display as text "Causal Forest Predictions"
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

    return local  forest_type "`forest_type'"
    return scalar N_train   = `n_train'
end
