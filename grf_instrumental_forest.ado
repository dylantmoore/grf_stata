*! grf_instrumental_forest.ado -- Instrumental Forest via grf C++ library
*! Version 0.1.0
*! Implements Athey, Tibshirani, Wager (2019) instrumental_forest()

program define grf_instrumental_forest, eclass
    version 14.0

    syntax varlist(min=4 numeric) [if] [in],  ///
        GENerate(name)                         ///
        [                                      ///
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
            CIGroupsize(integer 1)             ///
            NUMThreads(integer 0)              ///
            ESTIMATEVariance                   ///
            REPlace                            ///
            VARGenerate(name)                  ///
            REDucedformweight(real 0.0)        ///
            STABilizesplits                    ///
            NUISancetrees(integer 500)         ///
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

    /* ---- Parse estimate_variance ---- */
    local do_est_var 0
    if "`estimatevariance'" != "" {
        local do_est_var 1
        if `cigroupsize' < 2 {
            local cigroupsize 2
        }
    }

    /* ---- Parse stabilize_splits ---- */
    local do_stabilize 0
    if "`stabilizesplits'" != "" {
        local do_stabilize 1
    }

    /* ---- Handle replace ---- */
    if "`replace'" != "" {
        capture drop `generate'
        if `do_est_var' & "`vargenerate'" != "" {
            capture drop `vargenerate'
        }
    }
    confirm new variable `generate'

    /* ---- Parse varlist: Y W Z X1..Xp ---- */
    gettoken depvar rest : varlist
    gettoken treatment rest : rest
    gettoken instrument indepvars : rest
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

    /* ---- Display header ---- */
    display as text ""
    display as text "Generalized Random Forest: Instrumental Forest"
    display as text "{hline 55}"
    display as text "Dependent variable:    " as result "`depvar'"
    display as text "Treatment variable:    " as result "`treatment'"
    display as text "Instrument variable:   " as result "`instrument'"
    display as text "Predictors:            " as result "`indepvars'"
    display as text "Observations:          " as result `n_use'
    display as text "Trees:                 " as result `ntrees'
    display as text "Honesty:               " as result cond(`do_honesty', "yes", "no")
    display as text "Reduced form weight:   " as result `reducedformweight'
    if `do_stabilize' {
        display as text "Stabilize splits:      " as result "yes"
    }
    if `do_est_var' {
        display as text "Variance estimation:   " as result "yes (ci_group_size=`cigroupsize')"
    }
    display as text "{hline 55}"
    display as text ""

    /* ---- Nuisance pipeline: center Y, W, Z using regression forests ----
     *
     * Step 1: Y.hat = E[Y|X] via regression forest
     * Step 2: W.hat = E[W|X] via regression forest
     * Step 3: Z.hat = E[Z|X] via regression forest
     * Step 4: Y.c = Y - Y.hat, W.c = W - W.hat, Z.c = Z - Z.hat
     */

    display as text "Nuisance estimation (3 regression forests)..."

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

    /* ---- Step 1: Y.hat ---- */
    tempvar Y_hat
    quietly gen double `Y_hat' = .

    display as text "  Fitting Y ~ X..."
    plugin call grf_plugin `indepvars' `depvar' `Y_hat' ///
        if `touse',                                      ///
        "regression"                                     ///
        "`nuisancetrees'"                                ///
        "`seed'"                                         ///
        "`mtry'"                                         ///
        "`minnodesize'"                                  ///
        "`samplefrac'"                                   ///
        "`do_honesty'"                                   ///
        "`honestyfrac'"                                  ///
        "`do_honesty_prune'"                             ///
        "`alpha'"                                        ///
        "`imbalancepenalty'"                              ///
        "1"                                              ///
        "`numthreads'"                                   ///
        "0"                                              ///
        "0"                                              ///
        "`nindep'"                                       ///
        "1"                                              ///
        "0"                                              ///
        "0"                                              ///
        "1"

    /* ---- Step 2: W.hat ---- */
    tempvar W_hat
    quietly gen double `W_hat' = .

    display as text "  Fitting W ~ X..."
    plugin call grf_plugin `indepvars' `treatment' `W_hat' ///
        if `touse',                                         ///
        "regression"                                        ///
        "`nuisancetrees'"                                   ///
        "`seed'"                                            ///
        "`mtry'"                                            ///
        "`minnodesize'"                                     ///
        "`samplefrac'"                                      ///
        "`do_honesty'"                                      ///
        "`honestyfrac'"                                     ///
        "`do_honesty_prune'"                                ///
        "`alpha'"                                           ///
        "`imbalancepenalty'"                                 ///
        "1"                                                 ///
        "`numthreads'"                                      ///
        "0"                                                 ///
        "0"                                                 ///
        "`nindep'"                                          ///
        "1"                                                 ///
        "0"                                                 ///
        "0"                                                 ///
        "1"

    /* ---- Step 3: Z.hat ---- */
    tempvar Z_hat
    quietly gen double `Z_hat' = .

    display as text "  Fitting Z ~ X..."
    plugin call grf_plugin `indepvars' `instrument' `Z_hat' ///
        if `touse',                                          ///
        "regression"                                         ///
        "`nuisancetrees'"                                    ///
        "`seed'"                                             ///
        "`mtry'"                                             ///
        "`minnodesize'"                                      ///
        "`samplefrac'"                                       ///
        "`do_honesty'"                                       ///
        "`honestyfrac'"                                      ///
        "`do_honesty_prune'"                                 ///
        "`alpha'"                                            ///
        "`imbalancepenalty'"                                  ///
        "1"                                                  ///
        "`numthreads'"                                       ///
        "0"                                                  ///
        "0"                                                  ///
        "`nindep'"                                           ///
        "1"                                                  ///
        "0"                                                  ///
        "0"                                                  ///
        "1"

    /* ---- Step 4: Center Y, W, Z ---- */
    tempvar Y_centered W_centered Z_centered
    quietly gen double `Y_centered' = `depvar' - `Y_hat' if `touse'
    quietly gen double `W_centered' = `treatment' - `W_hat' if `touse'
    quietly gen double `Z_centered' = `instrument' - `Z_hat' if `touse'

    display as text "  Centering complete."
    display as text ""

    /* ---- Build output varlist ---- */
    local output_vars `generate'
    if `do_est_var' {
        local output_vars `generate' `vargenerate'
    }

    /* ---- Call plugin for instrumental forest ----
     *
     * Variable order: X1..Xp Y.c W.c Z.c out1 [out2]
     * argv: forest_type num_trees seed mtry min_node_size sample_fraction
     *       honesty honesty_fraction honesty_prune alpha imbalance_penalty
     *       ci_group_size num_threads estimate_variance compute_oob
     *       n_x n_y n_w n_z n_output reduced_form_weight stabilize_splits
     */
    plugin call grf_plugin `indepvars' `Y_centered' `W_centered' ///
        `Z_centered' `output_vars'                                ///
        if `touse',                                               ///
        "instrumental"                                            ///
        "`ntrees'"                                                ///
        "`seed'"                                                  ///
        "`mtry'"                                                  ///
        "`minnodesize'"                                           ///
        "`samplefrac'"                                            ///
        "`do_honesty'"                                            ///
        "`honestyfrac'"                                           ///
        "`do_honesty_prune'"                                      ///
        "`alpha'"                                                 ///
        "`imbalancepenalty'"                                       ///
        "`cigroupsize'"                                           ///
        "`numthreads'"                                            ///
        "`do_est_var'"                                            ///
        "0"                                                       ///
        "`nindep'"                                                ///
        "1"                                                       ///
        "1"                                                       ///
        "1"                                                       ///
        "`n_output'"                                              ///
        "`reducedformweight'"                                     ///
        "`do_stabilize'"

    /* ---- Store results ---- */
    ereturn clear
    ereturn scalar N                  = `n_use'
    ereturn scalar n_trees            = `ntrees'
    ereturn scalar seed               = `seed'
    ereturn scalar mtry               = `mtry'
    ereturn scalar min_node           = `minnodesize'
    ereturn scalar alpha              = `alpha'
    ereturn scalar honesty            = `do_honesty'
    ereturn scalar honesty_prune = `do_honesty_prune'
    ereturn scalar sample_fraction    = `samplefrac'
    ereturn scalar honesty_fraction   = `honestyfrac'
    ereturn scalar imbalance_penalty  = `imbalancepenalty'
    ereturn scalar ci_group_size      = `cigroupsize'
    ereturn scalar reduced_form_wt    = `reducedformweight'
    ereturn scalar stabilize_splits   = `do_stabilize'
    ereturn local  cmd                  "grf_instrumental_forest"
    ereturn local  forest_type          "instrumental"
    ereturn local  depvar               "`depvar'"
    ereturn local  treatment            "`treatment'"
    ereturn local  instrument           "`instrument'"
    ereturn local  indepvars            "`indepvars'"
    ereturn local  predict_var          "`generate'"
    if `do_est_var' {
        ereturn local variance_var "`vargenerate'"
    }

    /* ---- Summary stats ---- */
    quietly summarize `generate' if `touse'
    local n_pred = r(N)
    local pred_mean = r(mean)
    local pred_sd = r(sd)

    display as text ""
    display as text "LATE estimates written to: " as result "`generate'"
    display as text "  Non-missing:  " as result `n_pred'
    display as text "  Mean:         " as result %9.4f `pred_mean'
    display as text "  SD:           " as result %9.4f `pred_sd'
    if `do_est_var' {
        quietly summarize `vargenerate' if `touse'
        display as text "Variance estimates:     " as result "`vargenerate'"
        display as text "  Mean variance: " as result %9.6f r(mean)
    }
    display as text ""
end
