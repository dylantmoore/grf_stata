*! grf_expected_survival.ado -- Expected survival time E[T|X] from survival curves
*! Version 0.1.0
*! Computes E[T|X] = integral_0^max_t S(t|X) dt via trapezoidal integration
*! Pure Stata computation, no C++ plugin needed

program define grf_expected_survival, rclass
    version 14.0

    syntax , GENerate(name)        ///
        [                          ///
            REPlace                ///
            PRedictions(string)    ///
            GRID(numlist >0 sort)  ///
        ]

    /* ---- Verify prior estimation ---- */
    if "`e(cmd)'" != "grf_survival_forest" {
        display as error "grf_expected_survival requires prior estimation" ///
            " by grf_survival_forest"
        exit 301
    }

    /* ---- Handle replace ---- */
    if "`replace'" != "" {
        capture drop `generate'
    }
    confirm new variable `generate'

    /* ---- Determine prediction stub ---- */
    if "`predictions'" == "" {
        local predictions "`e(predict_stub)'"
    }
    if "`predictions'" == "" {
        display as error "could not determine survival curve variable stub"
        display as error "specify predictions(stub) where stub_s1, stub_s2, ... exist"
        exit 111
    }

    /* ---- Determine number of survival curve columns ---- */
    local n_output = e(n_output)
    if missing(`n_output') | `n_output' < 1 {
        /* Count how many stub_sN variables exist */
        local n_output 0
        local keep_looking 1
        while `keep_looking' {
            local next = `n_output' + 1
            capture confirm numeric variable `predictions'_s`next'
            if _rc {
                local keep_looking 0
            }
            else {
                local n_output = `next'
            }
        }
    }

    if `n_output' < 1 {
        display as error "no survival curve variables found"
        display as error "expected variables named `predictions'_s1, `predictions'_s2, ..."
        exit 111
    }

    /* ---- Confirm survival curve variables exist ---- */
    forvalues j = 1/`n_output' {
        confirm numeric variable `predictions'_s`j'
    }

    /* ---- Determine failure time grid ---- */
    if "`grid'" == "" {
        /* Try to read from e() results */
        local grid "`e(failure_times)'"
    }

    if "`grid'" == "" {
        /* Reconstruct from data: use unique sorted values of the time variable */
        local timevar "`e(timevar)'"
        local statusvar "`e(statusvar)'"

        if "`timevar'" == "" {
            display as error "could not determine failure time grid"
            display as error "specify grid(numlist) with the failure times"
            exit 111
        }

        /* Extract unique failure times (where event occurred) */
        tempvar is_event
        quietly gen byte `is_event' = (`statusvar' == 1) if !missing(`statusvar')

        /* Get unique sorted failure times */
        tempname grid_mat
        quietly tab `timevar' if `is_event' == 1, matrow(`grid_mat')
        local n_failures = rowsof(`grid_mat')

        if `n_failures' < 1 {
            display as error "no failure times found in the data"
            exit 198
        }

        /* Use up to n_output failure times */
        local n_grid = min(`n_failures', `n_output')
        local grid ""
        forvalues j = 1/`n_grid' {
            local t_j = `grid_mat'[`j', 1]
            local grid "`grid' `t_j'"
        }
    }

    /* ---- Count grid points and validate ---- */
    local n_grid 0
    foreach t of local grid {
        local n_grid = `n_grid' + 1
    }

    if `n_grid' != `n_output' {
        /* If grid has more points than columns, truncate grid */
        if `n_grid' > `n_output' {
            display as text "Warning: grid has `n_grid' points but only" ///
                " `n_output' survival curve columns"
            display as text "  Using first `n_output' grid points"
            local new_grid ""
            local gi 0
            foreach t of local grid {
                local gi = `gi' + 1
                if `gi' <= `n_output' {
                    local new_grid "`new_grid' `t'"
                }
            }
            local grid "`new_grid'"
            local n_grid = `n_output'
        }
        else if `n_grid' < `n_output' {
            display as error "grid has `n_grid' points but there are" ///
                " `n_output' survival curve columns"
            display as error "provide a grid with exactly `n_output' time points"
            exit 198
        }
    }

    /* ---- Store grid values in locals for indexing ---- */
    local gi 0
    foreach t of local grid {
        local gi = `gi' + 1
        local t_`gi' = `t'
    }

    /* ---- Identify sample ---- */
    tempvar touse
    quietly gen byte `touse' = !missing(`predictions'_s1)
    forvalues j = 2/`n_output' {
        quietly replace `touse' = 0 if missing(`predictions'_s`j')
    }
    quietly count if `touse'
    local n_use = r(N)

    if `n_use' < 1 {
        display as error "no observations with non-missing survival curves"
        exit 2000
    }

    /* ---- Display header ---- */
    display as text ""
    display as text "Expected Survival Time: E[T|X]"
    display as text "{hline 55}"
    display as text "Survival curves:       " as result ///
        "`predictions'_s1 ... `predictions'_s`n_output'"
    display as text "Grid points:           " as result `n_grid'
    display as text "Time range:            " as result ///
        "[`t_1', `t_`n_grid'']"
    display as text "Observations:          " as result `n_use'
    display as text "{hline 55}"
    display as text ""

    /* ---- Compute E[T|X] via trapezoidal integration ----
     *
     * E[T|X] = integral_0^t_max S(t|X) dt
     *
     * We approximate using the trapezoidal rule on the grid:
     *   E[T|X] ~ S(t1)*t1   (rectangle from 0 to t1, S=1 up to first time)
     *          + sum_{j=2}^{N} 0.5*(S(t_{j-1}) + S(t_j)) * (t_j - t_{j-1})
     *
     * Note: We assume S(0)=1, so the first interval [0, t_1] contributes
     *       0.5 * (1 + S(t_1)) * t_1
     */

    display as text "Computing expected survival times ..."

    /* Start with the interval [0, t_1]: S(0) = 1 */
    quietly gen double `generate' = 0.5 * (1 + `predictions'_s1) * `t_1' ///
        if `touse'

    /* Add trapezoidal areas for subsequent intervals */
    forvalues j = 2/`n_output' {
        local jm1 = `j' - 1
        local dt = `t_`j'' - `t_`jm1''
        quietly replace `generate' = `generate' ///
            + 0.5 * (`predictions'_s`jm1' + `predictions'_s`j') * `dt' ///
            if `touse'
    }

    label variable `generate' "Expected survival time E[T|X]"

    /* ---- Summary statistics ---- */
    quietly summarize `generate' if `touse'
    local n_computed = r(N)
    local est_mean = r(mean)
    local est_sd = r(sd)
    local est_min = r(min)
    local est_max = r(max)

    /* ---- Display results ---- */
    display as text ""
    display as text "Expected Survival Time Results"
    display as text "{hline 55}"
    display as text "Variable created:      " as result "`generate'"
    display as text "Observations:          " as result `n_computed'
    display as text ""
    display as text "Summary of E[T|X]:"
    display as text "  Mean:         " as result %12.4f `est_mean'
    display as text "  Std. Dev.:    " as result %12.4f `est_sd'
    display as text "  Min:          " as result %12.4f `est_min'
    display as text "  Max:          " as result %12.4f `est_max'
    display as text "{hline 55}"
    display as text ""

    /* ---- Store results ---- */
    return scalar N     = `n_computed'
    return scalar mean  = `est_mean'
    return scalar sd    = `est_sd'
    return scalar min   = `est_min'
    return scalar max   = `est_max'
    return local  generate "`generate'"
    return local  grid     "`grid'"
    return scalar n_grid = `n_grid'
end
