*! grf_forest_summary.ado -- Print summary metadata for latest GRF model
*! Version 0.2.0

program define grf_forest_summary, rclass
    version 14.0

    syntax [, ALL]

    if "`e(forest_type)'" == "" {
        di as error "grf_forest_summary requires a prior grf_* estimation command"
        exit 301
    }

    local cmd "`e(cmd)'"
    local forest_type "`e(forest_type)'"
    local depvar "`e(depvar)'"
    local indepvars "`e(indepvars)'"

    local N = .
    capture local N = e(N)
    if _rc local N = .

    local n_trees = .
    capture local n_trees = e(n_trees)
    if _rc local n_trees = .

    local seed = .
    capture local seed = e(seed)
    if _rc local seed = .

    local model_id = .
    capture local model_id = e(model_id)
    if _rc local model_id = .

    local escalars : e(scalars)
    local emacros : e(macros)

    di as text ""
    di as text "GRF Forest Summary"
    di as text "{hline 55}"
    di as text "Command:              " as result "`cmd'"
    di as text "Forest type:          " as result "`forest_type'"
    if !missing(`model_id') di as text "Model id:             " as result `model_id'
    if !missing(`N') di as text "Observations:         " as result `N'
    if !missing(`n_trees') di as text "Trees:                " as result `n_trees'
    if !missing(`seed') di as text "Seed:                 " as result `seed'
    if "`depvar'" != "" di as text "Outcome variable:     " as result "`depvar'"
    if "`indepvars'" != "" di as text "Predictors:           " as result "`indepvars'"
    if "`all'" != "" {
        di as text "{hline 55}"
        di as text "Stored e() scalars:   " as result cond("`escalars'" == "", "(none)", "`escalars'")
        foreach s of local escalars {
            capture noisily di as text "  e(`s') = " as result %12.6g e(`s')
        }
        di as text "Stored e() macros:    " as result cond("`emacros'" == "", "(none)", "`emacros'")
    }
    di as text "{hline 55}"

    return local cmd "`cmd'"
    return local forest_type "`forest_type'"
    return local e_scalars "`escalars'"
    return local e_macros "`emacros'"
    if "`depvar'" != "" return local depvar "`depvar'"
    if "`indepvars'" != "" return local indepvars "`indepvars'"
    if !missing(`model_id') return scalar model_id = `model_id'
    if !missing(`N') return scalar N = `N'
    if !missing(`n_trees') return scalar n_trees = `n_trees'
    if !missing(`seed') return scalar seed = `seed'
end
