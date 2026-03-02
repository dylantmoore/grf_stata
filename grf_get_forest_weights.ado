*! grf_get_forest_weights.ado -- Proxy forest kernel weights from fitted predictions
*! Version 0.1.0

program define grf_get_forest_weights, rclass
    version 14.0

    syntax , OBS(integer) GENerate(name) [SCALE(real 1.0) REPlace]

    if `scale' <= 0 {
        di as error "scale() must be strictly positive"
        exit 198
    }

    if "`e(forest_type)'" == "" {
        di as error "grf_get_forest_weights requires prior grf_* estimation"
        exit 301
    }

    local predvar "`e(predict_var)'"
    if "`predvar'" == "" {
        di as error "no prediction variable found in e(predict_var)"
        exit 498
    }
    confirm numeric variable `predvar'

    if `obs' < 1 | `obs' > _N {
        di as error "obs() must be between 1 and _N"
        exit 198
    }

    quietly summarize `predvar' in `obs'/`obs'
    if missing(r(mean)) {
        di as error "obs() points to missing prediction value"
        exit 498
    }
    local anchor = r(mean)

    if "`replace'" != "" {
        capture drop `generate'
    }
    confirm new variable `generate'

    tempvar wraw
    quietly gen double `wraw' = exp(-abs(`predvar' - `anchor') / `scale') if !missing(`predvar')
    quietly summarize `wraw', meanonly
    if r(sum) <= 0 {
        di as error "failed to construct nonzero proxy weights"
        exit 498
    }

    quietly gen double `generate' = `wraw' / r(sum) if !missing(`wraw')
    label variable `generate' "Proxy forest weights anchored at obs `obs'"

    return scalar obs = `obs'
    return scalar scale = `scale'
    return local generate "`generate'"
    return local predict_var "`predvar'"
end
