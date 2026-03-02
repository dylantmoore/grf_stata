*! grf_get_forest_weights.ado -- Proxy forest kernel weights from fitted predictions/features
*! Version 0.2.0

program define grf_get_forest_weights, rclass
    version 14.0

    syntax , OBS(integer) GENerate(name) ///
        [SCALE(real 1.0) XSCALE(real 1.0) PREDWEIGHT(real 1.0) XVARS(varlist numeric) REPlace]

    if `scale' <= 0 {
        di as error "scale() must be strictly positive"
        exit 198
    }
    if `xscale' <= 0 {
        di as error "xscale() must be strictly positive"
        exit 198
    }
    if `predweight' < 0 | `predweight' > 1 {
        di as error "predweight() must be in [0, 1]"
        exit 198
    }
    if "`xvars'" == "" & `predweight' < 1 {
        di as error "predweight() < 1 requires xvars()"
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

    tempvar touse
    quietly gen byte `touse' = !missing(`predvar')
    if "`xvars'" != "" {
        markout `touse' `xvars'
    }
    quietly count if `touse'
    if r(N) < 2 {
        di as error "need at least 2 complete observations for proxy-weight construction"
        exit 2000
    }

    if !`touse'[`obs'] {
        di as error "obs() points to missing prediction or missing xvars() values"
        exit 498
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

    tempvar pdist xdist dist wraw
    quietly gen double `pdist' = abs(`predvar' - `anchor') / `scale' if `touse'
    quietly gen double `xdist' = 0 if `touse'

    local uses_xvars = 0
    if "`xvars'" != "" {
        local uses_xvars = 1
        foreach x of local xvars {
            quietly summarize `x' if `touse'
            local sdx = r(sd)
            if `sdx' <= 1e-12 {
                local sdx = 1
            }
            local anchor_x = `x'[`obs']
            quietly replace `xdist' = `xdist' + ((`x' - `anchor_x') / `sdx')^2 if `touse'
        }
        quietly replace `xdist' = sqrt(`xdist') / `xscale' if `touse'
    }

    if "`xvars'" != "" {
        quietly gen double `dist' = (`predweight' * `pdist' + (1 - `predweight') * `xdist') if `touse'
    }
    else {
        quietly gen double `dist' = `pdist' if `touse'
    }
    quietly summarize `dist', meanonly
    local min_dist = r(min)
    quietly gen double `wraw' = exp(-(`dist' - `min_dist')) if `touse'
    quietly summarize `wraw', meanonly
    if r(sum) <= 0 {
        di as error "failed to construct nonzero proxy weights"
        exit 498
    }

    quietly gen double `generate' = `wraw' / r(sum) if !missing(`wraw')
    if "`xvars'" != "" {
        label variable `generate' "Proxy weights obs `obs' (predweight=`predweight', xvars)"
    }
    else {
        label variable `generate' "Proxy forest weights anchored at obs `obs'"
    }

    return scalar obs = `obs'
    return scalar scale = `scale'
    return scalar xscale = `xscale'
    return scalar predweight = `predweight'
    return scalar uses_xvars = `uses_xvars'
    return local generate "`generate'"
    return local predict_var "`predvar'"
    if "`xvars'" != "" return local xvars "`xvars'"
end
