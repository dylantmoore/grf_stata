*! grf_get_leaf_node.ado -- Pseudo leaf-node assignment from fitted predictions
*! Version 0.1.0

program define grf_get_leaf_node, rclass
    version 14.0

    syntax , GENerate(name) [GROUPS(integer 20) REPlace]

    if `groups' < 2 {
        di as error "groups() must be at least 2"
        exit 198
    }

    if "`e(forest_type)'" == "" {
        di as error "grf_get_leaf_node requires prior grf_* estimation"
        exit 301
    }

    local predvar "`e(predict_var)'"
    if "`predvar'" == "" {
        di as error "no prediction variable found in e(predict_var)"
        exit 498
    }
    confirm numeric variable `predvar'

    if "`replace'" != "" {
        capture drop `generate'
    }
    confirm new variable `generate'

    tempvar touse
    quietly gen byte `touse' = !missing(`predvar')
    quietly count if `touse'
    if r(N) < 2 {
        di as error "need at least 2 non-missing predictions"
        exit 2000
    }

    quietly xtile `generate' = `predvar' if `touse', nq(`groups')
    label variable `generate' "Pseudo leaf-node group from `predvar' (nq=`groups')"

    return scalar N = r(N)
    return scalar groups = `groups'
    return local generate "`generate'"
    return local predict_var "`predvar'"
end
