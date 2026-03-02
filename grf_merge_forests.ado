*! grf_merge_forests.ado -- Merge two prediction vectors into a composite score
*! Version 0.1.0

program define grf_merge_forests, rclass
    version 14.0

    syntax varlist(min=2 max=2 numeric) [if] [in], GENerate(name) [WEIghts(numlist min=2 max=2) REPlace]

    marksample touse

    gettoken p1 p2 : varlist
    confirm numeric variable `p1'
    confirm numeric variable `p2'

    local w1 = 0.5
    local w2 = 0.5
    if "`weights'" != "" {
        local n_w : word count `weights'
        if `n_w' != 2 {
            di as error "weights() must contain exactly 2 numbers"
            exit 198
        }
        local w1 : word 1 of `weights'
        local w2 : word 2 of `weights'
        if `w1' < 0 | `w2' < 0 {
            di as error "weights() must be non-negative"
            exit 198
        }
        if (`w1' + `w2') <= 0 {
            di as error "weights() sum must be positive"
            exit 198
        }
        local _sum = `w1' + `w2'
        local w1 = `w1' / `_sum'
        local w2 = `w2' / `_sum'
    }

    if "`replace'" != "" {
        capture drop `generate'
    }
    confirm new variable `generate'

    quietly gen double `generate' = `w1' * `p1' + `w2' * `p2' if `touse'
    label variable `generate' "Merged forest score: `w1'*`p1' + `w2'*`p2'"

    return scalar w1 = `w1'
    return scalar w2 = `w2'
    return local pred1 "`p1'"
    return local pred2 "`p2'"
    return local generate "`generate'"
end
