*! grf_tree_summary.ado -- Tree-level metadata summary for latest GRF model
*! Version 0.1.0

program define grf_tree_summary, rclass
    version 14.0

    syntax [, TREE(integer 1)]

    if "`e(forest_type)'" == "" {
        di as error "grf_tree_summary requires a prior grf_* estimation command"
        exit 301
    }

    local n_trees = .
    capture local n_trees = e(n_trees)
    if _rc local n_trees = .

    local model_id = .
    capture local model_id = e(model_id)
    if _rc local model_id = .

    if missing(`n_trees') {
        di as error "tree metadata unavailable: e(n_trees) not found"
        exit 498
    }

    if `tree' < 1 | `tree' > `n_trees' {
        di as error "tree() must be between 1 and `n_trees'"
        exit 198
    }

    di as text ""
    di as text "GRF Tree Summary"
    di as text "{hline 55}"
    di as text "Forest type:          " as result "`e(forest_type)'"
    if !missing(`model_id') di as text "Model id:             " as result `model_id'
    di as text "Tree index:           " as result `tree' " of " `n_trees'
    di as text "Note: detailed node-level persistence is not available in the current plugin architecture."
    di as text "{hline 55}"

    return scalar tree = `tree'
    return scalar n_trees = `n_trees'
    if !missing(`model_id') return scalar model_id = `model_id'
    return local forest_type "`e(forest_type)'"
end
