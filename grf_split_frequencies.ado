*! grf_split_frequencies.ado -- Split-frequency proxy via variable importance backend
*! Version 0.1.0

program define grf_split_frequencies, rclass
    version 14.0

    syntax varlist(min=2 numeric) [if] [in], [NTrees(integer 2000) SEED(integer 42) MAXDepth(integer 4) DECAYexponent(real 2.0)]

    quietly grf_variable_importance `varlist' `if' `in', ///
        ntrees(`ntrees') seed(`seed') maxdepth(`maxdepth') decayexponent(`decayexponent')

    tempname sf
    matrix `sf' = r(importance)
    matrix rownames `sf' = depth_aggregate

    di as text ""
    di as text "Split Frequencies (depth-aggregated proxy)"
    di as text "{hline 55}"
    di as text "This command currently returns a depth-aggregated split-frequency proxy"
    di as text "from variable importance weights."
    di as text "{hline 55}"

    return matrix split_frequencies = `sf'
    return scalar max_depth = `maxdepth'
    return scalar n_trees = `ntrees'
end
