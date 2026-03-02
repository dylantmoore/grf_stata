*! grf_plot_tree.ado -- Plot split-frequency proxy chart
*! Version 0.1.0

program define grf_plot_tree, rclass
    version 14.0

    syntax varlist(min=2 numeric) [if] [in], [NTrees(integer 500) SEED(integer 42) MAXDepth(integer 4) NAME(name) REPlace]

    quietly grf_variable_importance `varlist' `if' `in', ntrees(`ntrees') seed(`seed') maxdepth(`maxdepth')

    tempname imp
    matrix `imp' = r(importance)
    local cn : colnames `imp'
    local k : word count `cn'

    preserve
    clear
    set obs `k'
    gen str64 variable = ""
    gen double importance = .

    local j = 0
    foreach v of local cn {
        local ++j
        replace variable = "`v'" in `j'
        replace importance = `imp'[1, `j'] in `j'
    }

    if "`name'" == "" {
        local name "grf_tree_proxy_plot"
    }

    graph bar importance, over(variable, sort(1) descending) ///
        title("GRF Split-Frequency Proxy") ///
        ytitle("Importance") name(`name', `replace')
    restore

    return matrix importance = `imp'
    return local graph_name "`name'"
end
