*! grf_get_tree.ado -- Retrieve tree metadata from latest GRF model
*! Version 0.1.0

program define grf_get_tree, rclass
    version 14.0

    syntax [, TREE(integer 1)]

    quietly grf_tree_summary, tree(`tree')
    return add
end
