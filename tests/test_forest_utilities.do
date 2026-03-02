* test_forest_utilities.do -- Tests for forest introspection/utility commands
* Covers: grf_forest_summary, grf_tree_summary, grf_get_tree,
*         grf_get_leaf_node, grf_get_forest_weights, grf_split_frequencies,
*         grf_merge_forests, grf_plot_tree

clear all
set more off

local errors = 0

* ============================================================
* Setup data and baseline model
* ============================================================
clear
set obs 400
set seed 42
forvalues j = 1/5 {
    gen x`j' = rnormal()
}
gen y = 2*x1 + x2^2 + rnormal()

grf_regression_forest y x1-x5, gen(pred_base) ntrees(200) seed(42)
assert !missing(e(model_id))

* ---- Test 1: forest summary ----
capture noisily {
    grf_forest_summary, all
    assert "`r(forest_type)'" == "regression"
    assert r(N) == 400
    assert r(n_trees) == 200
    assert r(model_id) == e(model_id)
    assert strpos("`r(e_scalars)'", "N") > 0
    assert strpos("`r(e_macros)'", "forest_type") > 0
}
if _rc {
    display as error "FAIL: grf_forest_summary"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_forest_summary"
}

* ---- Test 2: tree summary / get_tree ----
capture noisily {
    grf_tree_summary, tree(1)
    assert r(tree) == 1
    assert r(n_trees) == 200
    assert r(model_id) == e(model_id)

    grf_get_tree, tree(1)
    assert r(tree) == 1
}
if _rc {
    display as error "FAIL: grf_tree_summary / grf_get_tree"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_tree_summary / grf_get_tree"
}

* ---- Test 3: leaf-node proxy assignment ----
capture noisily {
    grf_get_leaf_node, gen(leaf_id) groups(15)
    assert r(groups) == 15
    quietly count if !missing(leaf_id)
    assert r(N) > 0
    drop leaf_id
}
if _rc {
    display as error "FAIL: grf_get_leaf_node"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_get_leaf_node"
}

* ---- Test 4: forest-weights proxy ----
capture noisily {
    grf_get_forest_weights, obs(10) gen(fw_proxy) scale(0.5)
    quietly summarize fw_proxy
    assert r(N) > 0
    assert abs(r(sum) - 1) < 1e-8
    assert r(Var) > 0

    grf_get_forest_weights, obs(10) gen(fw_proxy_x) scale(0.5) ///
        xvars(x1 x2 x3) predweight(0.5) xscale(1.0)
    assert r(uses_xvars) == 1
    assert reldif(r(predweight), 0.5) < 1e-8
    quietly summarize fw_proxy_x
    assert r(N) > 0
    assert abs(r(sum) - 1) < 1e-8
    assert r(Var) > 0

    gen double fw_diff = fw_proxy - fw_proxy_x
    quietly summarize fw_diff
    assert r(sd) > 0

    drop fw_proxy fw_proxy_x fw_diff
}
if _rc {
    display as error "FAIL: grf_get_forest_weights"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_get_forest_weights"
}

* ---- Test 5: split frequencies proxy ----
capture noisily {
    grf_split_frequencies y x1-x5, ntrees(200) seed(42) maxdepth(3)
    matrix sf = r(split_frequencies)
    assert rowsof(sf) == 1
    assert colsof(sf) == 5
}
if _rc {
    display as error "FAIL: grf_split_frequencies"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_split_frequencies"
}

* ---- Test 6: merge forests ----
capture noisily {
    grf_regression_forest y x1-x5, gen(pred_alt) ntrees(200) seed(99) replace
    grf_merge_forests pred_base pred_alt, gen(pred_merge) weights(0.7 0.3)
    assert reldif(r(w1), 0.7) < 1e-8
    assert reldif(r(w2), 0.3) < 1e-8
    quietly count if !missing(pred_merge)
    assert r(N) > 0
    drop pred_alt pred_merge
}
if _rc {
    display as error "FAIL: grf_merge_forests"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_merge_forests"
}

* ---- Test 7: plot tree proxy ----
capture noisily {
    grf_plot_tree y x1-x5, ntrees(100) seed(42) maxdepth(3) ///
        name(grf_tree_plot_test) replace
    assert "`r(graph_name)'" == "grf_tree_plot_test"
}
if _rc {
    display as error "FAIL: grf_plot_tree"
    local errors = `errors' + 1
}
else {
    display as result "PASS: grf_plot_tree"
}

* ============================================================
* Summary
* ============================================================
display ""
if `errors' > 0 {
    display as error "FAILED: `errors' forest utility test(s) did not pass"
    exit 1
}
else {
    display as result "ALL FOREST UTILITY TESTS PASSED"
}
