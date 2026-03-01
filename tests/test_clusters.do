* test_clusters.do -- Tests for cluster-robust estimation
*
* Tests that the cluster() option enables cluster-aware splitting
* and standard errors in grf forests.
*
* Run: do tests/test_clusters.do

clear all
set more off

local errors = 0

* ============================================================
* 1. Causal forest with clustering succeeds
* ============================================================
display as text ""
display as text "=== Cluster Test 1: Causal forest with cluster() ==="

capture noisily {
    clear
    set obs 500
    set seed 42

    * Generate 50 clusters with 10 obs each
    gen cluster_id = ceil(_n / 10)
    gen x1 = rnormal()
    gen x2 = rnormal()

    * Cluster-level treatment assignment
    gen w = (cluster_id <= 25)

    * Cluster-correlated outcome
    bysort cluster_id: gen double cluster_effect = rnormal() if _n == 1
    bysort cluster_id (cluster_effect): replace cluster_effect = cluster_effect[1]
    gen y = 2*x1 + x1*w + cluster_effect + 0.5*rnormal()

    * Fit with clustering
    grf_causal_forest y w x1 x2, gen(cate_cl) ntrees(500) seed(42) cluster(cluster_id)
    assert e(N) == 500
}
if _rc {
    display as error "FAIL: causal forest with cluster()"
    local errors = `errors' + 1
}
else {
    display as result "PASS: causal forest with cluster()"
}

* ============================================================
* 2. e(cluster_var) stored correctly
* ============================================================
display as text ""
display as text "=== Cluster Test 2: e(cluster_var) stored ==="

capture noisily {
    assert "`e(cluster_var)'" == "cluster_id"
}
if _rc {
    display as error "FAIL: e(cluster_var) not stored"
    local errors = `errors' + 1
}
else {
    display as result "PASS: e(cluster_var) = `e(cluster_var)'"
}

* ============================================================
* 3. ATE SEs differ with and without clustering
* ============================================================
display as text ""
display as text "=== Cluster Test 3: Clustered vs unclustered ATE SEs ==="

capture noisily {
    * ATE with clustering
    grf_ate
    local se_cluster = r(se)

    * Re-fit without clustering
    grf_causal_forest y w x1 x2, gen(cate_nocl) ntrees(500) seed(42) replace
    grf_ate
    local se_nocluster = r(se)

    * SEs should differ when data has cluster structure
    assert `se_cluster' != `se_nocluster'
}
if _rc {
    display as error "FAIL: clustered vs unclustered ATE SEs"
    local errors = `errors' + 1
}
else {
    display as result "PASS: clustered SE=`se_cluster' vs unclustered SE=`se_nocluster'"
}

* ============================================================
* 4. Regression forest with clustering (smoke test)
* ============================================================
display as text ""
display as text "=== Cluster Test 4: Regression forest with cluster() ==="

capture noisily {
    grf_regression_forest y x1 x2, gen(pred_cl) ntrees(500) seed(42) ///
        cluster(cluster_id) replace
    assert e(N) == 500
    assert "`e(cluster_var)'" == "cluster_id"
}
if _rc {
    display as error "FAIL: regression forest with cluster()"
    local errors = `errors' + 1
}
else {
    display as result "PASS: regression forest with cluster()"
}

* ============================================================
* 5. Error: nonexistent cluster variable
* ============================================================
display as text ""
display as text "=== Cluster Test 5: Error on nonexistent cluster var ==="

capture {
    grf_regression_forest y x1 x2, gen(pred_badcl) ntrees(500) seed(42) ///
        cluster(nonexistent_var) replace
}
if _rc {
    display as result "PASS: error on nonexistent cluster variable"
}
else {
    display as error "FAIL: should have errored on nonexistent cluster var"
    local errors = `errors' + 1
}

* ============================================================
* Summary
* ============================================================
display ""
if `errors' > 0 {
    display as error "FAILED: `errors' cluster test(s) did not pass"
    exit 1
}
else {
    display as result "ALL CLUSTER TESTS PASSED"
}
