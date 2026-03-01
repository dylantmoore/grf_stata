* ATE Fidelity Tests - Stata
* All 19 tests for grf_ate
* Run after run_r.R has generated CSV data files

adopath ++ /tmp/grf_stata

local outdir "/tmp/grf_stata/tests/fidelity_reports/03_ate"

* Initialize results file
capture erase "`outdir'/stata_results.csv"
file open fh using "`outdir'/stata_results.csv", write replace
file write fh "test_id,test_name,ate,se,errored,note" _n

********************************************************************************
* TEST 01: ATE (all) — default AIPW
********************************************************************************
di as text "Running Test 01: ATE (all)..."
import delimited "`outdir'/test01_data.csv", clear

grf_causal_forest y w x1 x2 x3 x4 x5, gen(tau) ntrees(2000) seed(42)
grf_ate

local ate01 = r(ate)
local se01  = r(se)
file write fh `"1,ATE_all,`ate01',`se01',0,"' _n
di "  Stata Test 01: ATE=" `ate01' " SE=" `se01'

********************************************************************************
* TEST 02: ATT (treated)
********************************************************************************
di as text "Running Test 02: ATT (treated)..."
grf_ate, targetsample(treated)

local ate02 = r(ate)
local se02  = r(se)
file write fh `"2,ATT_treated,`ate02',`se02',0,"' _n
di "  Stata Test 02: ATE=" `ate02' " SE=" `se02'

********************************************************************************
* TEST 03: ATC (control)
********************************************************************************
di as text "Running Test 03: ATC (control)..."
grf_ate, targetsample(control)

local ate03 = r(ate)
local se03  = r(se)
file write fh `"3,ATC_control,`ate03',`se03',0,"' _n
di "  Stata Test 03: ATE=" `ate03' " SE=" `se03'

********************************************************************************
* TEST 04: Overlap weights
********************************************************************************
di as text "Running Test 04: Overlap..."
grf_ate, targetsample(overlap)

local ate04 = r(ate)
local se04  = r(se)
file write fh `"4,ATE_overlap,`ate04',`se04',0,"' _n
di "  Stata Test 04: ATE=" `ate04' " SE=" `se04'

********************************************************************************
* TEST 05: Debiasing weights
********************************************************************************
di as text "Running Test 05: Debiasing weights..."
import delimited "`outdir'/test05_data.csv", clear

grf_causal_forest y w x1 x2 x3 x4 x5, gen(tau) ntrees(2000) seed(42)
grf_ate, debiasingweights(dbw)

local ate05 = r(ate)
local se05  = r(se)
file write fh `"5,ATE_debias,`ate05',`se05',0,"' _n
di "  Stata Test 05: ATE=" `ate05' " SE=" `se05'

********************************************************************************
* TEST 06: Large sample n=2000
********************************************************************************
di as text "Running Test 06: Large n=2000..."
import delimited "`outdir'/test06_data.csv", clear

grf_causal_forest y w x1 x2 x3 x4 x5, gen(tau) ntrees(2000) seed(42)
grf_ate

local ate06 = r(ate)
local se06  = r(se)
file write fh `"6,ATE_large_n,`ate06',`se06',0,"' _n
di "  Stata Test 06: ATE=" `ate06' " SE=" `se06'

********************************************************************************
* TEST 07: TMLE ATE
********************************************************************************
di as text "Running Test 07: TMLE ATE..."
import delimited "`outdir'/test01_data.csv", clear

grf_causal_forest y w x1 x2 x3 x4 x5, gen(tau) ntrees(2000) seed(42)
grf_ate, method(TMLE)

local ate07 = r(ate)
local se07  = r(se)
file write fh `"7,TMLE_ATE,`ate07',`se07',0,"' _n
di "  Stata Test 07: ATE=" `ate07' " SE=" `se07'

********************************************************************************
* TEST 08: TMLE ATT
********************************************************************************
di as text "Running Test 08: TMLE ATT..."
grf_ate, method(TMLE) targetsample(treated)

local ate08 = r(ate)
local se08  = r(se)
file write fh `"8,TMLE_ATT,`ate08',`se08',0,"' _n
di "  Stata Test 08: ATE=" `ate08' " SE=" `se08'

********************************************************************************
* TEST 09: TMLE ATC
********************************************************************************
di as text "Running Test 09: TMLE ATC..."
grf_ate, method(TMLE) targetsample(control)

local ate09 = r(ate)
local se09  = r(se)
file write fh `"9,TMLE_ATC,`ate09',`se09',0,"' _n
di "  Stata Test 09: ATE=" `ate09' " SE=" `se09'

********************************************************************************
* TEST 10: TMLE + overlap → redirect to AIPW
********************************************************************************
di as text "Running Test 10: TMLE + overlap redirect..."
grf_ate, method(TMLE) targetsample(overlap)

local ate10 = r(ate)
local se10  = r(se)
file write fh `"10,TMLE_overlap_redirect,`ate10',`se10',0,TMLE_redirected_to_AIPW"' _n
di "  Stata Test 10: ATE=" `ate10' " SE=" `se10'

********************************************************************************
* TEST 11: TMLE + debiasing.weights → Stata errors rc=198
* Note: R (grf 2.5) silently ignores debiasing.weights with TMLE
********************************************************************************
di as text "Running Test 11: TMLE + debiasing.weights error..."
import delimited "`outdir'/test05_data.csv", clear
grf_causal_forest y w x1 x2 x3 x4 x5, gen(tau) ntrees(2000) seed(42)

capture grf_ate, method(TMLE) debiasingweights(dbw)
local rc11 = _rc
if `rc11' == 198 {
    file write fh `"11,TMLE_debias_error,.,.,1,Stata_rc198_R_silent_ignore"' _n
    di "  Stata Test 11: Correctly errored with rc=198 (R silently ignores)"
}
else if `rc11' == 0 {
    local ate11 = r(ate)
    local se11  = r(se)
    file write fh `"11,TMLE_debias_error,`ate11',`se11',0,no_error_unexpected"' _n
    di "  Stata Test 11: No error raised (unexpected)"
}
else {
    file write fh `"11,TMLE_debias_error,.,.,1,errored_rc`rc11'"' _n
    di "  Stata Test 11: Errored with rc=" `rc11'
}

********************************************************************************
* TEST 12: Constant treatment effect tau=2
********************************************************************************
di as text "Running Test 12: Constant tau=2..."
import delimited "`outdir'/test12_data.csv", clear

grf_causal_forest y w x1 x2 x3 x4 x5, gen(tau) ntrees(2000) seed(42)
grf_ate

local ate12 = r(ate)
local se12  = r(se)
file write fh `"12,const_tau2,`ate12',`se12',0,"' _n
di "  Stata Test 12: ATE=" `ate12' " SE=" `se12'

********************************************************************************
* TEST 13: Heterogeneous tau = 3*X1
********************************************************************************
di as text "Running Test 13: Heterogeneous tau=3*X1..."
import delimited "`outdir'/test13_data.csv", clear

grf_causal_forest y w x1 x2 x3 x4 x5, gen(tau) ntrees(2000) seed(42)
grf_ate

local ate13 = r(ate)
local se13  = r(se)
file write fh `"13,hetero_tau,`ate13',`se13',0,"' _n
di "  Stata Test 13: ATE=" `ate13' " SE=" `se13'

********************************************************************************
* TEST 14: Unbalanced treatment (80% treated)
********************************************************************************
di as text "Running Test 14: Unbalanced 80% treated..."
import delimited "`outdir'/test14_data.csv", clear

grf_causal_forest y w x1 x2 x3 x4 x5, gen(tau) ntrees(2000) seed(42)
grf_ate

local ate14 = r(ate)
local se14  = r(se)
file write fh `"14,unbalanced_80,`ate14',`se14',0,"' _n
di "  Stata Test 14: ATE=" `ate14' " SE=" `se14'

********************************************************************************
* TEST 15: With clusters (20 clusters)
********************************************************************************
di as text "Running Test 15: With clusters..."
import delimited "`outdir'/test15_data.csv", clear

grf_causal_forest y w x1 x2 x3 x4 x5, gen(tau) ntrees(2000) seed(42) ///
    cluster(cluster_id)
grf_ate

local ate15 = r(ate)
local se15  = r(se)
file write fh `"15,cluster_ATE,`ate15',`se15',0,"' _n
di "  Stata Test 15: ATE=" `ate15' " SE=" `se15'

********************************************************************************
* TEST 16: Clusters + TMLE → should error (not implemented)
********************************************************************************
di as text "Running Test 16: Clusters + TMLE..."
capture grf_ate, method(TMLE)
local rc16 = _rc
if `rc16' != 0 {
    file write fh `"16,cluster_TMLE,.,.,1,correctly_errored_rc`rc16'_not_implemented"' _n
    di "  Stata Test 16: Errored rc=" `rc16' " (matches R behavior)"
}
else {
    local ate16 = r(ate)
    local se16  = r(se)
    file write fh `"16,cluster_TMLE,`ate16',`se16',0,unexpected_no_error"' _n
    di "  Stata Test 16: No error (unexpected)"
}

********************************************************************************
* TEST 17: Clustered + equalize cluster weights
* R: equalize.cluster.weights not in grf 2.5, falls back to plain overlap ATE
* Stata: grf_ate, targetsample(overlap) — same fallback
********************************************************************************
di as text "Running Test 17: Clusters + equalize weights..."
* Try with equalizeclusterweights option (may not exist in Stata)
capture grf_ate, targetsample(overlap) equalizeclusterweights
local rc17 = _rc
if `rc17' == 0 {
    local ate17 = r(ate)
    local se17  = r(se)
    file write fh `"17,cluster_equalize,`ate17',`se17',0,equalizeclusterweights_supported"' _n
    di "  Stata Test 17: ATE=" `ate17' " SE=" `se17'
}
else {
    * Fall back to plain overlap (matching R fallback)
    grf_ate, targetsample(overlap)
    local ate17 = r(ate)
    local se17  = r(se)
    file write fh `"17,cluster_equalize,`ate17',`se17',0,equalize_not_supported_used_overlap"' _n
    di "  Stata Test 17: equalizeclusterweights not supported; used overlap. ATE=" `ate17' " SE=" `se17'
}

********************************************************************************
* TEST 18: Large treatment effect tau=10
********************************************************************************
di as text "Running Test 18: Large tau=10..."
import delimited "`outdir'/test18_data.csv", clear

grf_causal_forest y w x1 x2 x3 x4 x5, gen(tau) ntrees(2000) seed(42)
grf_ate

local ate18 = r(ate)
local se18  = r(se)
file write fh `"18,large_tau10,`ate18',`se18',0,"' _n
di "  Stata Test 18: ATE=" `ate18' " SE=" `se18'

********************************************************************************
* TEST 19: Near-zero treatment effect tau=0.01
********************************************************************************
di as text "Running Test 19: Near-zero tau=0.01..."
import delimited "`outdir'/test19_data.csv", clear

grf_causal_forest y w x1 x2 x3 x4 x5, gen(tau) ntrees(2000) seed(42)
grf_ate

local ate19 = r(ate)
local se19  = r(se)
file write fh `"19,nearzero_tau,`ate19',`se19',0,"' _n
di "  Stata Test 19: ATE=" `ate19' " SE=" `se19'

********************************************************************************
* Close file and finish
********************************************************************************
file close fh
di as text "Stata results saved to stata_results.csv"
di as text "All tests complete."
