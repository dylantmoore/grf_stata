/* ============================================================
   run_all_stata.do — Survival forest fidelity tests
   Runs all 17 test cases and exports comparison CSVs
   ============================================================ */

adopath ++ /tmp/grf_stata

local OUTDIR "/tmp/grf_stata/tests/fidelity_reports/07_survival"

/* ============================================================
   TEST 01: Default — noutput=20, Kaplan-Meier
   ============================================================ */
di ""
di "=== TEST 01: Default (noutput=20, KM) ==="
import delimited "`OUTDIR'/test01_default_data.csv", clear

grf_survival_forest time status x1 x2 x3 x4 x5, ///
    gen(pred) ntrees(500) seed(42) noutput(20) predtype(1) replace

keep pred_s1-pred_s20 r_s1-r_s20
export delimited "`OUTDIR'/test01_stata.csv", replace
di "TEST 01 DONE"

/* ============================================================
   TEST 02: noutput=50
   ============================================================ */
di ""
di "=== TEST 02: noutput=50 ==="
import delimited "`OUTDIR'/test02_noutput50_data.csv", clear

grf_survival_forest time status x1 x2 x3 x4 x5, ///
    gen(pred) ntrees(500) seed(42) noutput(50) predtype(1) replace

keep pred_s1-pred_s50 r_s1-r_s50
export delimited "`OUTDIR'/test02_stata.csv", replace
di "TEST 02 DONE"

/* ============================================================
   TEST 03: noutput=100
   ============================================================ */
di ""
di "=== TEST 03: noutput=100 ==="
import delimited "`OUTDIR'/test03_noutput100_data.csv", clear

grf_survival_forest time status x1 x2 x3 x4 x5, ///
    gen(pred) ntrees(500) seed(42) noutput(100) predtype(1) replace

keep pred_s1-pred_s100 r_s1-r_s100
export delimited "`OUTDIR'/test03_stata.csv", replace
di "TEST 03 DONE"

/* ============================================================
   TEST 04: predtype=0 (Nelson-Aalen)
   ============================================================ */
di ""
di "=== TEST 04: predtype=0 (Nelson-Aalen) ==="
import delimited "`OUTDIR'/test04_predtype0_data.csv", clear

grf_survival_forest time status x1 x2 x3 x4 x5, ///
    gen(pred) ntrees(500) seed(42) noutput(20) predtype(0) replace

keep pred_s1-pred_s20 r_s1-r_s20
export delimited "`OUTDIR'/test04_stata.csv", replace
di "TEST 04 DONE"

/* ============================================================
   TEST 05: predtype=1 (Kaplan-Meier, explicit)
   ============================================================ */
di ""
di "=== TEST 05: predtype=1 (KM explicit) ==="
import delimited "`OUTDIR'/test05_predtype1_data.csv", clear

grf_survival_forest time status x1 x2 x3 x4 x5, ///
    gen(pred) ntrees(500) seed(42) noutput(20) predtype(1) replace

keep pred_s1-pred_s20 r_s1-r_s20
export delimited "`OUTDIR'/test05_stata.csv", replace
di "TEST 05 DONE"

/* ============================================================
   TEST 06: nofastlogrank
   R default is fast.logrank=FALSE; Stata default may differ
   ============================================================ */
di ""
di "=== TEST 06: nofastlogrank ==="
import delimited "`OUTDIR'/test06_nofastlogrank_data.csv", clear

grf_survival_forest time status x1 x2 x3 x4 x5, ///
    gen(pred) ntrees(500) seed(42) noutput(20) predtype(1) nofastlogrank replace

keep pred_s1-pred_s20 r_s1-r_s20
export delimited "`OUTDIR'/test06_stata.csv", replace
di "TEST 06 DONE"

/* ============================================================
   TEST 07: cluster()
   ============================================================ */
di ""
di "=== TEST 07: cluster() ==="
import delimited "`OUTDIR'/test07_cluster_data.csv", clear

grf_survival_forest time status x1 x2 x3 x4 x5, ///
    gen(pred) ntrees(500) seed(42) noutput(20) predtype(1) ///
    cluster(cluster_id) replace

keep pred_s1-pred_s20 r_s1-r_s20
export delimited "`OUTDIR'/test07_stata.csv", replace
di "TEST 07 DONE"

/* ============================================================
   TEST 08: weights()
   ============================================================ */
di ""
di "=== TEST 08: weights() ==="
import delimited "`OUTDIR'/test08_weights_data.csv", clear

grf_survival_forest time status x1 x2 x3 x4 x5, ///
    gen(pred) ntrees(500) seed(42) noutput(20) predtype(1) ///
    weights(wt) replace

keep pred_s1-pred_s20 r_s1-r_s20
export delimited "`OUTDIR'/test08_stata.csv", replace
di "TEST 08 DONE"

/* ============================================================
   TEST 09: nohonesty
   ============================================================ */
di ""
di "=== TEST 09: nohonesty ==="
import delimited "`OUTDIR'/test09_nohonesty_data.csv", clear

grf_survival_forest time status x1 x2 x3 x4 x5, ///
    gen(pred) ntrees(500) seed(42) noutput(20) predtype(1) ///
    nohonesty replace

keep pred_s1-pred_s20 r_s1-r_s20
export delimited "`OUTDIR'/test09_stata.csv", replace
di "TEST 09 DONE"

/* ============================================================
   TEST 10: mtry=2
   ============================================================ */
di ""
di "=== TEST 10: mtry=2 ==="
import delimited "`OUTDIR'/test10_mtry2_data.csv", clear

grf_survival_forest time status x1 x2 x3 x4 x5, ///
    gen(pred) ntrees(500) seed(42) noutput(20) predtype(1) ///
    mtry(2) replace

keep pred_s1-pred_s20 r_s1-r_s20
export delimited "`OUTDIR'/test10_stata.csv", replace
di "TEST 10 DONE"

/* ============================================================
   TEST 11: minnodesize=20
   ============================================================ */
di ""
di "=== TEST 11: minnodesize=20 ==="
import delimited "`OUTDIR'/test11_minnodesize20_data.csv", clear

grf_survival_forest time status x1 x2 x3 x4 x5, ///
    gen(pred) ntrees(500) seed(42) noutput(20) predtype(1) ///
    minnodesize(20) replace

keep pred_s1-pred_s20 r_s1-r_s20
export delimited "`OUTDIR'/test11_stata.csv", replace
di "TEST 11 DONE"

/* ============================================================
   TEST 12: Heavy censoring (high event rate, rate=0.05 censoring)
   ============================================================ */
di ""
di "=== TEST 12: High event rate (near ~94%) ==="
import delimited "`OUTDIR'/test12_heavy_censoring_data.csv", clear

grf_survival_forest time status x1 x2 x3 x4 x5, ///
    gen(pred) ntrees(500) seed(42) noutput(20) predtype(1) replace

keep pred_s1-pred_s20 r_s1-r_s20
export delimited "`OUTDIR'/test12_stata.csv", replace
di "TEST 12 DONE"

/* ============================================================
   TEST 13: Light censoring (event rate ~34%)
   ============================================================ */
di ""
di "=== TEST 13: Lower event rate (~34%) ==="
import delimited "`OUTDIR'/test13_light_censoring_data.csv", clear

grf_survival_forest time status x1 x2 x3 x4 x5, ///
    gen(pred) ntrees(500) seed(42) noutput(20) predtype(1) replace

keep pred_s1-pred_s20 r_s1-r_s20
export delimited "`OUTDIR'/test13_stata.csv", replace
di "TEST 13 DONE"

/* ============================================================
   TEST 14: Expected survival
   Fit survival forest and compute expected survival
   ============================================================ */
di ""
di "=== TEST 14: Expected survival ==="
import delimited "`OUTDIR'/test14_expected_survival_data.csv", clear

/* Fit with all available columns to match R's full integration */
/* We use noutput=387 but R only uses the first K that exist */
/* To match R which integrates over all 387 failure times, use high noutput */
grf_survival_forest time status x1 x2 x3 x4 x5, ///
    gen(pred) ntrees(500) seed(42) noutput(387) predtype(1) replace

grf_expected_survival, gen(stata_esurv) predictions(pred)

keep stata_esurv r_esurv
export delimited "`OUTDIR'/test14_stata.csv", replace
di "TEST 14 DONE"

/* ============================================================
   TEST 15: Expected survival consistency
   E[T|X] should be negatively correlated with X1
   ============================================================ */
di ""
di "=== TEST 15: Expected survival consistency ==="
import delimited "`OUTDIR'/test15_esurv_consistency_data.csv", clear

grf_survival_forest time status x1 x2 x3 x4 x5, ///
    gen(pred) ntrees(500) seed(42) noutput(387) predtype(1) replace

grf_expected_survival, gen(stata_esurv) predictions(pred)

/* Check correlation with x1 (should be negative) */
correlate stata_esurv x1
di "Correlation (stata_esurv, x1) = " r(rho) " (should be negative)"

keep stata_esurv r_esurv x1
export delimited "`OUTDIR'/test15_stata.csv", replace
di "TEST 15 DONE"

/* ============================================================
   TEST 16: numfailures(50)
   Limit to 50 unique failure times
   ============================================================ */
di ""
di "=== TEST 16: numfailures(50) ==="
import delimited "`OUTDIR'/test16_numfailures50_data.csv", clear

grf_survival_forest time status x1 x2 x3 x4 x5, ///
    gen(pred) ntrees(500) seed(42) noutput(20) predtype(1) ///
    numfailures(50) replace

keep pred_s1-pred_s20 r_s1-r_s20
export delimited "`OUTDIR'/test16_stata.csv", replace
di "TEST 16 DONE"

/* ============================================================
   TEST 17: Combined noutput=50 + predtype=0 (Nelson-Aalen) + cluster
   ============================================================ */
di ""
di "=== TEST 17: Combined (noutput=50 + NA + cluster) ==="
import delimited "`OUTDIR'/test17_combined_data.csv", clear

grf_survival_forest time status x1 x2 x3 x4 x5, ///
    gen(pred) ntrees(500) seed(42) noutput(50) predtype(0) ///
    cluster(cluster_id) replace

keep pred_s1-pred_s50 r_s1-r_s50
export delimited "`OUTDIR'/test17_stata.csv", replace
di "TEST 17 DONE"

di ""
di "=== ALL STATA TESTS COMPLETE ==="
