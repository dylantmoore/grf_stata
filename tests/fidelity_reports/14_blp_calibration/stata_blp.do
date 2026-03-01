/* ============================================================
   stata_blp.do
   Run BLP, test_calibration, and get_scores tests
   R vs Stata fidelity comparison
   NOTE: Stata lowercases CSV variable names on import
   NOTE: grf_best_linear_projection replaces e() so we must
         re-invoke grf_causal_forest before each BLP call.
         We use yhatinput/whatinput to skip re-training nuisance.
   ============================================================ */

adopath ++ /tmp/grf_stata
set more off
set linesize 200
cd /tmp/grf_stata/tests/fidelity_reports/14_blp_calibration

/* ============================================================
   Open results collection file
   ============================================================ */
capture file close fres
file open fres using "stata_results.csv", write replace
file write fres "test,param,value" _n

/* ============================================================
   Prepare main dataset (saved as dta for fast reloading)
   ============================================================ */
import delimited "data.csv", clear varnames(1)
gen _n_obs = _n
save "/tmp/14blp_data.dta", replace

import delimited "predictions.csv", clear varnames(1)
gen _n_obs = _n
save "/tmp/14blp_preds.dta", replace

use "/tmp/14blp_data.dta", clear
merge 1:1 _n_obs using "/tmp/14blp_preds.dta", nogenerate
drop _n_obs
save "/tmp/14blp_combined.dta", replace

/* ============================================================
   Macro: restore causal forest e() context
   We pass pre-computed yhat/what from R so training is fast
   (the forest still runs but nuisance estimates are fixed).
   ============================================================ */

/* ---- Test 1: BLP HC3 (default), all covariates ---- */
di _n "=== Test 1: BLP HC3 default, all covariates ==="
use "/tmp/14blp_combined.dta", clear
grf_causal_forest y w x1 x2 x3 x4 x5, ///
    generate(_tau1) yhatinput(yhat) whatinput(what) ntrees(2000) seed(42) replace
grf_best_linear_projection x1 x2 x3 x4 x5, vcovtype(HC3)
file write fres "blp_hc3_all,coef__cons," (string(_b[_cons], "%20.12f")) _n
file write fres "blp_hc3_all,se__cons," (string(_se[_cons], "%20.12f")) _n
file write fres "blp_hc3_all,coef_x1," (string(_b[x1], "%20.12f")) _n
file write fres "blp_hc3_all,se_x1," (string(_se[x1], "%20.12f")) _n
file write fres "blp_hc3_all,coef_x2," (string(_b[x2], "%20.12f")) _n
file write fres "blp_hc3_all,se_x2," (string(_se[x2], "%20.12f")) _n
file write fres "blp_hc3_all,coef_x3," (string(_b[x3], "%20.12f")) _n
file write fres "blp_hc3_all,se_x3," (string(_se[x3], "%20.12f")) _n
file write fres "blp_hc3_all,coef_x4," (string(_b[x4], "%20.12f")) _n
file write fres "blp_hc3_all,se_x4," (string(_se[x4], "%20.12f")) _n
file write fres "blp_hc3_all,coef_x5," (string(_b[x5], "%20.12f")) _n
file write fres "blp_hc3_all,se_x5," (string(_se[x5], "%20.12f")) _n

/* ---- Test 2: BLP HC0 ---- */
di _n "=== Test 2: BLP HC0 ==="
use "/tmp/14blp_combined.dta", clear
grf_causal_forest y w x1 x2 x3 x4 x5, ///
    generate(_tau2) yhatinput(yhat) whatinput(what) ntrees(2000) seed(42) replace
grf_best_linear_projection x1 x2 x3 x4 x5, vcovtype(HC0)
file write fres "blp_hc0,coef__cons," (string(_b[_cons], "%20.12f")) _n
file write fres "blp_hc0,se__cons," (string(_se[_cons], "%20.12f")) _n
file write fres "blp_hc0,coef_x1," (string(_b[x1], "%20.12f")) _n
file write fres "blp_hc0,se_x1," (string(_se[x1], "%20.12f")) _n
file write fres "blp_hc0,coef_x2," (string(_b[x2], "%20.12f")) _n
file write fres "blp_hc0,se_x2," (string(_se[x2], "%20.12f")) _n
file write fres "blp_hc0,coef_x3," (string(_b[x3], "%20.12f")) _n
file write fres "blp_hc0,se_x3," (string(_se[x3], "%20.12f")) _n
file write fres "blp_hc0,coef_x4," (string(_b[x4], "%20.12f")) _n
file write fres "blp_hc0,se_x4," (string(_se[x4], "%20.12f")) _n
file write fres "blp_hc0,coef_x5," (string(_b[x5], "%20.12f")) _n
file write fres "blp_hc0,se_x5," (string(_se[x5], "%20.12f")) _n

/* ---- Test 3: BLP HC1 ---- */
di _n "=== Test 3: BLP HC1 ==="
use "/tmp/14blp_combined.dta", clear
grf_causal_forest y w x1 x2 x3 x4 x5, ///
    generate(_tau3) yhatinput(yhat) whatinput(what) ntrees(2000) seed(42) replace
grf_best_linear_projection x1 x2 x3 x4 x5, vcovtype(HC1)
file write fres "blp_hc1,coef__cons," (string(_b[_cons], "%20.12f")) _n
file write fres "blp_hc1,se__cons," (string(_se[_cons], "%20.12f")) _n
file write fres "blp_hc1,coef_x1," (string(_b[x1], "%20.12f")) _n
file write fres "blp_hc1,se_x1," (string(_se[x1], "%20.12f")) _n
file write fres "blp_hc1,coef_x2," (string(_b[x2], "%20.12f")) _n
file write fres "blp_hc1,se_x2," (string(_se[x2], "%20.12f")) _n
file write fres "blp_hc1,coef_x3," (string(_b[x3], "%20.12f")) _n
file write fres "blp_hc1,se_x3," (string(_se[x3], "%20.12f")) _n
file write fres "blp_hc1,coef_x4," (string(_b[x4], "%20.12f")) _n
file write fres "blp_hc1,se_x4," (string(_se[x4], "%20.12f")) _n
file write fres "blp_hc1,coef_x5," (string(_b[x5], "%20.12f")) _n
file write fres "blp_hc1,se_x5," (string(_se[x5], "%20.12f")) _n

/* ---- Test 4: BLP HC2 ---- */
di _n "=== Test 4: BLP HC2 ==="
use "/tmp/14blp_combined.dta", clear
grf_causal_forest y w x1 x2 x3 x4 x5, ///
    generate(_tau4) yhatinput(yhat) whatinput(what) ntrees(2000) seed(42) replace
grf_best_linear_projection x1 x2 x3 x4 x5, vcovtype(HC2)
file write fres "blp_hc2,coef__cons," (string(_b[_cons], "%20.12f")) _n
file write fres "blp_hc2,se__cons," (string(_se[_cons], "%20.12f")) _n
file write fres "blp_hc2,coef_x1," (string(_b[x1], "%20.12f")) _n
file write fres "blp_hc2,se_x1," (string(_se[x1], "%20.12f")) _n
file write fres "blp_hc2,coef_x2," (string(_b[x2], "%20.12f")) _n
file write fres "blp_hc2,se_x2," (string(_se[x2], "%20.12f")) _n
file write fres "blp_hc2,coef_x3," (string(_b[x3], "%20.12f")) _n
file write fres "blp_hc2,se_x3," (string(_se[x3], "%20.12f")) _n
file write fres "blp_hc2,coef_x4," (string(_b[x4], "%20.12f")) _n
file write fres "blp_hc2,se_x4," (string(_se[x4], "%20.12f")) _n
file write fres "blp_hc2,coef_x5," (string(_b[x5], "%20.12f")) _n
file write fres "blp_hc2,se_x5," (string(_se[x5], "%20.12f")) _n

/* ---- Test 5: BLP subset x1, x2 ---- */
di _n "=== Test 5: BLP subset x1 x2 ==="
use "/tmp/14blp_combined.dta", clear
grf_causal_forest y w x1 x2 x3 x4 x5, ///
    generate(_tau5) yhatinput(yhat) whatinput(what) ntrees(2000) seed(42) replace
grf_best_linear_projection x1 x2, vcovtype(HC3)
file write fres "blp_subset_x1x2,coef__cons," (string(_b[_cons], "%20.12f")) _n
file write fres "blp_subset_x1x2,se__cons," (string(_se[_cons], "%20.12f")) _n
file write fres "blp_subset_x1x2,coef_x1," (string(_b[x1], "%20.12f")) _n
file write fres "blp_subset_x1x2,se_x1," (string(_se[x1], "%20.12f")) _n
file write fres "blp_subset_x1x2,coef_x2," (string(_b[x2], "%20.12f")) _n
file write fres "blp_subset_x1x2,se_x2," (string(_se[x2], "%20.12f")) _n

/* ---- Test 6: BLP target.sample=treated (Stata-only) ---- */
di _n "=== Test 6: BLP target.sample=treated (Stata-only) ==="
use "/tmp/14blp_combined.dta", clear
grf_causal_forest y w x1 x2 x3 x4 x5, ///
    generate(_tau6) yhatinput(yhat) whatinput(what) ntrees(2000) seed(42) replace
grf_best_linear_projection x1 x2 x3 x4 x5, vcovtype(HC3) targetsample(treated)
file write fres "blp_treated,coef__cons," (string(_b[_cons], "%20.12f")) _n
file write fres "blp_treated,se__cons," (string(_se[_cons], "%20.12f")) _n
file write fres "blp_treated,coef_x1," (string(_b[x1], "%20.12f")) _n
file write fres "blp_treated,se_x1," (string(_se[x1], "%20.12f")) _n
file write fres "blp_treated,coef_x2," (string(_b[x2], "%20.12f")) _n
file write fres "blp_treated,se_x2," (string(_se[x2], "%20.12f")) _n
file write fres "blp_treated,coef_x3," (string(_b[x3], "%20.12f")) _n
file write fres "blp_treated,se_x3," (string(_se[x3], "%20.12f")) _n
file write fres "blp_treated,coef_x4," (string(_b[x4], "%20.12f")) _n
file write fres "blp_treated,se_x4," (string(_se[x4], "%20.12f")) _n
file write fres "blp_treated,coef_x5," (string(_b[x5], "%20.12f")) _n
file write fres "blp_treated,se_x5," (string(_se[x5], "%20.12f")) _n

/* ---- Test 7: BLP target.sample=control (Stata-only) ---- */
di _n "=== Test 7: BLP target.sample=control (Stata-only) ==="
use "/tmp/14blp_combined.dta", clear
grf_causal_forest y w x1 x2 x3 x4 x5, ///
    generate(_tau7) yhatinput(yhat) whatinput(what) ntrees(2000) seed(42) replace
grf_best_linear_projection x1 x2 x3 x4 x5, vcovtype(HC3) targetsample(control)
file write fres "blp_control,coef__cons," (string(_b[_cons], "%20.12f")) _n
file write fres "blp_control,se__cons," (string(_se[_cons], "%20.12f")) _n
file write fres "blp_control,coef_x1," (string(_b[x1], "%20.12f")) _n
file write fres "blp_control,se_x1," (string(_se[x1], "%20.12f")) _n
file write fres "blp_control,coef_x2," (string(_b[x2], "%20.12f")) _n
file write fres "blp_control,se_x2," (string(_se[x2], "%20.12f")) _n
file write fres "blp_control,coef_x3," (string(_b[x3], "%20.12f")) _n
file write fres "blp_control,se_x3," (string(_se[x3], "%20.12f")) _n
file write fres "blp_control,coef_x4," (string(_b[x4], "%20.12f")) _n
file write fres "blp_control,se_x4," (string(_se[x4], "%20.12f")) _n
file write fres "blp_control,coef_x5," (string(_b[x5], "%20.12f")) _n
file write fres "blp_control,se_x5," (string(_se[x5], "%20.12f")) _n

/* ---- Test 8: BLP target.sample=overlap ---- */
di _n "=== Test 8: BLP target.sample=overlap ==="
use "/tmp/14blp_combined.dta", clear
grf_causal_forest y w x1 x2 x3 x4 x5, ///
    generate(_tau8) yhatinput(yhat) whatinput(what) ntrees(2000) seed(42) replace
grf_best_linear_projection x1 x2 x3 x4 x5, vcovtype(HC3) targetsample(overlap)
file write fres "blp_overlap,coef__cons," (string(_b[_cons], "%20.12f")) _n
file write fres "blp_overlap,se__cons," (string(_se[_cons], "%20.12f")) _n
file write fres "blp_overlap,coef_x1," (string(_b[x1], "%20.12f")) _n
file write fres "blp_overlap,se_x1," (string(_se[x1], "%20.12f")) _n
file write fres "blp_overlap,coef_x2," (string(_b[x2], "%20.12f")) _n
file write fres "blp_overlap,se_x2," (string(_se[x2], "%20.12f")) _n
file write fres "blp_overlap,coef_x3," (string(_b[x3], "%20.12f")) _n
file write fres "blp_overlap,se_x3," (string(_se[x3], "%20.12f")) _n
file write fres "blp_overlap,coef_x4," (string(_b[x4], "%20.12f")) _n
file write fres "blp_overlap,se_x4," (string(_se[x4], "%20.12f")) _n
file write fres "blp_overlap,coef_x5," (string(_b[x5], "%20.12f")) _n
file write fres "blp_overlap,se_x5," (string(_se[x5], "%20.12f")) _n

/* ---- Test 9: BLP with debiasing weights ---- */
di _n "=== Test 9: BLP with debiasing weights ==="
use "/tmp/14blp_combined.dta", clear
grf_causal_forest y w x1 x2 x3 x4 x5, ///
    generate(_tau9) yhatinput(yhat) whatinput(what) ntrees(2000) seed(42) replace
grf_best_linear_projection x1 x2 x3 x4 x5, vcovtype(HC3) debiasingweights(dbw)
file write fres "blp_dbw,coef__cons," (string(_b[_cons], "%20.12f")) _n
file write fres "blp_dbw,se__cons," (string(_se[_cons], "%20.12f")) _n
file write fres "blp_dbw,coef_x1," (string(_b[x1], "%20.12f")) _n
file write fres "blp_dbw,se_x1," (string(_se[x1], "%20.12f")) _n
file write fres "blp_dbw,coef_x2," (string(_b[x2], "%20.12f")) _n
file write fres "blp_dbw,se_x2," (string(_se[x2], "%20.12f")) _n
file write fres "blp_dbw,coef_x3," (string(_b[x3], "%20.12f")) _n
file write fres "blp_dbw,se_x3," (string(_se[x3], "%20.12f")) _n
file write fres "blp_dbw,coef_x4," (string(_b[x4], "%20.12f")) _n
file write fres "blp_dbw,se_x4," (string(_se[x4], "%20.12f")) _n
file write fres "blp_dbw,coef_x5," (string(_b[x5], "%20.12f")) _n
file write fres "blp_dbw,se_x5," (string(_se[x5], "%20.12f")) _n

/* ---- Test 10: BLP with clusters ---- */
di _n "=== Test 10: BLP with clusters ==="
import delimited "data.csv", clear varnames(1)
gen _n_obs = _n
save "/tmp/14blp_data_cl.dta", replace
import delimited "predictions_cl.csv", clear varnames(1)
gen _n_obs = _n
save "/tmp/14blp_preds_cl.dta", replace
use "/tmp/14blp_data_cl.dta", clear
merge 1:1 _n_obs using "/tmp/14blp_preds_cl.dta", nogenerate
drop _n_obs
save "/tmp/14blp_combined_cl.dta", replace

use "/tmp/14blp_combined_cl.dta", clear
grf_causal_forest y w x1 x2 x3 x4 x5, ///
    generate(_tau10) yhatinput(yhat) whatinput(what) ///
    cluster(clusters) ntrees(2000) seed(42) replace
grf_best_linear_projection x1 x2 x3 x4 x5, vcovtype(HC3)
file write fres "blp_clusters,coef__cons," (string(_b[_cons], "%20.12f")) _n
file write fres "blp_clusters,se__cons," (string(_se[_cons], "%20.12f")) _n
file write fres "blp_clusters,coef_x1," (string(_b[x1], "%20.12f")) _n
file write fres "blp_clusters,se_x1," (string(_se[x1], "%20.12f")) _n
file write fres "blp_clusters,coef_x2," (string(_b[x2], "%20.12f")) _n
file write fres "blp_clusters,se_x2," (string(_se[x2], "%20.12f")) _n
file write fres "blp_clusters,coef_x3," (string(_b[x3], "%20.12f")) _n
file write fres "blp_clusters,se_x3," (string(_se[x3], "%20.12f")) _n
file write fres "blp_clusters,coef_x4," (string(_b[x4], "%20.12f")) _n
file write fres "blp_clusters,se_x4," (string(_se[x4], "%20.12f")) _n
file write fres "blp_clusters,coef_x5," (string(_b[x5], "%20.12f")) _n
file write fres "blp_clusters,se_x5," (string(_se[x5], "%20.12f")) _n

/* ============================================================
   Test Calibration Tests
   Note: grf_test_calibration uses e() from causal forest
         but it issues a regress internally (which clears e()).
         Subsequent BLP calls need re-fitting. Since calibration
         uses r() not e(), it is safe to call in sequence.
   ============================================================ */

/* ---- Test 13: Basic calibration ---- */
di _n "=== Test 13: Basic calibration ==="
use "/tmp/14blp_combined.dta", clear
grf_causal_forest y w x1 x2 x3 x4 x5, ///
    generate(_tau13) yhatinput(yhat) whatinput(what) ntrees(2000) seed(42) replace
grf_test_calibration
file write fres "calibration_basic,b_mean," (string(r(b_mean), "%20.12f")) _n
file write fres "calibration_basic,se_mean," (string(r(se_mean), "%20.12f")) _n
file write fres "calibration_basic,t_mean," (string(r(t_mean), "%20.12f")) _n
file write fres "calibration_basic,p_mean," (string(r(p_mean), "%20.12f")) _n
file write fres "calibration_basic,b_diff," (string(r(b_diff), "%20.12f")) _n
file write fres "calibration_basic,se_diff," (string(r(se_diff), "%20.12f")) _n
file write fres "calibration_basic,t_diff," (string(r(t_diff), "%20.12f")) _n
file write fres "calibration_basic,p_diff," (string(r(p_diff), "%20.12f")) _n

/* ---- Test 14: Strong heterogeneity ---- */
di _n "=== Test 14: Calibration strong heterogeneity ==="
import delimited "data_strong.csv", clear varnames(1)
gen _n_obs = _n
save "/tmp/14blp_data_s.dta", replace
import delimited "predictions_strong.csv", clear varnames(1)
gen _n_obs = _n
save "/tmp/14blp_preds_s.dta", replace
use "/tmp/14blp_data_s.dta", clear
merge 1:1 _n_obs using "/tmp/14blp_preds_s.dta", nogenerate
drop _n_obs

grf_causal_forest y_strong w x1 x2 x3 x4 x5, ///
    generate(_tau_s) yhatinput(yhat) whatinput(what) ntrees(2000) seed(42) replace
grf_test_calibration
file write fres "calibration_strong_het,b_mean," (string(r(b_mean), "%20.12f")) _n
file write fres "calibration_strong_het,se_mean," (string(r(se_mean), "%20.12f")) _n
file write fres "calibration_strong_het,t_mean," (string(r(t_mean), "%20.12f")) _n
file write fres "calibration_strong_het,p_mean," (string(r(p_mean), "%20.12f")) _n
file write fres "calibration_strong_het,b_diff," (string(r(b_diff), "%20.12f")) _n
file write fres "calibration_strong_het,se_diff," (string(r(se_diff), "%20.12f")) _n
file write fres "calibration_strong_het,t_diff," (string(r(t_diff), "%20.12f")) _n
file write fres "calibration_strong_het,p_diff," (string(r(p_diff), "%20.12f")) _n

/* ---- Test 15: No heterogeneity ---- */
di _n "=== Test 15: Calibration no heterogeneity ==="
import delimited "data_const.csv", clear varnames(1)
gen _n_obs = _n
save "/tmp/14blp_data_c.dta", replace
import delimited "predictions_const.csv", clear varnames(1)
gen _n_obs = _n
save "/tmp/14blp_preds_c.dta", replace
use "/tmp/14blp_data_c.dta", clear
merge 1:1 _n_obs using "/tmp/14blp_preds_c.dta", nogenerate
drop _n_obs

grf_causal_forest y_const w x1 x2 x3 x4 x5, ///
    generate(_tau_c) yhatinput(yhat) whatinput(what) ntrees(2000) seed(42) replace
grf_test_calibration
file write fres "calibration_no_het,b_mean," (string(r(b_mean), "%20.12f")) _n
file write fres "calibration_no_het,se_mean," (string(r(se_mean), "%20.12f")) _n
file write fres "calibration_no_het,t_mean," (string(r(t_mean), "%20.12f")) _n
file write fres "calibration_no_het,p_mean," (string(r(p_mean), "%20.12f")) _n
file write fres "calibration_no_het,b_diff," (string(r(b_diff), "%20.12f")) _n
file write fres "calibration_no_het,se_diff," (string(r(se_diff), "%20.12f")) _n
file write fres "calibration_no_het,t_diff," (string(r(t_diff), "%20.12f")) _n
file write fres "calibration_no_het,p_diff," (string(r(p_diff), "%20.12f")) _n

/* ============================================================
   Get Scores Tests
   ============================================================ */
di _n "=== Test 17: Get DR scores ==="
use "/tmp/14blp_combined.dta", clear
grf_causal_forest y w x1 x2 x3 x4 x5, ///
    generate(_tau17) yhatinput(yhat) whatinput(what) ntrees(2000) seed(42) replace
grf_get_scores, gen(dr_scores_stata)

file write fres "dr_scores,mean," (string(r(mean), "%20.12f")) _n
file write fres "dr_scores,sd," (string(r(sd), "%20.12f")) _n
file write fres "dr_scores,min," (string(r(min), "%20.12f")) _n
file write fres "dr_scores,max," (string(r(max), "%20.12f")) _n
file write fres "dr_scores,N," (string(r(N), "%20.0f")) _n

keep dr_scores_stata
export delimited using "stata_dr_scores.csv", replace

file close fres
di _n "Stata do-file completed successfully."
di "Results: stata_results.csv"
di "DR scores: stata_dr_scores.csv"
