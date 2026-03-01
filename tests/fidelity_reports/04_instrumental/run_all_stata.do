/* ===================================================================
   Instrumental Forest Fidelity Tests -- Stata side
   Runs all 18 tests, reads R-generated CSVs, runs grf_instrumental_forest,
   saves correlations to a results CSV.
   =================================================================== */

adopath ++ /tmp/grf_stata
local outdir "/tmp/grf_stata/tests/fidelity_reports/04_instrumental"

/* Output CSV for results */
tempname fh
file open `fh' using "`outdir'/stata_results.csv", write replace
file write `fh' "test_id,test_name,n,r_mean,r_sd,stata_mean,stata_sd,correlation,status,notes" _n

/* Helper: compute correlation, write row */
cap program drop write_result
program define write_result
    args fh test_id test_name corr status notes
    quietly summarize r_pred
    local r_mean = r(mean)
    local r_sd   = r(sd)
    quietly summarize stata_pred
    local s_mean = r(mean)
    local s_sd   = r(sd)
    local n_obs  = r(N)
    file write `fh' "`test_id',`test_name',`n_obs'," ///
        "`r_mean',`r_sd',`s_mean',`s_sd'," ///
        "`corr',`status',`notes'" _n
end

/* ============================================================
   Test 01: Default options
   ============================================================ */
display as text _n "=== Test 01: Default options ==="
import delimited "`outdir'/test01_data.csv", clear
grf_instrumental_forest y w z x1 x2 x3 x4 x5, ///
    gen(stata_pred) ntrees(500) seed(42) replace
correlate r_pred stata_pred
local corr01 = r(rho)
local stat01 = cond(`corr01' > 0.90, "PASS", "FAIL")
write_result `fh' "01" "Default_options" `corr01' "`stat01'" ""
display "Correlation: `corr01' -- `stat01'"

/* ============================================================
   Test 02: nostabilizesplits
   ============================================================ */
display as text _n "=== Test 02: nostabilizesplits ==="
import delimited "`outdir'/test02_data.csv", clear
grf_instrumental_forest y w z x1 x2 x3 x4 x5, ///
    gen(stata_pred) ntrees(500) seed(42) nostabilizesplits replace
correlate r_pred stata_pred
local corr02 = r(rho)
local stat02 = cond(`corr02' > 0.90, "PASS", "FAIL")
write_result `fh' "02" "nostabilizesplits" `corr02' "`stat02'" ""
display "Correlation: `corr02' -- `stat02'"

/* ============================================================
   Test 03: reducedformweight(0.5)
   ============================================================ */
display as text _n "=== Test 03: reducedformweight(0.5) ==="
import delimited "`outdir'/test03_data.csv", clear
grf_instrumental_forest y w z x1 x2 x3 x4 x5, ///
    gen(stata_pred) ntrees(500) seed(42) reducedformweight(0.5) replace
correlate r_pred stata_pred
local corr03 = r(rho)
local stat03 = cond(`corr03' > 0.90, "PASS", "FAIL")
write_result `fh' "03" "reducedformweight_0.5" `corr03' "`stat03'" ""
display "Correlation: `corr03' -- `stat03'"

/* ============================================================
   Test 04: reducedformweight(1.0) -- pure reduced form
   ============================================================ */
display as text _n "=== Test 04: reducedformweight(1.0) ==="
import delimited "`outdir'/test04_data.csv", clear
grf_instrumental_forest y w z x1 x2 x3 x4 x5, ///
    gen(stata_pred) ntrees(500) seed(42) reducedformweight(1.0) replace
correlate r_pred stata_pred
local corr04 = r(rho)
local stat04 = cond(`corr04' > 0.90, "PASS", "FAIL")
write_result `fh' "04" "reducedformweight_1.0" `corr04' "`stat04'" ""
display "Correlation: `corr04' -- `stat04'"

/* ============================================================
   Test 05: User-supplied Y.hat (all three nuisance supplied)
   ============================================================ */
display as text _n "=== Test 05: User-supplied Y.hat (all three) ==="
import delimited "`outdir'/test05_data.csv", clear
grf_instrumental_forest y w z x1 x2 x3 x4 x5, ///
    gen(stata_pred) ntrees(500) seed(42) ///
    yhatinput(yhat_input) whatinput(what_input) zhatinput(zhat_input) replace
correlate r_pred stata_pred
local corr05 = r(rho)
local stat05 = cond(`corr05' > 0.90, "PASS", "FAIL")
write_result `fh' "05" "User_Yhat_supplied" `corr05' "`stat05'" "all_three_nuisance"
display "Correlation: `corr05' -- `stat05'"

/* ============================================================
   Test 06: User-supplied W.hat focus (all three, different seeds)
   ============================================================ */
display as text _n "=== Test 06: User-supplied W.hat (all three) ==="
import delimited "`outdir'/test06_data.csv", clear
grf_instrumental_forest y w z x1 x2 x3 x4 x5, ///
    gen(stata_pred) ntrees(500) seed(42) ///
    yhatinput(yhat_input) whatinput(what_input) zhatinput(zhat_input) replace
correlate r_pred stata_pred
local corr06 = r(rho)
local stat06 = cond(`corr06' > 0.90, "PASS", "FAIL")
write_result `fh' "06" "User_What_supplied" `corr06' "`stat06'" "all_three_nuisance"
display "Correlation: `corr06' -- `stat06'"

/* ============================================================
   Test 07: User-supplied Z.hat focus (all three, different seeds)
   ============================================================ */
display as text _n "=== Test 07: User-supplied Z.hat (all three) ==="
import delimited "`outdir'/test07_data.csv", clear
grf_instrumental_forest y w z x1 x2 x3 x4 x5, ///
    gen(stata_pred) ntrees(500) seed(42) ///
    yhatinput(yhat_input) whatinput(what_input) zhatinput(zhat_input) replace
correlate r_pred stata_pred
local corr07 = r(rho)
local stat07 = cond(`corr07' > 0.90, "PASS", "FAIL")
write_result `fh' "07" "User_Zhat_supplied" `corr07' "`stat07'" "all_three_nuisance"
display "Correlation: `corr07' -- `stat07'"

/* ============================================================
   Test 08: All three nuisance supplied (lm-based)
   ============================================================ */
display as text _n "=== Test 08: All three nuisance (lm-based) ==="
import delimited "`outdir'/test08_data.csv", clear
grf_instrumental_forest y w z x1 x2 x3 x4 x5, ///
    gen(stata_pred) ntrees(500) seed(42) ///
    yhatinput(yhat_input) whatinput(what_input) zhatinput(zhat_input) replace
correlate r_pred stata_pred
local corr08 = r(rho)
local stat08 = cond(`corr08' > 0.90, "PASS", "FAIL")
write_result `fh' "08" "All_nuisance_lm" `corr08' "`stat08'" "lm_nuisance"
display "Correlation: `corr08' -- `stat08'"

/* ============================================================
   Test 09: estimate.variance
   ============================================================ */
display as text _n "=== Test 09: estimate.variance ==="
import delimited "`outdir'/test09_data.csv", clear
grf_instrumental_forest y w z x1 x2 x3 x4 x5, ///
    gen(stata_pred) ntrees(500) seed(42) estimatevariance vargenerate(stata_var) replace
correlate r_pred stata_pred
local corr09 = r(rho)
correlate r_var stata_var
local corr09v = r(rho)
local stat09 = cond(`corr09' > 0.90 & `corr09v' > 0.85, "PASS", "FAIL")
write_result `fh' "09" "estimate_variance" `corr09' "`stat09'" "var_corr=`corr09v'"
display "Pred correlation: `corr09' | Var correlation: `corr09v' -- `stat09'"

/* ============================================================
   Test 10: cluster()
   ============================================================ */
display as text _n "=== Test 10: cluster() ==="
import delimited "`outdir'/test10_data.csv", clear
grf_instrumental_forest y w z x1 x2 x3 x4 x5, ///
    gen(stata_pred) ntrees(500) seed(42) cluster(cluster_id) replace
correlate r_pred stata_pred
local corr10 = r(rho)
local stat10 = cond(`corr10' > 0.90, "PASS", "FAIL")
write_result `fh' "10" "cluster" `corr10' "`stat10'" "50_clusters"
display "Correlation: `corr10' -- `stat10'"

/* ============================================================
   Test 11: weights()
   ============================================================ */
display as text _n "=== Test 11: weights() ==="
import delimited "`outdir'/test11_data.csv", clear
grf_instrumental_forest y w z x1 x2 x3 x4 x5, ///
    gen(stata_pred) ntrees(500) seed(42) weights(obs_wt) replace
correlate r_pred stata_pred
local corr11 = r(rho)
local stat11 = cond(`corr11' > 0.90, "PASS", "FAIL")
write_result `fh' "11" "observation_weights" `corr11' "`stat11'" "uniform_0.5_2.0"
display "Correlation: `corr11' -- `stat11'"

/* ============================================================
   Test 12: nohonesty
   ============================================================ */
display as text _n "=== Test 12: nohonesty ==="
import delimited "`outdir'/test12_data.csv", clear
grf_instrumental_forest y w z x1 x2 x3 x4 x5, ///
    gen(stata_pred) ntrees(500) seed(42) nohonesty replace
correlate r_pred stata_pred
local corr12 = r(rho)
local stat12 = cond(`corr12' > 0.90, "PASS", "FAIL")
write_result `fh' "12" "nohonesty" `corr12' "`stat12'" ""
display "Correlation: `corr12' -- `stat12'"

/* ============================================================
   Test 13: mtry=2 + minnodesize=20
   ============================================================ */
display as text _n "=== Test 13: mtry=2 + minnodesize=20 ==="
import delimited "`outdir'/test13_data.csv", clear
grf_instrumental_forest y w z x1 x2 x3 x4 x5, ///
    gen(stata_pred) ntrees(500) seed(42) mtry(2) minnodesize(20) replace
correlate r_pred stata_pred
local corr13 = r(rho)
local stat13 = cond(`corr13' > 0.90, "PASS", "FAIL")
write_result `fh' "13" "mtry2_minnodesize20" `corr13' "`stat13'" ""
display "Correlation: `corr13' -- `stat13'"

/* ============================================================
   Test 14: Strong instrument
   ============================================================ */
display as text _n "=== Test 14: Strong instrument ==="
import delimited "`outdir'/test14_data.csv", clear
grf_instrumental_forest y w z x1 x2 x3 x4 x5, ///
    gen(stata_pred) ntrees(500) seed(42) replace
correlate r_pred stata_pred
local corr14 = r(rho)
local stat14 = cond(`corr14' > 0.90, "PASS", "FAIL")
write_result `fh' "14" "strong_instrument" `corr14' "`stat14'" "compliance_0.95"
display "Correlation: `corr14' -- `stat14'"

/* ============================================================
   Test 15: Weak instrument
   ============================================================ */
display as text _n "=== Test 15: Weak instrument ==="
import delimited "`outdir'/test15_data.csv", clear
grf_instrumental_forest y w z x1 x2 x3 x4 x5, ///
    gen(stata_pred) ntrees(500) seed(42) replace
correlate r_pred stata_pred
local corr15 = r(rho)
/* Weak instrument: threshold 0.50 */
local stat15 = cond(`corr15' > 0.50, "PASS", "FAIL")
write_result `fh' "15" "weak_instrument" `corr15' "`stat15'" "compliance_0.05_threshold_0.50"
display "Correlation: `corr15' -- `stat15' (threshold=0.50)"

/* ============================================================
   Test 16: No heterogeneity (constant tau)
   ============================================================ */
display as text _n "=== Test 16: No heterogeneity (constant tau=1.5) ==="
import delimited "`outdir'/test16_data.csv", clear
grf_instrumental_forest y w z x1 x2 x3 x4 x5, ///
    gen(stata_pred) ntrees(500) seed(42) replace
correlate r_pred stata_pred
local corr16 = r(rho)
local stat16 = cond(`corr16' > 0.90, "PASS", "FAIL")
write_result `fh' "16" "no_heterogeneity" `corr16' "`stat16'" "constant_tau_1.5"
display "Correlation: `corr16' -- `stat16'"

/* ============================================================
   Test 17: nuisancetrees=100
   ============================================================ */
display as text _n "=== Test 17: nuisancetrees(100) ==="
import delimited "`outdir'/test17_data.csv", clear
grf_instrumental_forest y w z x1 x2 x3 x4 x5, ///
    gen(stata_pred) ntrees(500) seed(42) ///
    yhatinput(yhat_input) whatinput(what_input) zhatinput(zhat_input) replace
correlate r_pred stata_pred
local corr17 = r(rho)
local stat17 = cond(`corr17' > 0.90, "PASS", "FAIL")
write_result `fh' "17" "nuisancetrees_100" `corr17' "`stat17'" "100tree_nuisance"
display "Correlation: `corr17' -- `stat17'"

/* ============================================================
   Test 18: yhatgen + whatgen + zhatgen (save nuisance estimates)
   ============================================================ */
display as text _n "=== Test 18: Save nuisance estimates (yhatgen/whatgen/zhatgen) ==="
import delimited "`outdir'/test18_data.csv", clear
grf_instrumental_forest y w z x1 x2 x3 x4 x5, ///
    gen(stata_pred) ntrees(500) seed(42) ///
    yhatgenerate(s_yhat) whatgenerate(s_what) zhatgenerate(s_zhat) replace

/* Check nuisance vars exist and are non-missing */
local yhat_ok = 0
local what_ok = 0
local zhat_ok = 0
capture confirm variable s_yhat
if !_rc {
    quietly count if !missing(s_yhat)
    local yhat_ok = (r(N) == 500)
}
capture confirm variable s_what
if !_rc {
    quietly count if !missing(s_what)
    local what_ok = (r(N) == 500)
}
capture confirm variable s_zhat
if !_rc {
    quietly count if !missing(s_zhat)
    local zhat_ok = (r(N) == 500)
}

correlate r_pred stata_pred
local corr18 = r(rho)
correlate r_yhat s_yhat
local corr18y = r(rho)
correlate r_what s_what
local corr18w = r(rho)
correlate r_zhat s_zhat
local corr18z = r(rho)

local nuis_all_ok = (`yhat_ok' & `what_ok' & `zhat_ok')
local stat18 = cond(`corr18' > 0.90 & `nuis_all_ok', "PASS", "FAIL")
local notes18 = "yhat_corr=`corr18y'_what_corr=`corr18w'_zhat_corr=`corr18z'"
write_result `fh' "18" "save_nuisance_estimates" `corr18' "`stat18'" "`notes18'"
display "Pred corr: `corr18' | Yhat corr: `corr18y' | What corr: `corr18w' | Zhat corr: `corr18z' -- `stat18'"

/* ============================================================
   Close results file
   ============================================================ */
file close `fh'
display as text _n "=== All Stata tests complete ==="
display as text "Results written to: `outdir'/stata_results.csv"
