* Fix Test 1 - save the prediction for Test 1 (it ran successfully, just export issue)
adopath ++ /tmp/grf_stata
set more off

local outdir "/tmp/grf_stata/tests/fidelity_reports/08_causal_survival"
local horizon = 0.532911631814285

import delimited "`outdir'/test_data.csv", clear varnames(1) case(lower)

* Test 1 re-run
capture drop tau _grf_cs_what
grf_causal_survival_forest time status w x1 x2 x3 x4 x5, ///
    gen(tau) ntrees(500) seed(42) horizon(`horizon') target(1)

* Save regardless of _rc using capture
capture export delimited tau using "`outdir'/s_pred_01.csv", replace
if _rc != 0 {
    * try differently
    outsheet tau using "`outdir'/s_pred_01.csv", comma replace
}

di "Done. tau has " _N " obs"
sum tau
