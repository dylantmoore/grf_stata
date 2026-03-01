* Test 04: nclasses specified (explicit nclasses(2)) vs auto-detect
adopath ++ /tmp/grf_stata
import delimited "/tmp/grf_stata/tests/fidelity_reports/06_probability/test04_data.csv", clear

* Run with explicit nclasses(2)
grf_probability_forest y x1 x2 x3 x4 x5, gen(spec_pred) ntrees(500) seed(42) nclasses(2)
rename spec_pred_c0 stata_spec_c0
rename spec_pred_c1 stata_spec_c1

* Run with auto-detect (default nclasses=0 => auto)
grf_probability_forest y x1 x2 x3 x4 x5, gen(auto_pred) ntrees(500) seed(42)
rename auto_pred_c0 stata_auto_c0
rename auto_pred_c1 stata_auto_c1

keep stata_spec_c0 stata_spec_c1 stata_auto_c0 stata_auto_c1 r_pred_c0 r_pred_c1
export delimited "/tmp/grf_stata/tests/fidelity_reports/06_probability/test04_stata.csv", replace
