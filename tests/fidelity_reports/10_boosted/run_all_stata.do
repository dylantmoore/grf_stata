* run_all_stata.do -- Run all 16 boosted regression forest fidelity tests
set more off

local testdir "/tmp/grf_stata/tests/fidelity_reports/10_boosted"

local tests "test01_default test02_booststeps1 test03_booststeps3 test04_booststeps5 test05_boostmaxsteps10 test06_errreduction090 test07_errreduction099 test08_boosttreestune50 test09_nostabilize test10_cluster test11_weights test12_nohonesty test13_mtry2 test14_vs_regression test15_linear_data test16_highly_nonlinear"

local errors = 0
local total  = 0

foreach t of local tests {
    local total = `total' + 1
    display as text _newline "======================================"
    display as text "Running: `t'"
    display as text "======================================"
    capture noisily do "`testdir'/`t'.do"
    if _rc {
        display as error "FAILED: `t' (rc=`_rc')"
        local errors = `errors' + 1
    }
    else {
        display as result "PASSED: `t'"
    }
}

display _newline
display "======================================"
display "SUMMARY: `errors' errors out of `total' tests"
display "======================================"
if `errors' > 0 {
    exit 1
}
