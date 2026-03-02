* run_all.do -- Standard test runner for grf_stata
* Run from project root: do tests/run_all.do

clear all
set more off

local errors = 0

* Resolve test path prefix for both invocation modes:
*   1) repo root:   do tests/run_all.do
*   2) tests dir:   do run_all.do
local test_prefix ""
capture confirm file "tests/test_generate_data.do"
if !_rc {
    local test_prefix "tests/"
}
else {
    capture confirm file "test_generate_data.do"
    if _rc {
        display as error "Could not locate test files from current working directory."
        exit 601
    }
    adopath ++ ".."
}

local test_files ///
    "`test_prefix'test_parity_scope.do" ///
    "`test_prefix'test_generate_data.do" ///
    "`test_prefix'test_seed_reproducibility.do" ///
    "`test_prefix'test_model_id.do" ///
    "`test_prefix'test_options_survival.do" ///
    "`test_prefix'test_options_causal_survival.do" ///
    "`test_prefix'test_get_scores.do" ///
    "`test_prefix'test_forest_utilities.do" ///
    "`test_prefix'test_rate.do" ///
    "`test_prefix'test_fidelity.do" ///
    "`test_prefix'test_fidelity_gaps.do" ///
    "`test_prefix'test_full_pipeline.do"

foreach tf of local test_files {
    display as text ""
    display as text "============================================================"
    display as text "Running: `tf'"
    display as text "============================================================"

    capture noisily do `tf'
    if _rc {
        display as error "FAIL: `tf'"
        local errors = `errors' + 1
    }
    else {
        display as result "PASS: `tf'"
    }
}

display as text ""
if `errors' > 0 {
    display as error "FAILED: `errors' test file(s) failed"
    exit 1
}
else {
    display as result "ALL STANDARD TESTS PASSED"
}
