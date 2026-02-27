* test_generate_data.do -- Tests for grf_generate_causal_data and
*                         grf_generate_causal_survival_data
* Tests all DGPs, seed reproducibility, and error handling

clear all
set more off

local errors = 0

* ============================================================
* grf_generate_causal_data: test all 11 DGPs
* ============================================================

local causal_dgps "simple aw1 aw2 aw3 ai1 ai2 kunzel nw1 nw2 nw3 nw4"
foreach dgp of local causal_dgps {

    capture noisily {
        clear
        grf_generate_causal_data, n(200) p(5) dgp(`dgp') seed(42)

        * Save r() results immediately before any other r-class command
        local saved_N = r(N)
        local saved_p = r(p)
        local saved_dgp "`r(dgp)'"

        * Verify r() results
        assert `saved_N' == 200
        assert `saved_p' == 5
        assert "`saved_dgp'" == "`dgp'"

        * Verify variables exist
        confirm numeric variable X1
        confirm numeric variable X2
        confirm numeric variable X3
        confirm numeric variable X4
        confirm numeric variable X5
        confirm numeric variable W
        confirm numeric variable Y
        confirm numeric variable tau_true

        * W is binary
        quietly count if !inlist(W, 0, 1)
        assert r(N) == 0

        * Y is non-missing
        quietly count if missing(Y)
        assert r(N) == 0

        * tau_true is non-missing
        quietly count if missing(tau_true)
        assert r(N) == 0
    }
    if _rc {
        display as error "FAIL: grf_generate_causal_data dgp(`dgp')"
        local errors = `errors' + 1
    }
    else {
        display as result "PASS: grf_generate_causal_data dgp(`dgp')"
    }
}

* ============================================================
* Causal data: seed reproducibility
* ============================================================

* ---- Test 12: seed reproducibility (causal) ----
capture noisily {
    clear
    grf_generate_causal_data, n(200) p(5) dgp(aw1) seed(99)
    tempfile run1
    quietly save `run1'

    clear
    grf_generate_causal_data, n(200) p(5) dgp(aw1) seed(99)
    * Compare all variables to first run
    quietly cf _all using `run1'
}
if _rc {
    display as error "FAIL: causal data seed reproducibility"
    local errors = `errors' + 1
}
else {
    display as result "PASS: causal data seed reproducibility"
}

* ============================================================
* Causal data: error cases
* ============================================================

* ---- Test 13: minimum p enforcement -- dgp(nw1) with p(2) should error ----
capture noisily {
    clear
    capture grf_generate_causal_data, n(200) p(2) dgp(nw1) seed(42)
    assert _rc != 0
}
if _rc {
    display as error "FAIL: causal data minimum p enforcement (nw1 with p=2)"
    local errors = `errors' + 1
}
else {
    display as result "PASS: causal data minimum p enforcement (nw1 with p=2)"
}

* ---- Test 14: minimum p enforcement -- dgp(simple) with p(2) should error ----
capture noisily {
    clear
    capture grf_generate_causal_data, n(200) p(2) dgp(simple) seed(42)
    assert _rc != 0
}
if _rc {
    display as error "FAIL: causal data minimum p enforcement (simple with p=2)"
    local errors = `errors' + 1
}
else {
    display as result "PASS: causal data minimum p enforcement (simple with p=2)"
}

* ---- Test 15: invalid DGP name should error ----
capture noisily {
    clear
    capture grf_generate_causal_data, n(200) p(5) dgp(bogus) seed(42)
    assert _rc != 0
}
if _rc {
    display as error "FAIL: causal data invalid DGP name"
    local errors = `errors' + 1
}
else {
    display as result "PASS: causal data invalid DGP name"
}

* ============================================================
* grf_generate_causal_survival_data: test all 6 DGPs
* ============================================================

local surv_dgps "simple1 type1 type2 type3 type4 type5"
foreach dgp of local surv_dgps {

    capture noisily {
        clear
        grf_generate_causal_survival_data, n(200) p(5) dgp(`dgp') seed(42)

        * Save r() results immediately before count overwrites them
        local saved_n_events = r(n_events)
        local saved_n_censored = r(n_censored)

        * Verify variables exist
        confirm numeric variable X1
        confirm numeric variable X2
        confirm numeric variable X3
        confirm numeric variable X4
        confirm numeric variable X5
        confirm numeric variable W
        confirm numeric variable T
        confirm numeric variable D
        confirm numeric variable tau_true

        * T > 0 for all observations
        quietly count if T <= 0
        assert r(N) == 0

        * D is binary (0/1)
        quietly count if !inlist(D, 0, 1)
        assert r(N) == 0

        * tau_true is non-missing
        quietly count if missing(tau_true)
        assert r(N) == 0

        * There should be events
        assert `saved_n_events' > 0
    }
    if _rc {
        display as error "FAIL: grf_generate_causal_survival_data dgp(`dgp')"
        local errors = `errors' + 1
    }
    else {
        display as result "PASS: grf_generate_causal_survival_data dgp(`dgp')"
    }
}

* ============================================================
* Survival data: seed reproducibility
* ============================================================

* ---- Test 22: seed reproducibility (survival) ----
capture noisily {
    clear
    grf_generate_causal_survival_data, n(200) p(5) dgp(type1) seed(99)
    tempfile srun1
    quietly save `srun1'

    clear
    grf_generate_causal_survival_data, n(200) p(5) dgp(type1) seed(99)
    quietly cf _all using `srun1'
}
if _rc {
    display as error "FAIL: survival data seed reproducibility"
    local errors = `errors' + 1
}
else {
    display as result "PASS: survival data seed reproducibility"
}

* ============================================================
* Summary
* ============================================================
display ""
if `errors' > 0 {
    display as error "FAILED: `errors' errors in data generation tests"
    exit 1
}
else {
    display as result "ALL DATA GENERATION TESTS PASSED"
}
