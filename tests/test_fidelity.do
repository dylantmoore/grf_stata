* test_fidelity.do -- R-vs-Stata numerical fidelity tests
* Compares Stata grf results to R grf reference data in tests/ref/
*
* Each comparison gracefully skips if reference files are missing,
* allowing the test suite to pass even if R has not been run yet.
*
* Run: do tests/test_fidelity.do

clear all
set more off

local errors = 0

* ============================================================
* 1. Regression OOB predictions (correlation > 0.95)
* ============================================================
capture confirm file "ref/regression_input.csv"
if _rc {
    display as text "SKIP: regression predictions (ref file missing)"
}
else {
    capture noisily {
        import delimited using "ref/regression_input.csv", clear
        grf_regression_forest y x1-x5, gen(pred_stata) ntrees(2000) seed(42)

        * Load R predictions
        preserve
        import delimited using "ref/regression_output.csv", clear
        rename prediction pred_r
        tempfile r_preds
        save `r_preds'
        restore

        merge 1:1 _n using `r_preds', nogenerate
        correlate pred_stata pred_r
        assert r(rho) > 0.95
    }
    if _rc {
        display as error "FAIL: regression predictions (correlation with R)"
        local errors = `errors' + 1
    }
    else {
        display as result "PASS: regression predictions (correlation with R)"
    }
}

* ============================================================
* 2. Causal ATE (all) -- z-test < 3
* ============================================================
capture confirm file "ref/causal_input.csv"
if _rc {
    display as text "SKIP: causal ATE all (ref file missing)"
}
else {
    capture noisily {
        import delimited using "ref/causal_input.csv", clear
        grf_causal_forest y w x1-x5, gen(cate_fid) ntrees(2000) seed(42)
        grf_ate
        local stata_ate = r(ate)
        local stata_se  = r(se)

        * Load R ATE
        preserve
        import delimited using "ref/causal_ate.csv", clear
        local r_ate = estimate[1]
        local r_se  = std_err[1]
        restore

        * z-test: |stata - R| / sqrt(se_stata^2 + se_R^2)
        local z = abs(`stata_ate' - `r_ate') / sqrt(`stata_se'^2 + `r_se'^2)
        assert `z' < 3
    }
    if _rc {
        display as error "FAIL: causal ATE (all)"
        local errors = `errors' + 1
    }
    else {
        display as result "PASS: causal ATE (all)"
    }
}

* ============================================================
* 3. Causal ATE (treated) -- z-test < 3
* ============================================================
capture confirm file "ref/causal_ate_treated.csv"
if _rc {
    display as text "SKIP: causal ATE treated (ref file missing)"
}
else {
    capture noisily {
        * Reuse data from previous test; re-fit if needed
        capture confirm variable cate_fid
        if _rc {
            import delimited using "ref/causal_input.csv", clear
            grf_causal_forest y w x1-x5, gen(cate_fid) ntrees(2000) seed(42)
        }
        grf_ate, targetsample(treated)
        local stata_ate = r(ate)
        local stata_se  = r(se)

        preserve
        import delimited using "ref/causal_ate_treated.csv", clear
        local r_ate = estimate[1]
        local r_se  = std_err[1]
        restore

        local z = abs(`stata_ate' - `r_ate') / sqrt(`stata_se'^2 + `r_se'^2)
        assert `z' < 3
    }
    if _rc {
        display as error "FAIL: causal ATE (treated)"
        local errors = `errors' + 1
    }
    else {
        display as result "PASS: causal ATE (treated)"
    }
}

* ============================================================
* 4. Causal ATE (control) -- z-test < 3
* ============================================================
capture confirm file "ref/causal_ate_control.csv"
if _rc {
    display as text "SKIP: causal ATE control (ref file missing)"
}
else {
    capture noisily {
        capture confirm variable cate_fid
        if _rc {
            import delimited using "ref/causal_input.csv", clear
            grf_causal_forest y w x1-x5, gen(cate_fid) ntrees(2000) seed(42)
        }
        grf_ate, targetsample(control)
        local stata_ate = r(ate)
        local stata_se  = r(se)

        preserve
        import delimited using "ref/causal_ate_control.csv", clear
        local r_ate = estimate[1]
        local r_se  = std_err[1]
        restore

        local z = abs(`stata_ate' - `r_ate') / sqrt(`stata_se'^2 + `r_se'^2)
        assert `z' < 3
    }
    if _rc {
        display as error "FAIL: causal ATE (control)"
        local errors = `errors' + 1
    }
    else {
        display as result "PASS: causal ATE (control)"
    }
}

* ============================================================
* 5. Causal ATE (overlap) -- z-test < 3
* ============================================================
capture confirm file "ref/causal_ate_overlap.csv"
if _rc {
    display as text "SKIP: causal ATE overlap (ref file missing)"
}
else {
    capture noisily {
        capture confirm variable cate_fid
        if _rc {
            import delimited using "ref/causal_input.csv", clear
            grf_causal_forest y w x1-x5, gen(cate_fid) ntrees(2000) seed(42)
        }
        grf_ate, targetsample(overlap)
        local stata_ate = r(ate)
        local stata_se  = r(se)

        preserve
        import delimited using "ref/causal_ate_overlap.csv", clear
        local r_ate = estimate[1]
        local r_se  = std_err[1]
        restore

        local z = abs(`stata_ate' - `r_ate') / sqrt(`stata_se'^2 + `r_se'^2)
        assert `z' < 3
    }
    if _rc {
        display as error "FAIL: causal ATE (overlap)"
        local errors = `errors' + 1
    }
    else {
        display as result "PASS: causal ATE (overlap)"
    }
}

* ============================================================
* 6. Causal BLP coefficients -- per-coefficient z-test < 3
* ============================================================
capture confirm file "ref/causal_blp.csv"
if _rc {
    display as text "SKIP: causal BLP coefficients (ref file missing)"
}
else {
    capture noisily {
        * Re-fit causal forest (BLP replaces e())
        import delimited using "ref/causal_input.csv", clear
        grf_causal_forest y w x1-x5, gen(cate_blp) ntrees(2000) seed(42)
        grf_best_linear_projection x1-x5

        * Save Stata coefficients (intercept + x1..x5)
        matrix b_stata = e(b)
        matrix V_stata = e(V)

        * Load R BLP results
        preserve
        import delimited using "ref/causal_blp.csv", clear
        * R file has: variable, estimate, std.error, t.value, p.value
        * Row 1 = intercept, rows 2-6 = x1..x5
        * Stata column order: x1 x2 x3 x4 x5 _cons
        local blp_ok 1
        forvalues j = 1/5 {
            local r_est = estimate[`j' + 1]
            local r_se  = std_error[`j' + 1]
            local s_est = b_stata[1, `j']
            local s_se  = sqrt(V_stata[`j', `j'])
            local z = abs(`s_est' - `r_est') / sqrt(`s_se'^2 + `r_se'^2)
            if `z' >= 3 {
                local blp_ok 0
            }
        }
        * Also check intercept (last column in Stata = _cons, first row in R)
        local ncols = colsof(b_stata)
        local r_est = estimate[1]
        local r_se  = std_error[1]
        local s_est = b_stata[1, `ncols']
        local s_se  = sqrt(V_stata[`ncols', `ncols'])
        local z = abs(`s_est' - `r_est') / sqrt(`s_se'^2 + `r_se'^2)
        if `z' >= 3 {
            local blp_ok 0
        }
        restore

        assert `blp_ok' == 1
    }
    if _rc {
        display as error "FAIL: causal BLP coefficients"
        local errors = `errors' + 1
    }
    else {
        display as result "PASS: causal BLP coefficients"
    }
}

* ============================================================
* 7. Variable importance rank correlation (Spearman > 0.70)
* ============================================================
capture confirm file "ref/causal_variable_importance.csv"
if _rc {
    display as text "SKIP: variable importance rank correlation (ref file missing)"
}
else {
    capture noisily {
        import delimited using "ref/causal_input.csv", clear
        grf_variable_importance y x1-x5, ntrees(2000) seed(42)
        matrix vi_stata = r(importance)

        * Load R variable importance
        preserve
        import delimited using "ref/causal_variable_importance.csv", clear

        * Create a dataset with both R and Stata importance for Spearman rank
        gen vi_stata = .
        forvalues j = 1/5 {
            replace vi_stata = vi_stata[1, `j'] if _n == `j'
        }

        * Compute Spearman rank correlation manually
        * Rank R importance
        egen rank_r = rank(importance)
        * Rank Stata importance
        egen rank_s = rank(vi_stata)
        correlate rank_r rank_s
        local spearman = r(rho)
        restore

        assert `spearman' > 0.70
    }
    if _rc {
        display as error "FAIL: variable importance rank correlation"
        local errors = `errors' + 1
    }
    else {
        display as result "PASS: variable importance rank correlation"
    }
}

* ============================================================
* 8. DR scores correlation (> 0.90)
* ============================================================
capture confirm file "ref/causal_scores.csv"
if _rc {
    display as text "SKIP: DR scores correlation (ref file missing)"
}
else {
    capture noisily {
        import delimited using "ref/causal_input.csv", clear
        grf_causal_forest y w x1-x5, gen(cate_sc) ntrees(2000) seed(42)
        grf_get_scores, gen(scores_stata)

        * Load R DR scores
        preserve
        import delimited using "ref/causal_scores.csv", clear
        rename score score_r
        tempfile r_scores
        save `r_scores'
        restore

        merge 1:1 _n using `r_scores', nogenerate
        correlate scores_stata score_r
        assert r(rho) > 0.90
    }
    if _rc {
        display as error "FAIL: DR scores correlation"
        local errors = `errors' + 1
    }
    else {
        display as result "PASS: DR scores correlation"
    }
}

* ============================================================
* Summary
* ============================================================
display ""
if `errors' > 0 {
    display as error "FAILED: `errors' fidelity test(s) did not pass"
    exit 1
}
else {
    display as result "ALL FIDELITY TESTS PASSED"
}
