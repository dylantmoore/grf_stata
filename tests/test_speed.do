* test_speed.do -- Speed comparison tests for GRF Stata plugin
* Run from the test-c-plugin-skill-grf/ directory
*
* Benchmarks all 12 forest types at different dataset sizes.
* Reports wall-clock time for each.

clear all
set more off

display as text ""
display as text "============================================================"
display as text " GRF Speed Benchmarks"
display as text "============================================================"
display as text " Date: `c(current_date)' `c(current_time)'"
display as text " Stata: `c(stata_version)' (`c(processors)' processors)"
display as text "============================================================"

* ---- Helper program to time a command ----
capture program drop benchmark
program define benchmark
    syntax, CMD(string) NAME(string) [N(integer 0)]

    timer clear 1
    timer on 1
    `cmd'
    timer off 1
    quietly timer list 1
    local elapsed = r(t1)
    if `n' > 0 {
        display as text "  " _col(5) as result "`name'" _col(45) as text %8.2f `elapsed' " sec" _col(60) "  (n=" as result `n' as text ")"
    }
    else {
        display as text "  " _col(5) as result "`name'" _col(45) as text %8.2f `elapsed' " sec"
    }
end

* ================================================================
* PART 1: Forest type comparison (n=2000, p=10, trees=500)
* ================================================================
display as text ""
display as text "--- Part 1: Forest Type Comparison (n=2000, p=10, trees=500) ---"
display as text ""

clear
set obs 2000
set seed 42

* Generate covariates
forvalues j = 1/10 {
    gen x`j' = rnormal()
}

* Outcome with nonlinear DGP
gen y = 2*x1 + x2^2 + 0.5*x3*x4 + rnormal()

* Treatment
gen w = (x1 + rnormal() > 0)

* Instrument
gen z = (x2 + rnormal() > 0)

* Treatment effect
replace y = y + w*(1 + x5) + z*0.3

* Second outcome for multi-regression
gen y2 = x1 + x3 + rnormal()

* Survival data
gen double time_surv = exp(0.5*x1 + rnormal())
gen status_surv = (runiform() > 0.3)

* Multi-arm treatment
gen w_multi = floor(3*runiform())

local xvars "x1 x2 x3 x4 x5 x6 x7 x8 x9 x10"

display as text _col(5) "Forest Type" _col(45) "Time" _col(60) "  Obs"
display as text _col(5) "{hline 55}"

* 1. Regression forest
benchmark, cmd("grf_regression_forest y `xvars', gen(p_reg) ntrees(500) seed(42)") ///
    name("Regression Forest") n(2000)

* 2. Causal forest
benchmark, cmd("grf_causal_forest y w `xvars', gen(p_caus) ntrees(500) seed(42)") ///
    name("Causal Forest") n(2000)

* 3. Quantile forest (inline timing — commas in quantiles() break cmd() parsing)
timer clear 1
timer on 1
grf_quantile_forest y x1 x2 x3 x4 x5 x6 x7 x8 x9 x10, gen(p_q) ntrees(500) seed(42) quantiles(0.25 0.5 0.75)
timer off 1
quietly timer list 1
display as text "  " _col(5) as result "Quantile Forest (3 quantiles)" _col(45) as text %8.2f r(t1) " sec" _col(60) "  (n=" as result "2000" as text ")"

* 4. Instrumental forest
benchmark, cmd("grf_instrumental_forest y w z `xvars', gen(p_inst) ntrees(500) seed(42)") ///
    name("Instrumental Forest") n(2000)

* 5. Probability forest
gen w_binary = (w > 0)
benchmark, cmd("grf_probability_forest w_binary `xvars', gen(p_prob) ntrees(500) seed(42)") ///
    name("Probability Forest") n(2000)

* 6. Survival forest
benchmark, cmd("grf_survival_forest time_surv status_surv `xvars', gen(p_surv) ntrees(500) seed(42)") ///
    name("Survival Forest") n(2000)

* 7. Causal survival forest
benchmark, cmd("grf_causal_survival_forest time_surv status_surv w `xvars', gen(p_csurv) ntrees(500) seed(42)") ///
    name("Causal Survival Forest") n(2000)

* 8. Multi-arm causal forest
benchmark, cmd("grf_multi_arm_causal_forest y w_multi `xvars', gen(p_ma) ntrees(500) seed(42) ntreat(3)") ///
    name("Multi-arm Causal Forest") n(2000)

* 9. Multi-regression forest
benchmark, cmd("grf_multi_regression_forest y y2 `xvars', gen(p_mr) ndep(2) ntrees(500) seed(42)") ///
    name("Multi-regression Forest") n(2000)

* 10. Local linear regression forest
benchmark, cmd("grf_ll_regression_forest y `xvars', gen(p_ll) ntrees(500) seed(42)") ///
    name("LL Regression Forest") n(2000)

* 11. LL regression with LL splits
benchmark, cmd("grf_ll_regression_forest y `xvars', gen(p_lls) ntrees(500) seed(42) llsplit replace") ///
    name("LL Regression (LL splits)") n(2000)

* 12. Boosted regression forest
benchmark, cmd("grf_boosted_regression_forest y `xvars', gen(p_boost) ntrees(500) seed(42) booststeps(3)") ///
    name("Boosted Regression (3 steps)") n(2000)

* 13. LM forest
benchmark, cmd("grf_lm_forest y w, gen(p_lm) xvars(`xvars') ntrees(500) seed(42)") ///
    name("LM Forest (1 regressor)") n(2000)


* ================================================================
* PART 2: Scaling with dataset size (regression forest)
* ================================================================
display as text ""
display as text "--- Part 2: Scaling with N (Regression Forest, p=10, trees=500) ---"
display as text ""
display as text _col(5) "N" _col(25) "Time (sec)" _col(40) "Time/N (ms)"
display as text _col(5) "{hline 50}"

foreach nn in 500 1000 2000 5000 10000 {
    clear
    set obs `nn'
    set seed 42

    forvalues j = 1/10 {
        gen x`j' = rnormal()
    }
    gen y = 2*x1 + x2^2 + rnormal()

    timer clear 1
    timer on 1
    grf_regression_forest y x1-x10, gen(pred) ntrees(500) seed(42)
    timer off 1
    quietly timer list 1
    local elapsed = r(t1)
    local per_obs = `elapsed' / `nn' * 1000
    display as text _col(5) as result %6.0f `nn' _col(25) as text %8.3f `elapsed' _col(40) %8.3f `per_obs'
}


* ================================================================
* PART 3: Scaling with number of trees
* ================================================================
display as text ""
display as text "--- Part 3: Scaling with Trees (n=2000, p=10) ---"
display as text ""
display as text _col(5) "Trees" _col(25) "Time (sec)" _col(40) "Time/Tree (ms)"
display as text _col(5) "{hline 50}"

clear
set obs 2000
set seed 42
forvalues j = 1/10 {
    gen x`j' = rnormal()
}
gen y = 2*x1 + x2^2 + rnormal()

foreach nt in 100 500 1000 2000 {
    timer clear 1
    timer on 1
    grf_regression_forest y x1-x10, gen(pred) ntrees(`nt') seed(42) replace
    timer off 1
    quietly timer list 1
    local elapsed = r(t1)
    local per_tree = `elapsed' / `nt' * 1000
    display as text _col(5) as result %6.0f `nt' _col(25) as text %8.3f `elapsed' _col(40) %8.3f `per_tree'
}


* ================================================================
* PART 4: Scaling with number of covariates
* ================================================================
display as text ""
display as text "--- Part 4: Scaling with p (n=2000, trees=500) ---"
display as text ""
display as text _col(5) "p" _col(25) "Time (sec)" _col(40) "mtry"
display as text _col(5) "{hline 50}"

foreach pp in 5 10 20 50 {
    clear
    set obs 2000
    set seed 42

    forvalues j = 1/`pp' {
        gen x`j' = rnormal()
    }
    gen y = 2*x1 + rnormal()

    timer clear 1
    timer on 1
    grf_regression_forest y x1-x`pp', gen(pred) ntrees(500) seed(42)
    timer off 1
    quietly timer list 1
    local elapsed = r(t1)
    local used_mtry = e(mtry)
    display as text _col(5) as result %6.0f `pp' _col(25) as text %8.3f `elapsed' _col(40) %6.0f `used_mtry'
}


* ================================================================
* PART 5: Causal forest with/without nuisance estimation
* ================================================================
display as text ""
display as text "--- Part 5: Causal Forest Breakdown (n=2000, p=10, trees=500) ---"
display as text ""

clear
set obs 2000
set seed 42
forvalues j = 1/10 {
    gen x`j' = rnormal()
}
gen w = (x1 + rnormal() > 0)
gen y = 2*x1 + w*(1 + x5) + rnormal()

* Full causal forest (includes nuisance estimation)
timer clear 1
timer on 1
grf_causal_forest y w x1-x10, gen(tau_full) ntrees(500) seed(42)
timer off 1
quietly timer list 1
display as text "  Full causal forest (with nuisance):   " as result %8.3f r(t1) " sec"

* Causal forest with fewer nuisance trees
timer clear 1
timer on 1
grf_causal_forest y w x1-x10, gen(tau_fast) ntrees(500) seed(42) nuisancetrees(100) replace
timer off 1
quietly timer list 1
display as text "  With nuisancetrees(100):              " as result %8.3f r(t1) " sec"

* Summary
display as text ""
display as text "============================================================"
display as text " Speed benchmarks completed"
display as text "============================================================"
display as text ""
