* debug_survival.do -- Diagnose survival fidelity mismatch
clear all
set more off

* Load same data R used
import delimited "tests/ref/survival_input.csv", clear
describe
summarize time status

* Run survival forest matching R parameters
grf_survival_forest time status x1 x2 x3 x4 x5, gen(stata_surv) ntrees(2000) ///
    seed(42) noutput(10) replace

* Show Stata predictions (first 5 obs)
list stata_surv_s1 stata_surv_s2 stata_surv_s3 in 1/5

* Load R predictions
preserve
import delimited "tests/ref/survival_output.csv", clear
list t1 t2 t3 in 1/5
* Summary of R columns
summarize t1 t2 t3
restore

* Merge and compare
preserve
import delimited "tests/ref/survival_output.csv", clear
gen n = _n
tempfile rref
save `rref'
restore

gen n = _n
merge 1:1 n using `rref', nogenerate

* Direct comparison
display "Stata s1 vs R t1:"
summarize stata_surv_s1 t1
corr stata_surv_s1 t1
list stata_surv_s1 t1 stata_surv_s2 t2 in 1/10

* Check ranges
display "Range of Stata s1:"
summarize stata_surv_s1
display "Range of R t1:"
summarize t1
