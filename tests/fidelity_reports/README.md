# Fidelity Test Suites: R `grf` 2.5.0 vs Stata `grf_stata`

This directory contains 15 independent test suites that numerically compare every
`grf_stata` command against its R `grf` equivalent. Each suite generates shared
data in R, runs the matching Stata command, and computes Pearson correlations
between R and Stata output.

**Total: 295 tests across 15 suites. See `CONSOLIDATED_DISCREPANCIES.md` for the
full discrepancy report.**

---

## Quick Start

```bash
# Run all 15 suites end-to-end (takes ~30-60 minutes)
cd /tmp/grf_stata/tests/fidelity_reports
bash run_all_suites.sh

# Run specific suites by ID
bash run_all_suites.sh 01 07 14

# Set custom Stata path
STATA=/path/to/stata-mp bash run_all_suites.sh
```

### Prerequisites

- R 4.5+ with `grf` 2.5.0 installed (`install.packages("grf")`)
- Stata MP (StataNow 19.5+)
- `grf_stata` package at `/tmp/grf_stata` with compiled plugin

---

## Suite Index

| ID | Directory | Command(s) Tested | Tests | Pass Rate | Report |
|----|-----------|-------------------|:-----:|:---------:|--------|
| 01 | `01_regression/` | regression_forest | 21 | 100% | `01_regression_forest.md` |
| 02 | `02_causal_fit/` | causal_forest (fitting) | 19 | 95% | `02_causal_forest_fit.md` |
| 03 | `03_ate/` | average_treatment_effect (AIPW+TMLE) | 19 | 84% | `03_ate.md` |
| 04 | `04_instrumental/` | instrumental_forest | 18 | 94% | `04_instrumental_forest.md` |
| 05 | `05_quantile/` | quantile_forest | 17 | 88% | `05_quantile_forest.md` |
| 06 | `06_probability/` | probability_forest | 13 | 100% | `06_probability_forest.md` |
| 07 | `07_survival/` | survival_forest, expected_survival | 17 | 41%* | `07_survival_expected.md` |
| 08 | `08_causal_survival/` | causal_survival_forest | 16 | 19%** | `08_causal_survival_forest.md` |
| 09 | `09_multi_arm/` | multi_arm_causal_forest | 36 | 92% | `09_multi_arm_causal_forest.md` |
| 10 | `10_boosted/` | boosted_regression_forest | 17 | 100% | `10_boosted_regression_forest.md` |
| 11 | `11_ll_regression/` | ll_regression_forest | 18 | 100% | `11_ll_regression_forest.md` |
| 12 | `12_lm_forest/` | lm_forest | 17 | 100% | `12_lm_forest.md` |
| 13 | `13_multi_regression/` | multi_regression_forest | 14 | 100% | `13_multi_regression_forest.md` |
| 14 | `14_blp_calibration/` | BLP, test_calibration, get_scores | 29 | 93% | `14_blp_calibration_scores.md` |
| 15 | `15_vi_rate_tune/` | variable_importance, RATE, tune | 24 | 58% | `15_vi_rate_tune.md` |

\* Suite 07: 9 "failures" are a narrow-variance statistical artifact at early survival times.
\** Suite 08: Low pass rate reflects a known architectural gap (simplified IPCW pipeline).

---

## Pipeline Architecture

Each suite follows a 3-step pipeline:

```
Step 1: R generates data + R predictions → CSV files
        (shared DGP with set.seed(42), all test configurations)

Step 2: Stata loads CSVs, runs matching commands → exports Stata predictions
        (adopath ++ /tmp/grf_stata, loads plugin, writes *_stata.csv)

Step 3: R computes Pearson correlations between R and Stata outputs
        (reads paired CSVs, produces summary table)
```

### File naming conventions

| Pattern | Contents |
|---------|----------|
| `testNN_data.csv` or `testNN_*.csv` | Shared data with R predictions (X, Y, W, r_pred, ...) |
| `testNN_stata.csv` | Stata predictions exported for comparison |
| `*_results.csv` | Machine-readable correlation results |
| `run_all_r.R` or `run_r_tests.R` | R reference script (Step 1) |
| `run_all_stata.do` or `*_stata.do` | Stata test script (Step 2) |
| `compare*.R` or `compute_correlations.R` | Comparison script (Step 3) |

---

## Suite Details

### Suite 01: regression_forest (21 tests)

**Scripts:** Individual R/do pairs per test + `compare_all.R`

| Test | Option(s) | DGP |
|------|-----------|-----|
| 01 | Default | Y = X1 + 2*X2 + eps |
| 02 | nohonesty | same |
| 03 | mtry=2 | same |
| 04 | minnodesize=20 | same |
| 05 | samplefrac=0.3 | same |
| 06 | honestyfrac=0.7 | same |
| 07 | alpha=0.15 | same |
| 08 | imbalancepenalty=1.0 | same |
| 09 | estimatevariance + cigroupsize=2 | same |
| 10 | cluster() | same + 50 clusters |
| 11 | weights() | same + Uniform weights |
| 12 | equalizeclusterweights | same + clusters |
| 13 | nomia | same (complete data) |
| 14 | cigroupsize=2 | same |
| 15 | Combined: cluster+weights+nohonesty | same |
| 16 | Combined: mtry=2+minnodesize=20+samplefrac=0.3 | same |
| 17 | Large p (p=20, n=1000) | Y = X1+X2+eps |
| 18a | ntrees=100 | same |
| 18b | ntrees=2000 | same |
| 19 | Missing data | same + 5% MCAR |
| 20 | Seed reproducibility | same, seed=42 vs seed=42 |

---

### Suite 02: causal_forest fitting (19 tests)

**Scripts:** `run_r_tests.R`, `run_stata_tests.do`

| Test | Option(s) | Notes |
|------|-----------|-------|
| 01 | Default | tau = X1+X2 |
| 02 | nostabilizesplits | |
| 03 | nuisancetrees=100 | R: pre-computed 100-tree nuisance |
| 04 | yhatinput + whatinput | Full user-supplied nuisance |
| 05 | nohonesty | |
| 06 | mtry=2 | |
| 07 | minnodesize=20 | |
| 08 | samplefrac=0.3 | |
| 09 | cluster() | 50 clusters |
| 10 | weights() | Uniform(0.5, 1.5) |
| 11 | equalizeclusterweights | |
| 12 | Combined: nohonesty+mtry=2+minnodesize=20 | |
| 13 | ntrees=2000 | |
| 14 | estimatevariance + cigroupsize=2 | |
| 15 | Large p (p=20, n=1000) | |
| 16 | Strong heterogeneity tau=5*X1 | |
| 17 | Homogeneous tau=2 | |
| 18 | Binary+continuous covariates | |
| 19 | Null effect tau=0 | |

---

### Suite 03: average_treatment_effect (19 tests)

**Scripts:** `run_r.R`, `run_stata.do`, `compare.R`

| Test | Method | Option(s) | Notes |
|------|--------|-----------|-------|
| 01 | AIPW | ATE (all) | |
| 02 | AIPW | ATT (treated) | |
| 03 | AIPW | ATC (control) | |
| 04 | AIPW | overlap | |
| 05 | AIPW | debiasing.weights | **D-3: SE differs** |
| 06 | AIPW | clusters | |
| 07 | AIPW | equalizeclusterweights | |
| 08 | AIPW | ntrees=4000 | |
| 09 | AIPW | Strong heterogeneity | |
| 10 | AIPW | Null effect | |
| 11 | TMLE | ATE (all) | |
| 12 | TMLE | ATT (treated) | |
| 13 | TMLE | ATC (control) | |
| 14 | TMLE | overlap | |
| 15 | TMLE | Strong heterogeneity | |
| 16 | TMLE | clusters | **D-1: Stata should error** |
| 17 | TMLE | debiasing.weights | **D-18: Stata correctly errors** |
| 18 | AIPW | Heteroscedastic DGP | |
| 19 | AIPW | High-dimensional (p=20) | |

---

### Suite 04: instrumental_forest (18 tests)

**Scripts:** `run_all_r.R`, `run_all_stata.do`

| Test | Option(s) | Notes |
|------|-----------|-------|
| 01 | Default | Binary instrument |
| 02 | nohonesty | |
| 03 | mtry=2 | |
| 04 | minnodesize=20 | |
| 05 | samplefrac=0.3 | |
| 06 | cluster() | |
| 07 | weights() | |
| 08 | nostabilizesplits | |
| 09 | reducedformweight=0.5 | |
| 10 | reducedformweight=1.0 | |
| 11 | estimatevariance + cigroupsize=2 | |
| 12 | Combined options | |
| 13 | All nuisance inputs (yhat+what+zhat) | |
| 14 | ntrees=2000 | |
| 15 | Large p (p=20) | |
| 16 | Weak instrument | |
| 17 | Strong instrument | |
| 18 | Continuous instrument | |

---

### Suite 05: quantile_forest (17 tests)

**Scripts:** Individual R/do pairs + `analyze_results.R`

| Test | Option(s) | Notes |
|------|-----------|-------|
| 01 | Default (q=0.5) | |
| 02 | Single quantile (q=0.25) | |
| 03 | Many quantiles (0.1,0.25,0.5,0.75,0.9) | |
| 04 | Extreme quantiles (0.01,0.99) | **D-22: statistical limitation** |
| 05 | regression.splitting | |
| 06 | regression.splitting + multi-quantile | |
| 07 | Heteroscedastic DGP | |
| 08 | Skewed DGP (lognormal) | |
| 09 | cluster() | |
| 10 | weights() | Stata extension |
| 11 | nohonesty | |
| 12 | mtry=2 | |
| 13 | minnodesize=20 | |
| 14 | samplefrac=0.3 | |
| 15 | Combined options | |
| 16 | Large p (p=20) | |
| 17 | Crossing check (quantile ordering) | |

---

### Suite 06: probability_forest (13 tests)

**Scripts:** Individual R/do pairs + `compare_all.R`

| Test | Option(s) | Notes |
|------|-----------|-------|
| 01 | Binary (2 classes) | |
| 02 | 3 classes | |
| 03 | 5 classes | |
| 04 | nclasses explicit | |
| 05 | Unbalanced classes | |
| 06 | cluster() | |
| 07 | weights() | |
| 08 | nohonesty | |
| 09 | mtry=2 | |
| 10 | minnodesize=20 | |
| 11 | Combined options | |
| 12 | Probability sum = 1 check | |
| 13 | Large sample (n=2000) | |

---

### Suite 07: survival_forest + expected_survival (17 tests)

**Scripts:** `run_all_r.R`, `run_all_stata.do`, `analyze_results.R`

| Test | Option(s) | Notes |
|------|-----------|-------|
| 01 | Default (KM, noutput=20) | |
| 02 | noutput=50 | |
| 03 | noutput=100 | |
| 04 | Nelson-Aalen (predtype=0) | |
| 05 | KM explicit (predtype=1) | |
| 06 | nofastlogrank | Clean PASS (matches R default) |
| 07 | cluster() | |
| 08 | weights() | |
| 09 | nohonesty | |
| 10 | mtry=2 | Best results (min r=0.96) |
| 11 | minnodesize=20 | |
| 12 | High event rate (~94%) | |
| 13 | Low event rate (~34%) | |
| 14 | Expected survival (base) | r=0.9898 |
| 15 | Expected survival consistency | Directional check |
| 16 | numfailures(50) | **D-13: API design diff** |
| 17 | Combined (noutput=50+NA+cluster) | |

---

### Suite 08: causal_survival_forest (16 tests)

**Scripts:** `run_r_tests.R`, `run_stata_tests.do`, `compute_correlations.R`

| Test | Option(s) | Notes |
|------|-----------|-------|
| 01 | Default (RMST) | |
| 02 | horizon = 75th percentile | |
| 03 | horizon = 90th percentile | |
| 04 | nohonesty | |
| 05 | mtry=2 | |
| 06 | cluster() | **D-4: cluster not propagated** |
| 07 | weights() | |
| 08 | samplefrac=0.3 | |
| 09 | ntrees=2000 | |
| 10 | nuisancetrees=200 | |
| 11 | Combined options | |
| 12 | Null treatment effect | |
| 13 | Strong treatment effect | |
| 14 | High censoring (~60%) | |
| 15 | Low censoring (~10%) | |
| 16 | target = survival.probability | |

---

### Suite 09: multi_arm_causal_forest (36 checks across ~12 configs)

**Scripts:** `run_r_tests.R`, `run_stata_tests.do`

| Config | Arms | Option(s) | Checks |
|--------|:----:|-----------|:------:|
| 01 | 2 | Default | 2 |
| 02 | 3 | Default | 4 |
| 03 | 4 | Default | 6 |
| 04 | 2 | nohonesty | 2 |
| 05 | 2 | mtry=2 | 2 |
| 06 | 2 | cluster() | 2 |
| 07 | 2 | weights() | 2 |
| 08 | 2 | whatinput (K vars) | 2 |
| 09 | 2 | estimatevariance | 4 |
| 10 | 3 | Combined options | 4 |
| 11 | 2 | Homogeneous effect | 2 |
| 12 | 2 | Large p (p=20) | 2 |

---

### Suite 10: boosted_regression_forest (17 tests)

**Scripts:** Individual R/do pairs + `run_all_stata.do`, `compute_correlations.R`

| Test | Option(s) | DGP | Notes |
|------|-----------|-----|-------|
| 01 | Default (auto-tune steps) | Nonlinear | **D-23: step count may differ** |
| 02 | booststeps=1 | same | |
| 03 | booststeps=3 | same | |
| 04 | booststeps=5 | same | |
| 05 | boostmaxsteps=10 | same | |
| 06 | boosterrorreduction=0.90 | same | |
| 07 | boosterrorreduction=0.99 | same | |
| 08 | boosttreestune=50 | same | |
| 09 | nostabilizesplits | same | |
| 10 | cluster() | same | |
| 11 | weights() | same | |
| 12 | nohonesty | same | |
| 13 | mtry=2 | same | |
| 14 | BRF vs plain RF comparison | same | MSE comparison |
| 15 | Linear DGP | Y=X1+X2+eps | |
| 16 | Highly nonlinear DGP | Y=sin(3X1)*cos(2X2) | Marginal pass |

---

### Suite 11: ll_regression_forest (18 tests)

**Scripts:** `run_all_r.R`, `run_all_stata.do`, `compute_correlations.R`

| Test | Option(s) | DGP | Notes |
|------|-----------|-----|-------|
| 01 | Default LL (all vars, lambda=0.1) | Linear+interaction | |
| 02 | llenable | same | |
| 03 | llvars(x1 x2) subset | same | |
| 04 | llsplitvars(x1) | same | |
| 05 | lllambda=0.01 | same | |
| 06 | lllambda=1.0 | same | |
| 07 | lllambda=10.0 | same | |
| 08 | llweightpenalty | same | |
| 09 | llcutoff=3 | same | |
| 10 | cluster() | same | |
| 11 | weights() | same | **D-10: Stata extension** |
| 12 | nohonesty | same | |
| 13 | mtry=2 | same | r=0.9990 (highest) |
| 14 | minnodesize=20 | same | |
| 15 | Combined: llvars+lambda+weightpenalty | same | |
| 16 | LL vs non-LL MSE comparison | same | 52% MSE reduction |
| 17 | Nonlinear DGP | Y=sin(X1)+X2^2 | |
| 18 | Pure linear DGP | Y=sum(X) | |

---

### Suite 12: lm_forest (17 tests)

**Scripts:** `01_gen_data_and_r_tests.R`, `02_stata_tests.do`, `03_compare_and_report.R`

| Test | Option(s) | Notes |
|------|-----------|-------|
| 01 | K=1 (single W) | |
| 02 | K=2 (two W) | |
| 03 | K=3 (three W) | |
| 04 | gradient.weights c(0.7, 0.3) | |
| 05 | stabilize.splits = TRUE | |
| 06 | User-supplied Y.hat | |
| 07 | User-supplied W.hat | |
| 08 | Both Y.hat and W.hat | |
| 09 | cluster() | |
| 10 | weights() | |
| 11 | nohonesty | |
| 12 | mtry=2 | |
| 13 | minnodesize=20 | |
| 14 | nuisancetrees=100 | |
| 15 | estimatevariance + cigroupsize=2 | |
| 16 | Homogeneous beta=2 | |
| 17 | Strong heterogeneity beta=5*I(X1>0) | |

---

### Suite 13: multi_regression_forest (14 tests)

**Scripts:** `run_all_r.R`, `run_all_stata.do`, `compute_correlations.R`

| Test | Option(s) | Notes |
|------|-----------|-------|
| 01 | K=2 outcomes | |
| 02 | K=3 outcomes | |
| 03 | K=5 outcomes | |
| 04 | Correlated outcomes (r~0.98) | |
| 05 | Independent outcomes | |
| 06 | cluster() | |
| 07 | weights() | |
| 08 | nohonesty | |
| 09 | mtry=2 | r=1.0000 (exact match) |
| 10 | minnodesize=20 | |
| 11 | samplefrac=0.3 | |
| 12 | Combined: cluster+weights+nohonesty | |
| 13 | Linear outcomes (low noise) | |
| 14 | Nonlinear outcomes (sin/cos) | |

---

### Suite 14: BLP, test_calibration, get_scores (29 checks)

**Scripts:** `r_reference.R`, `stata_blp.do`

| Test | Command | Option(s) | Notes |
|------|---------|-----------|-------|
| 1-2 | BLP | HC3 (default), all covariates | Coefs + SEs |
| 3-4 | BLP | HC0 | |
| 5-6 | BLP | HC1 | |
| 7-8 | BLP | HC2 | |
| 9-10 | BLP | Subset covariates (X1, X2) | |
| 11-12 | BLP | target.sample=overlap | |
| 13-14 | BLP | clusters | |
| 15-16 | BLP | debiasing.weights | **D-2: formula mismatch** |
| 17 | BLP | target.sample=treated | Stata extension |
| 18 | BLP | target.sample=control | Stata extension |
| 19-20 | BLP | z-test + SE ratio summary | |
| 21-23 | test_calibration | Basic / Strong het / No het | t-stats |
| 24-26 | test_calibration | Same 3 configs — p-values | |
| 27 | get_scores | DR score summary stats | |
| 28 | get_scores | DR score correlation (r=0.9999) | |
| 29 | get_scores | instrumental_forest | N/A |
| 30 | get_scores | mean(DR) ~ ATE check | |

---

### Suite 15: variable_importance, RATE, tune (24 tests)

**Scripts:** `test15_r.R`, `test15_stata.do`, `compare15.R`

| Test | Command | Condition | Notes |
|------|---------|-----------|-------|
| 1 | VI | Default (decay=2, depth=4) | **D-6: re-trains forest** |
| 2 | VI | decay=1.0 | |
| 3 | VI | decay=3.0 | |
| 4 | VI | decay=0.5 | |
| 5 | VI | max.depth=2 | |
| 6 | VI | max.depth=8 | |
| 7 | VI | cluster() | Only VI test that passes Spearman |
| 8 | VI | weights() | |
| 9 | VI | ntrees=1000 | |
| 10 | VI | Ranking correctness | Top-2 always correct |
| 11 | VI | p=20 irrelevant vars | |
| 12 | RATE | AUTOC, CATE priorities | |
| 13 | RATE | QINI, CATE priorities | Near-perfect |
| 14 | RATE | AUTOC, X1 priorities | |
| 15 | RATE | AUTOC, random priorities | Null result |
| 16 | RATE | AUTOC, bootstrap=500 | |
| 17 | RATE | AUTOC, debiasing.weights | |
| 18 | Tune | mtry + minnodesize (regression) | |
| 19 | Tune | samplefrac | |
| 20 | Tune | All parameters | |
| 21 | Tune | causal_forest | |
| 22 | Tune | tunenumtrees=100 | |
| 23 | Tune | tunenumreps=10 | |
| 24 | Tune | MSE improvement check | |

---

## Pass Criteria by Suite

| Suite | Primary Metric | Threshold | Notes |
|-------|---------------|-----------|-------|
| 01 | Pearson r (predictions) | > 0.90 (> 0.95 = STRONG) | Variance: > 0.85 |
| 02 | Pearson r (CATE) | > 0.90 | Variance: > 0.85 |
| 03 | z-test: \|ATE_R - ATE_S\| / SE | < 3.0 | SE ratio < 0.50 |
| 04 | Pearson r (LATE) | > 0.90 | Weak IV: > 0.50 |
| 05 | Pearson r (quantile preds) | > 0.90 | |
| 06 | Pearson r (class probs) | > 0.90 | + prob sum = 1 |
| 07 | Pearson r per column | > 0.85 per column | + monotonicity |
| 08 | Pearson r (CATE on survival) | > 0.85 | surv prob: > 0.75 |
| 09 | Pearson r per arm | > 0.85 | Variance: > 0.80 |
| 10 | Pearson r (BRF predictions) | > 0.90 | |
| 11 | Pearson r (LL predictions) | > 0.90 | Partial: > 0.80 |
| 12 | Pearson r per coefficient | > 0.85 | Variance: > 0.80 |
| 13 | Pearson r per outcome | > 0.90 | |
| 14 | z-test for BLP coefs | < 3.0 | SE ratio [0.5, 2.0] |
| 15 | Spearman rho (VI) | > 0.70 | RATE: \|diff\|/SE < 3 |

---

## Adding New Tests

To add a test to an existing suite:

1. **R script:** Add the test configuration, generate R predictions, append to the
   suite's data CSV with a column like `r_pred_testNN` or write a new `testNN_data.csv`
2. **Stata do-file:** Load the CSV, run the matching Stata command with the same
   options, export Stata predictions to `testNN_stata.csv`
3. **Comparison:** Add the new test to the suite's comparison R script
4. **Report:** Update the suite's markdown report with results

To add a new suite:

1. Create a directory `NN_command_name/`
2. Follow the 3-step pipeline pattern (R → Stata → Compare)
3. Add the suite to `run_all_suites.sh` SUITES array
4. Write a markdown report `NN_command_name.md`
5. Update `CONSOLIDATED_DISCREPANCIES.md` with any new findings

---

## Data Generating Processes

Most suites share a base DGP unless testing specific scenarios:

| Suite | Default DGP | n | p |
|-------|-------------|---|---|
| 01 | Y = X1 + 2*X2 + eps; X ~ U(0,1) | 500 | 5 |
| 02-03 | Y = X1 + (X1+X2)*W + eps; X ~ N(0,1); W ~ Bern(0.5) | 500 | 5 |
| 04 | Y = X1 + (X1+X2)*W + 0.5*X3 + eps; Z ~ Bern(0.5) | 500 | 5 |
| 05 | Y = X1 + 2*X2 + eps; X ~ N(0,1) | 500 | 5 |
| 06 | Y ~ Multinomial(softmax(X beta)) | 500 | 5 |
| 07 | T ~ Exp(exp(0.5*X1)); C ~ Exp(0.3) | 500 | 5 |
| 08 | T ~ Exp with treatment effect on hazard | 500 | 5 |
| 09 | Y = X1 + tau_k*W_k + eps; multi-arm | 500 | 5 |
| 10 | Y = sin(2X1) + X2^2 + 0.5X3 + eps | 500 | 5 |
| 11 | Y = 3X1 + 2X2 - X3 + 0.5*X1*X2 + eps | 500 | 5 |
| 12 | Y = X1 + (X1+X2)*W1 + (-X1+0.5X3)*W2 + eps | 500 | 5 |
| 13 | Y1 = X1+X2+eps; Y2 = -X1+X3+eps; Y3 = X2*X3+eps | 500 | 5 |
| 14 | Y = X1 + (X1+X2)*W + eps; n=1000 | 1000 | 5 |
| 15 | Y = 3X1 + 2X2 + eps; X ~ N(0,1) | 500 | 10 |

All DGPs use `set.seed(42)` / `seed(42)` for reproducibility.

---

## Interpreting Results

**Pearson r > 0.95**: Excellent — implementations are functionally identical.
Minor differences arise from PRNG divergence between R and Stata plugin wrappers.

**Pearson r 0.85-0.95**: Good — implementations agree on signal but differ in
noise. Common with options that increase forest variance (samplefrac=0.3, mtry=2)
or when nuisance estimation paths differ.

**Pearson r < 0.85**: Investigate. Check if:
- The predictions have near-zero variance (narrow range → low correlation artifact)
- An API design difference causes column misalignment (e.g., numfailures)
- A genuine bug exists (e.g., cluster not propagated)

**Variance correlation < 0.30**: Expected. Variance estimation uses jackknife
groups that are assigned by a different RNG path in R vs Stata. Per-observation
variance correlations are inherently low even within R across seeds. Focus on
distributional agreement (mean ratio within 20%) instead.
