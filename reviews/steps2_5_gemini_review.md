Loaded cached credentials.
Hook registry initialized with 0 hook entries
Here is the review of Steps 2-5 of the gap-closing plan. While the architecture and Stata syntax correctly expose the new functionality, there are several critical bugs in how the C++ plugin utilizes the data and how the `.ado` files construct the nuisance pipelines. 

Here are the specific issues that need to be addressed:

### 1. `sample_weights` are entirely ignored in the C++ backend
**File:** `grf_plugin.cpp` (Lines ~337-353)
You successfully parse `weight_col_idx` and populate the `std::vector<double> sample_weights`. However, this vector is **never passed** to the `grf::Data` object, `ForestOptions`, or the `trainer`. Consequently, the forest executes completely unweighted. *(Note: `test_weights.do` likely passes due to multithreading non-determinism altering the random split sequences, rather than actual weight applications).*
*   **Fix:** Apply the weights to the data object (e.g., `data.set_sample_weights(sample_weights);` or `data.set_weight_index(...)` depending on your `grf` C++ library version). 

### 2. Nuisance models do not respect clustering or sample weights
**Files:** `grf_causal_forest.ado`, `grf_instrumental_forest.ado`, `grf_multi_arm_causal_forest.ado`, `grf_causal_survival_forest.ado`, `grf_lm_forest.ado` (Various plugin calls)
In all double-machine-learning wrappers, the first-stage regression forests (e.g., $Y \sim X$ and $W \sim X$) omit the `extra_vars` from the `plugin call` variable list and hardcode `"0"` for both `cluster_col_idx` and `weight_col_idx`. 
*   **Fix:** Nuisance estimation must use the exact same sample weights and clustering as the main forest for orthogonalization to be valid. You must pass `` `extra_vars' ``, `` `cluster_col_idx' ``, and `` `weight_col_idx' `` to the nuisance regression plugin calls.

### 3. Missing `noMIA`, `CLuster()`, and `WEIghts()` in Tuning 
**File:** `grf_tune.ado` (Lines ~16-33, and all `plugin call` blocks)
*   The `syntax` block completely omits the `noMIA`, `CLuster()`, and `WEIghts()` options.
*   In every `plugin call` within the random search loop, `allow_missing_x` is hardcoded to `"1"`, and cluster/weight indices are hardcoded to `"0"`.

### 4. `grf_predict.ado` unconditionally uses MIA
**File:** `grf_predict.ado` (Lines ~335, 381, 417, etc.)
`grf_predict` passes `"1"` for `allow_missing_x` to the plugin regardless of how the forest was trained. If a user trained a forest with `nomia`, the test data predictions will silently use MIA instead of casewise deletion.
*   **Fix:** The training scripts should save `ereturn local mia = "\`allow_missing_x'"` so `grf_predict` can inherit the boolean natively, or `grf_predict` should accept its own `noMIA` syntax.

### 5. Incorrect DR score math for Multi-Arm Causal
**File:** `grf_get_scores.ado` (Lines ~69-72)
For `multi_causal`, the $Y$ residual is computed iteratively for each arm $j$ as:
```stata
quietly gen double `y_resid_`j'' = `depvar' - `yhatvar' - `tau_j' * `w_resid_`j''
```
This is mathematically incorrect for a joint linear model ($Y = c(X) + \sum \tau_k W_k$). You must subtract the treatment effects of **all** arms simultaneously to isolate the residual: $Y - \hat{Y} - \sum_{k=1}^{K} \hat{\tau}_k (W_k - \hat{W}_k)$. This summed residual should be computed once and reused for all arm scores.

### 6. Unfinished DR score for Causal Survival
**File:** `grf_get_scores.ado` (Lines ~144-155)
The code calculates `w_resid` and `w_resid_var`, but then leaves them completely unused. It simply sets `generate` equal to `tauvar`. 
*   **Fix:** Implement the propensity correction (`tauvar + (w_resid / w_resid_var) * ...`) if you intend to offer a simplified DR score, or remove the dead variable generations if standard output is intended.
