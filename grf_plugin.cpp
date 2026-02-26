/*
 * grf_plugin.cpp -- Stata plugin wrapping the grf C++ library
 *
 * Generalized Random Forests (Athey, Tibshirani, Wager 2019)
 * Wraps the grf C++ backend for all forest types:
 *   regression, causal, quantile, instrumental, probability,
 *   survival, causal_survival, multi_arm_causal, multi_regression,
 *   ll_regression, boosted_regression, lm_forest
 *
 * Dispatched by forest type via argv[0]:
 *   "regression", "causal", "quantile", "instrumental", "probability",
 *   "survival", "causal_survival", "multi_arm_causal", "multi_regression",
 *   "ll_regression", "boosted_regression", "lm_forest"
 *
 * Copyright: GPL-3.0 (following grf license)
 */

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>
#include <numeric>
#include <fstream>
#include <sstream>
#include <set>

/* Eigen linear algebra (for local linear regression) */
#include <Eigen/Dense>

/* Stata plugin interface -- must be included with C linkage */
extern "C" {
#include "stplugin.h"
}

/* grf C++ library headers */
#include "commons/Data.h"
#include "commons/globals.h"
#include "forest/Forest.h"
#include "forest/ForestOptions.h"
#include "forest/ForestTrainer.h"
#include "forest/ForestTrainers.h"
#include "forest/ForestPredictor.h"
#include "forest/ForestPredictors.h"
#include "prediction/Prediction.h"
#include "analysis/SplitFrequencyComputer.h"
#include "tree/Tree.h"

/* ================================================================
 * Helper: parse integer from argv with default
 * ================================================================ */
static int parse_int(const char* s, int def) {
    if (!s || !*s) return def;
    int v = atoi(s);
    return (v != 0 || strcmp(s, "0") == 0) ? v : def;
}

static double parse_double(const char* s, double def) {
    if (!s || !*s) return def;
    double v = atof(s);
    return v;
}

/* ================================================================
 * Helper: split column-major data into train and test arrays.
 *
 * train_vec: first n_train rows, all n_cols columns
 * test_vec:  remaining rows, only first n_x columns (X only)
 *
 * predict() only needs X columns in test_data; the prediction
 * strategies read Y/W/Z from train_data, not test_data.
 * ================================================================ */
static void split_train_test(
    const std::vector<double>& all_data,
    int n_total, int n_cols, int n_train, int n_x,
    std::vector<double>& train_vec,
    std::vector<double>& test_vec)
{
    int n_test = n_total - n_train;

    train_vec.resize((size_t)n_cols * n_train);
    for (int col = 0; col < n_cols; col++)
        for (int row = 0; row < n_train; row++)
            train_vec[(size_t)col * n_train + row] =
                all_data[(size_t)col * n_total + row];

    test_vec.resize((size_t)n_x * n_test);
    for (int col = 0; col < n_x; col++)
        for (int row = 0; row < n_test; row++)
            test_vec[(size_t)col * n_test + row] =
                all_data[(size_t)col * n_total + (n_train + row)];
}

/* Helper: set standard indices on a grf::Data object */
static void set_data_indices(grf::Data& d, int y_start, int n_y,
                             int w_start, int n_w, int z_start, int n_z)
{
    if (n_y == 1) {
        d.set_outcome_index((size_t)y_start);
    } else {
        std::vector<size_t> idx(n_y);
        for (int i = 0; i < n_y; i++) idx[i] = (size_t)(y_start + i);
        d.set_outcome_index(idx);
    }
    if (n_w > 0) {
        if (n_w == 1) {
            d.set_treatment_index((size_t)w_start);
        } else {
            std::vector<size_t> idx(n_w);
            for (int i = 0; i < n_w; i++) idx[i] = (size_t)(w_start + i);
            d.set_treatment_index(idx);
        }
    }
    if (n_z > 0) {
        d.set_instrument_index((size_t)z_start);
    }
}

/* ================================================================
 * Main entry point
 * ================================================================
 *
 * argv[0] = forest_type: "regression", "causal", "quantile",
 *           "instrumental", "probability", "survival",
 *           "causal_survival", "multi_arm_causal", "multi_regression",
 *           "ll_regression", "boosted_regression", "lm_forest",
 *           "variable_importance", "split_frequencies"
 *
 * Common args (argv[1..]):
 *   [1]  num_trees (int, default 2000)
 *   [2]  seed (int, default 42)
 *   [3]  mtry (int, 0=auto sqrt(p))
 *   [4]  min_node_size (int, 5)
 *   [5]  sample_fraction (double, 0.5)
 *   [6]  honesty (int: 0/1, default 1)
 *   [7]  honesty_fraction (double, 0.5)
 *   [8]  honesty_prune_leaves (int: 0/1, default 1)
 *   [9]  alpha (double, 0.05)
 *   [10] imbalance_penalty (double, 0.0)
 *   [11] ci_group_size (int, default 2 for variance, 1 otherwise)
 *   [12] num_threads (int, 0=auto)
 *   [13] estimate_variance (int: 0/1, default 0)
 *   [14] n_train (int: 0=OOB on all data, >0=first n_train rows are train, rest are test)
 *   [15] n_x (int, number of X variables)
 *   [16] n_y (int, number of outcome variables)
 *   [17] n_w (int, number of treatment variables)
 *   [18] n_z (int, number of instrument variables)
 *   [19] n_output (int, number of output variables)
 *
 * Forest-specific args (argv[20+]):
 *   Causal: [20]=stabilize_splits (int: 0/1)
 *   Quantile: [20]=quantiles (comma-separated, e.g. "0.1,0.5,0.9")
 *   Instrumental: [20]=reduced_form_weight (double), [21]=stabilize_splits
 *   Probability: [20]=num_classes (int)
 *   Survival: [20]=num_failures (int), [21]=prediction_type (int)
 *   Causal Survival: [20]=stabilize_splits (int), [21-23]=col indices,
 *                    [24]=target (int: 1=RMST, 2=survival probability)
 *   Multi-arm Causal: [20]=stabilize_splits (int)
 *   LL Regression: [20]=enable_ll_split (int), [21]=ll_lambda (double),
 *                  [22]=ll_weight_penalty (int), [23]=ll_split_cutoff (int)
 *   Boosted Regression: [20]=boost_steps (int), [21]=boost_error_reduction (double),
 *                       [22]=boost_max_steps (int), [23]=boost_trees_tune (int),
 *                       [24]=stabilize_splits (int: 0/1)
 *   LM Forest: [20]=stabilize_splits (int)
 *   Variable Importance: [20]=max_depth (int)
 *   Split Frequencies: [20]=max_depth (int)
 */
extern "C" STDLL stata_call(int argc, char *argv[])
{
    char msg[1024];

    if (argc < 20) {
        snprintf(msg, sizeof(msg),
                 "GRF error: expected at least 20 arguments, got %d\n", argc);
        SF_error(msg);
        return 198;
    }

    /* ----------------------------------------------------------
     * Step 0: Parse common arguments
     * ---------------------------------------------------------- */
    std::string forest_type = argv[0];

    int num_trees           = parse_int(argv[1], 2000);
    int seed                = parse_int(argv[2], 42);
    int mtry                = parse_int(argv[3], 0);
    int min_node_size       = parse_int(argv[4], 5);
    double sample_fraction  = parse_double(argv[5], 0.5);
    int honesty             = parse_int(argv[6], 1);
    double honesty_fraction = parse_double(argv[7], 0.5);
    int honesty_prune       = parse_int(argv[8], 1);
    double alpha            = parse_double(argv[9], 0.05);
    double imbalance_pen    = parse_double(argv[10], 0.0);
    int ci_group_size       = parse_int(argv[11], 1);
    int num_threads         = parse_int(argv[12], 0);
    int estimate_variance   = parse_int(argv[13], 0);
    int n_train             = parse_int(argv[14], 0);
    int n_x                 = parse_int(argv[15], 0);
    int n_y                 = parse_int(argv[16], 1);
    int n_w                 = parse_int(argv[17], 0);
    int n_z                 = parse_int(argv[18], 0);
    int n_output            = parse_int(argv[19], 1);

    /* Validate */
    if (num_trees <= 0) num_trees = 2000;
    if (min_node_size <= 0) min_node_size = 5;
    if (sample_fraction <= 0.0 || sample_fraction > 1.0) sample_fraction = 0.5;
    if (alpha < 0.0 || alpha > 0.25) alpha = 0.05;
    if (honesty_fraction <= 0.0 || honesty_fraction >= 1.0) honesty_fraction = 0.5;
    if (ci_group_size <= 0) ci_group_size = 1;
    if (estimate_variance && ci_group_size < 2) ci_group_size = 2;
    if (seed <= 0) seed = 42;
    bool predict_mode = (n_train > 0);

    /* ----------------------------------------------------------
     * Step 1: Read data from Stata
     * ---------------------------------------------------------- */
    int nvar = SF_nvars();
    int n_data_cols = nvar - n_output;

    if (n_x <= 0) {
        /* Auto-detect: n_x = total data cols - n_y - n_w - n_z */
        n_x = n_data_cols - n_y - n_w - n_z;
    }

    if (n_x < 1) {
        SF_error("GRF error: need at least 1 predictor variable.\n");
        return 198;
    }

    /* Count observations */
    ST_int obs1 = SF_in1();
    ST_int obs2 = SF_in2();
    int n = 0;
    for (ST_int i = obs1; i <= obs2; i++) {
        if (SF_ifobs(i)) n++;
    }

    if (n < 2) {
        SF_error("GRF error: need at least 2 non-missing observations.\n");
        return 2000;
    }

    /* Auto mtry */
    if (mtry <= 0) {
        mtry = (int)std::ceil(std::sqrt((double)n_x));
        if (mtry < 1) mtry = 1;
    }
    if (mtry > n_x) mtry = n_x;

    snprintf(msg, sizeof(msg),
             "GRF %s forest: n=%d, p=%d, trees=%d, mtry=%d, "
             "min_node=%d, honesty=%d\n",
             forest_type.c_str(), n, n_x, num_trees, mtry,
             min_node_size, honesty);
    SF_display(msg);

    /* Read data into column-major array
     *
     * Variable ordering in Stata (set by .ado):
     *   X1..Xp  Y1..Yn_y  W1..Wn_w  Z1..Zn_z  (extra cols)  out1..outn_output
     *
     * For grf::Data, we arrange columns as:
     *   X1..Xp  Y1..Yn_y  W1..Wn_w  Z1..Zn_z  (extra cols)
     * (same order minus outputs)
     *
     * Two-pass approach: first pass counts valid obs, second reads data.
     */
    std::vector<int> obs_map;
    obs_map.reserve(n);

    /* Pass 1: identify valid (non-missing) observations */
    for (ST_int i = obs1; i <= obs2; i++) {
        if (!SF_ifobs(i)) continue;

        bool has_missing = false;
        for (int j = 0; j < n_data_cols; j++) {
            double val;
            ST_retcode rc = SF_vdata(j + 1, i, &val);
            if (rc || SF_is_missing(val)) {
                has_missing = true;
                break;
            }
        }
        if (!has_missing) {
            obs_map.push_back((int)i);
        }
    }

    n = (int)obs_map.size();
    if (n < 2) {
        SF_error("GRF error: fewer than 2 complete observations.\n");
        return 2000;
    }

    /* Validate predict mode */
    if (predict_mode) {
        if (n_train >= n) {
            snprintf(msg, sizeof(msg),
                     "GRF error: n_train=%d >= n_obs=%d, no test data.\n",
                     n_train, n);
            SF_error(msg);
            return 198;
        }
        if (n_train < 2) {
            SF_error("GRF error: need at least 2 training observations.\n");
            return 2000;
        }
    }

    /* Pass 2: read data in column-major order */
    std::vector<double> data_vec((size_t)n_data_cols * n);
    for (int idx = 0; idx < n; idx++) {
        ST_int i = obs_map[idx];
        for (int j = 0; j < n_data_cols; j++) {
            double val;
            SF_vdata(j + 1, i, &val);
            data_vec[(size_t)j * n + idx] = val;
        }
    }

    /* ----------------------------------------------------------
     * Step 2: Create grf::Data and set indices
     * ---------------------------------------------------------- */
    grf::Data data(data_vec.data(), (size_t)n, (size_t)n_data_cols);

    /* Column indices: X is 0..n_x-1, Y starts at n_x, W at n_x+n_y, Z at n_x+n_y+n_w */
    int y_start = n_x;
    int w_start = n_x + n_y;
    int z_start = n_x + n_y + n_w;

    /* Set outcome index(es) */
    if (n_y == 1) {
        data.set_outcome_index((size_t)y_start);
    } else {
        std::vector<size_t> outcome_idx(n_y);
        for (int i = 0; i < n_y; i++) outcome_idx[i] = (size_t)(y_start + i);
        data.set_outcome_index(outcome_idx);
    }

    /* Set treatment index(es) if applicable */
    if (n_w > 0) {
        if (n_w == 1) {
            data.set_treatment_index((size_t)w_start);
        } else {
            std::vector<size_t> treat_idx(n_w);
            for (int i = 0; i < n_w; i++) treat_idx[i] = (size_t)(w_start + i);
            data.set_treatment_index(treat_idx);
        }
    }

    /* Set instrument index if applicable */
    if (n_z > 0) {
        data.set_instrument_index((size_t)z_start);
    }

    /* For survival forests, set censor index
     * Layout: X Y(time) censor W (treatment for causal survival)
     * We handle special indexing per forest type below */

    /* ----------------------------------------------------------
     * Step 3: Create ForestOptions
     * ---------------------------------------------------------- */
    std::vector<size_t> clusters;  /* empty = no clustering */
    grf::uint samples_per_cluster = 0;
    bool legacy_seed = false;

    grf::ForestOptions options(
        (grf::uint)num_trees,
        (size_t)ci_group_size,
        sample_fraction,
        (grf::uint)mtry,
        (grf::uint)min_node_size,
        (honesty != 0),
        honesty_fraction,
        (honesty_prune != 0),
        alpha,
        imbalance_pen,
        (grf::uint)num_threads,
        (grf::uint)seed,
        legacy_seed,
        clusters,
        samples_per_cluster
    );

    /* ----------------------------------------------------------
     * Step 4: Create trainer and predictor based on forest type,
     *         train, and predict
     * ---------------------------------------------------------- */
    std::vector<grf::Prediction> predictions;
    grf::uint resolved_threads = grf::ForestOptions::validate_num_threads((grf::uint)num_threads);

    try {

    if (forest_type == "regression") {
        /* ---- Regression Forest ---- */
        grf::ForestTrainer trainer = grf::regression_trainer();
        grf::ForestPredictor predictor = grf::regression_predictor(resolved_threads);
        bool est_var = (estimate_variance != 0);
        int out_col_pred = nvar - n_output + 1;  /* 1-indexed Stata var */
        int out_col_var = (n_output >= 2 && est_var) ? out_col_pred + 1 : 0;
        int n_written = 0;

        if (predict_mode) {
            int n_test = n - n_train;
            std::vector<double> train_vec, test_vec;
            split_train_test(data_vec, n, n_data_cols, n_train, n_x, train_vec, test_vec);
            grf::Data train_data(train_vec.data(), (size_t)n_train, (size_t)n_data_cols);
            set_data_indices(train_data, y_start, n_y, w_start, n_w, z_start, n_z);
            grf::Data test_data(test_vec.data(), (size_t)n_test, (size_t)n_x);

            SF_display("  Training regression forest...\n");
            grf::Forest forest = trainer.train(train_data, options);
            SF_display("  Forest trained. Predicting on new data...\n");
            predictions = predictor.predict(forest, train_data, test_data, est_var);

            for (int i = 0; i < n_test; i++) {
                const auto& pred = predictions[i].get_predictions();
                if (!pred.empty() && std::isfinite(pred[0])) {
                    SF_vstore(out_col_pred, obs_map[n_train + i], pred[0]);
                    n_written++;
                }
                if (out_col_var > 0 && predictions[i].contains_variance_estimates()) {
                    const auto& var_est = predictions[i].get_variance_estimates();
                    if (!var_est.empty()) {
                        SF_vstore(out_col_var, obs_map[n_train + i], var_est[0]);
                    }
                }
            }
        } else {
            SF_display("  Training regression forest...\n");
            grf::Forest forest = trainer.train(data, options);
            SF_display("  Forest trained.\n");

            SF_display("  Computing predictions...\n");
            predictions = predictor.predict_oob(forest, data, est_var);

            for (int i = 0; i < n; i++) {
                const auto& pred = predictions[i].get_predictions();
                if (!pred.empty() && std::isfinite(pred[0])) {
                    SF_vstore(out_col_pred, obs_map[i], pred[0]);
                    n_written++;
                }
                if (out_col_var > 0 && predictions[i].contains_variance_estimates()) {
                    const auto& var_est = predictions[i].get_variance_estimates();
                    if (!var_est.empty()) {
                        SF_vstore(out_col_var, obs_map[i], var_est[0]);
                    }
                }
            }
        }

        snprintf(msg, sizeof(msg), "  Wrote %d predictions.\n", n_written);
        SF_display(msg);

    } else if (forest_type == "causal") {
        /* ---- Causal Forest ----
         * Uses multi_causal_trainer(num_treatments=1, num_outcomes=1, stabilize_splits)
         * Data layout: X Y.centered W.centered
         * The .ado handles nuisance estimation (Y.hat, W.hat) and centering.
         */
        int stabilize = (argc > 20) ? parse_int(argv[20], 1) : 1;

        grf::ForestTrainer trainer = grf::multi_causal_trainer(
            (size_t)n_w, (size_t)n_y, (stabilize != 0));
        grf::ForestPredictor predictor = grf::multi_causal_predictor(
            resolved_threads, (size_t)n_w, (size_t)n_y);
        bool est_var = (estimate_variance != 0);
        int out_col_pred = nvar - n_output + 1;
        int out_col_var = (n_output >= 2 && est_var) ? out_col_pred + 1 : 0;
        int n_written = 0;

        if (predict_mode) {
            int n_test = n - n_train;
            std::vector<double> train_vec, test_vec;
            split_train_test(data_vec, n, n_data_cols, n_train, n_x, train_vec, test_vec);
            grf::Data train_data(train_vec.data(), (size_t)n_train, (size_t)n_data_cols);
            set_data_indices(train_data, y_start, n_y, w_start, n_w, z_start, n_z);
            grf::Data test_data(test_vec.data(), (size_t)n_test, (size_t)n_x);

            SF_display("  Training causal forest...\n");
            grf::Forest forest = trainer.train(train_data, options);
            SF_display("  Forest trained. Predicting on new data...\n");
            predictions = predictor.predict(forest, train_data, test_data, est_var);

            for (int i = 0; i < n_test; i++) {
                const auto& pred = predictions[i].get_predictions();
                if (!pred.empty() && std::isfinite(pred[0])) {
                    SF_vstore(out_col_pred, obs_map[n_train + i], pred[0]);
                    n_written++;
                }
                if (out_col_var > 0 && predictions[i].contains_variance_estimates()) {
                    const auto& var_est = predictions[i].get_variance_estimates();
                    if (!var_est.empty()) {
                        SF_vstore(out_col_var, obs_map[n_train + i], var_est[0]);
                    }
                }
            }
        } else {
            SF_display("  Training causal forest...\n");
            grf::Forest forest = trainer.train(data, options);
            SF_display("  Forest trained.\n");

            SF_display("  Computing predictions...\n");
            predictions = predictor.predict_oob(forest, data, est_var);

            for (int i = 0; i < n; i++) {
                const auto& pred = predictions[i].get_predictions();
                if (!pred.empty() && std::isfinite(pred[0])) {
                    SF_vstore(out_col_pred, obs_map[i], pred[0]);
                    n_written++;
                }
                if (out_col_var > 0 && predictions[i].contains_variance_estimates()) {
                    const auto& var_est = predictions[i].get_variance_estimates();
                    if (!var_est.empty()) {
                        SF_vstore(out_col_var, obs_map[i], var_est[0]);
                    }
                }
            }
        }

        snprintf(msg, sizeof(msg), "  Wrote %d predictions (CATE estimates).\n", n_written);
        SF_display(msg);

    } else if (forest_type == "quantile") {
        /* ---- Quantile Forest ----
         * argv[20] = quantiles (comma-separated, e.g. "0.1,0.5,0.9")
         */
        std::vector<double> quantiles;
        if (argc > 20 && argv[20] && *argv[20]) {
            std::string qstr(argv[20]);
            std::stringstream ss(qstr);
            std::string token;
            while (std::getline(ss, token, ',')) {
                double q = atof(token.c_str());
                if (q > 0.0 && q < 1.0) {
                    quantiles.push_back(q);
                }
            }
        }
        if (quantiles.empty()) {
            quantiles = {0.1, 0.5, 0.9};
        }

        grf::ForestTrainer trainer = grf::quantile_trainer(quantiles);
        grf::ForestPredictor predictor = grf::quantile_predictor(resolved_threads, quantiles);
        int out_col_start = nvar - n_output + 1;
        int n_written = 0;

        if (predict_mode) {
            int n_test = n - n_train;
            std::vector<double> train_vec, test_vec;
            split_train_test(data_vec, n, n_data_cols, n_train, n_x, train_vec, test_vec);
            grf::Data train_data(train_vec.data(), (size_t)n_train, (size_t)n_data_cols);
            set_data_indices(train_data, y_start, n_y, w_start, n_w, z_start, n_z);
            grf::Data test_data(test_vec.data(), (size_t)n_test, (size_t)n_x);

            snprintf(msg, sizeof(msg), "  Training quantile forest (%zu quantiles)...\n",
                     quantiles.size());
            SF_display(msg);
            grf::Forest forest = trainer.train(train_data, options);
            SF_display("  Forest trained. Predicting on new data...\n");
            predictions = predictor.predict(forest, train_data, test_data, false);

            for (int i = 0; i < n_test; i++) {
                const auto& pred = predictions[i].get_predictions();
                for (size_t q = 0; q < quantiles.size() && q < (size_t)n_output; q++) {
                    if (q < pred.size() && std::isfinite(pred[q])) {
                        SF_vstore(out_col_start + (int)q, obs_map[n_train + i], pred[q]);
                    }
                }
                if (!pred.empty()) n_written++;
            }
        } else {
            snprintf(msg, sizeof(msg), "  Training quantile forest (%zu quantiles)...\n",
                     quantiles.size());
            SF_display(msg);
            grf::Forest forest = trainer.train(data, options);
            SF_display("  Forest trained.\n");

            SF_display("  Computing predictions...\n");
            predictions = predictor.predict_oob(forest, data, false);

            for (int i = 0; i < n; i++) {
                const auto& pred = predictions[i].get_predictions();
                for (size_t q = 0; q < quantiles.size() && q < (size_t)n_output; q++) {
                    if (q < pred.size() && std::isfinite(pred[q])) {
                        SF_vstore(out_col_start + (int)q, obs_map[i], pred[q]);
                    }
                }
                if (!pred.empty()) n_written++;
            }
        }

        snprintf(msg, sizeof(msg), "  Wrote %d predictions (%zu quantiles per obs).\n",
                 n_written, quantiles.size());
        SF_display(msg);

    } else if (forest_type == "instrumental") {
        /* ---- Instrumental Forest ----
         * argv[20] = reduced_form_weight (double, default 0)
         * argv[21] = stabilize_splits (int, default 1)
         */
        double reduced_form_weight = (argc > 20) ? parse_double(argv[20], 0.0) : 0.0;
        int stabilize = (argc > 21) ? parse_int(argv[21], 1) : 1;

        grf::ForestTrainer trainer = grf::instrumental_trainer(
            reduced_form_weight, (stabilize != 0));
        grf::ForestPredictor predictor = grf::instrumental_predictor(resolved_threads);
        bool est_var = (estimate_variance != 0);
        int out_col_pred = nvar - n_output + 1;
        int out_col_var = (n_output >= 2 && est_var) ? out_col_pred + 1 : 0;
        int n_written = 0;

        if (predict_mode) {
            int n_test = n - n_train;
            std::vector<double> train_vec, test_vec;
            split_train_test(data_vec, n, n_data_cols, n_train, n_x, train_vec, test_vec);
            grf::Data train_data(train_vec.data(), (size_t)n_train, (size_t)n_data_cols);
            set_data_indices(train_data, y_start, n_y, w_start, n_w, z_start, n_z);
            grf::Data test_data(test_vec.data(), (size_t)n_test, (size_t)n_x);

            SF_display("  Training instrumental forest...\n");
            grf::Forest forest = trainer.train(train_data, options);
            SF_display("  Forest trained. Predicting on new data...\n");
            predictions = predictor.predict(forest, train_data, test_data, est_var);

            for (int i = 0; i < n_test; i++) {
                const auto& pred = predictions[i].get_predictions();
                if (!pred.empty() && std::isfinite(pred[0])) {
                    SF_vstore(out_col_pred, obs_map[n_train + i], pred[0]);
                    n_written++;
                }
                if (out_col_var > 0 && predictions[i].contains_variance_estimates()) {
                    const auto& var_est = predictions[i].get_variance_estimates();
                    if (!var_est.empty()) {
                        SF_vstore(out_col_var, obs_map[n_train + i], var_est[0]);
                    }
                }
            }
        } else {
            SF_display("  Training instrumental forest...\n");
            grf::Forest forest = trainer.train(data, options);
            SF_display("  Forest trained.\n");

            SF_display("  Computing predictions...\n");
            predictions = predictor.predict_oob(forest, data, est_var);

            for (int i = 0; i < n; i++) {
                const auto& pred = predictions[i].get_predictions();
                if (!pred.empty() && std::isfinite(pred[0])) {
                    SF_vstore(out_col_pred, obs_map[i], pred[0]);
                    n_written++;
                }
                if (out_col_var > 0 && predictions[i].contains_variance_estimates()) {
                    const auto& var_est = predictions[i].get_variance_estimates();
                    if (!var_est.empty()) {
                        SF_vstore(out_col_var, obs_map[i], var_est[0]);
                    }
                }
            }
        }

        snprintf(msg, sizeof(msg), "  Wrote %d predictions (LATE estimates).\n", n_written);
        SF_display(msg);

    } else if (forest_type == "probability") {
        /* ---- Probability (Classification) Forest ----
         * argv[20] = num_classes (int)
         */
        int num_classes = (argc > 20) ? parse_int(argv[20], 2) : 2;
        if (num_classes < 2) num_classes = 2;

        grf::ForestTrainer trainer = grf::probability_trainer((size_t)num_classes);
        grf::ForestPredictor predictor = grf::probability_predictor(
            resolved_threads, (size_t)num_classes);
        int out_col_start = nvar - n_output + 1;
        int n_written = 0;

        if (predict_mode) {
            int n_test = n - n_train;
            std::vector<double> train_vec, test_vec;
            split_train_test(data_vec, n, n_data_cols, n_train, n_x, train_vec, test_vec);
            grf::Data train_data(train_vec.data(), (size_t)n_train, (size_t)n_data_cols);
            set_data_indices(train_data, y_start, n_y, w_start, n_w, z_start, n_z);
            grf::Data test_data(test_vec.data(), (size_t)n_test, (size_t)n_x);

            snprintf(msg, sizeof(msg), "  Training probability forest (%d classes)...\n", num_classes);
            SF_display(msg);
            grf::Forest forest = trainer.train(train_data, options);
            SF_display("  Forest trained. Predicting on new data...\n");
            predictions = predictor.predict(forest, train_data, test_data, false);

            for (int i = 0; i < n_test; i++) {
                const auto& pred = predictions[i].get_predictions();
                for (int c = 0; c < num_classes && c < n_output; c++) {
                    if ((size_t)c < pred.size() && std::isfinite(pred[c])) {
                        SF_vstore(out_col_start + c, obs_map[n_train + i], pred[c]);
                    }
                }
                if (!pred.empty()) n_written++;
            }
        } else {
            snprintf(msg, sizeof(msg), "  Training probability forest (%d classes)...\n", num_classes);
            SF_display(msg);
            grf::Forest forest = trainer.train(data, options);
            SF_display("  Forest trained.\n");

            SF_display("  Computing predictions...\n");
            predictions = predictor.predict_oob(forest, data, false);

            for (int i = 0; i < n; i++) {
                const auto& pred = predictions[i].get_predictions();
                for (int c = 0; c < num_classes && c < n_output; c++) {
                    if ((size_t)c < pred.size() && std::isfinite(pred[c])) {
                        SF_vstore(out_col_start + c, obs_map[i], pred[c]);
                    }
                }
                if (!pred.empty()) n_written++;
            }
        }

        snprintf(msg, sizeof(msg), "  Wrote %d predictions (%d class probabilities).\n",
                 n_written, num_classes);
        SF_display(msg);

    } else if (forest_type == "survival") {
        /* ---- Survival Forest ----
         * Data layout: X Y(time) censor
         * argv[20] = num_failures (int, auto-detected from data)
         * argv[21] = prediction_type (int, 0=Kaplan-Meier, 1=Nelson-Aalen)
         *
         * The censor variable is the last data column before outputs.
         * We need to set the censor index.
         */
        int num_failures_arg = (argc > 20) ? parse_int(argv[20], 0) : 0;
        int prediction_type = (argc > 21) ? parse_int(argv[21], 0) : 0;
        bool fast_logrank = true;

        /* Set censor index: it's at position n_x + n_y (assuming n_y=1 for time) */
        /* Layout: X(0..n_x-1) time(n_x) censor(n_x+1)
         * But we've already set outcome to n_x, so censor is n_x + 1 */
        int censor_col = y_start + n_y;  /* After outcome(s) */
        if (censor_col < n_data_cols) {
            data.set_censor_index((size_t)censor_col);
        }

        /* Relabel survival times to integer indices.
         * grf C++ expects outcome = findInterval(Y, sort(unique(Y[D==1]))):
         *   - Collect sorted unique failure times (where censor > 0)
         *   - For each obs, replace Y with its interval index
         *   - This maps continuous times to 0, 1, ..., num_failures
         *
         * In predict mode, only TRAINING rows (0..n_train-1) are used to
         * build the failure time set, but ALL rows are relabeled.
         */
        int n_relabel_src = predict_mode ? n_train : n;
        std::vector<double> failure_times_vec;
        {
            std::set<double> ft_set;
            for (int i = 0; i < n_relabel_src; i++) {
                double t = data_vec[(size_t)y_start * n + i];
                double c = data_vec[(size_t)censor_col * n + i];
                if (c > 0.0) {
                    ft_set.insert(t);
                }
            }
            failure_times_vec.assign(ft_set.begin(), ft_set.end());
        }
        int num_failures = (int)failure_times_vec.size();
        if (num_failures <= 0) num_failures = 1;
        /* If user requested fewer failure times, subsample evenly */
        if (num_failures_arg > 0 && num_failures_arg < num_failures) {
            std::vector<double> subsampled;
            subsampled.reserve(num_failures_arg);
            for (int k = 0; k < num_failures_arg; k++) {
                int idx = (int)((double)k / num_failures_arg * num_failures);
                if (idx >= num_failures) idx = num_failures - 1;
                subsampled.push_back(failure_times_vec[idx]);
            }
            failure_times_vec = subsampled;
            num_failures = num_failures_arg;
        }

        /* Replace raw times with interval indices in data_vec (ALL rows) */
        for (int i = 0; i < n; i++) {
            double t = data_vec[(size_t)y_start * n + i];
            /* findInterval: count how many failure_times <= t */
            int idx = (int)(std::upper_bound(failure_times_vec.begin(),
                                             failure_times_vec.end(), t)
                            - failure_times_vec.begin());
            data_vec[(size_t)y_start * n + i] = (double)idx;
        }

        grf::ForestTrainer trainer = grf::survival_trainer(fast_logrank);
        grf::ForestPredictor predictor = grf::survival_predictor(
            resolved_threads, (size_t)num_failures, prediction_type);
        int out_col_start = nvar - n_output + 1;
        int n_written = 0;

        if (predict_mode) {
            int n_test = n - n_train;
            std::vector<double> train_vec, test_vec;
            split_train_test(data_vec, n, n_data_cols, n_train, n_x, train_vec, test_vec);
            grf::Data train_data(train_vec.data(), (size_t)n_train, (size_t)n_data_cols);
            set_data_indices(train_data, y_start, n_y, w_start, n_w, z_start, n_z);
            train_data.set_censor_index((size_t)censor_col);
            grf::Data test_data(test_vec.data(), (size_t)n_test, (size_t)n_x);

            snprintf(msg, sizeof(msg), "  Training survival forest (failures=%d)...\n", num_failures);
            SF_display(msg);
            grf::Forest forest = trainer.train(train_data, options);
            SF_display("  Forest trained. Predicting on new data...\n");
            predictions = predictor.predict(forest, train_data, test_data, false);

            for (int i = 0; i < n_test; i++) {
                const auto& pred = predictions[i].get_predictions();
                for (int k = 0; k < n_output && (size_t)k < pred.size(); k++) {
                    if (std::isfinite(pred[k])) {
                        SF_vstore(out_col_start + k, obs_map[n_train + i], pred[k]);
                    }
                }
                if (!pred.empty()) n_written++;
            }
        } else {
            /* Recreate the Data object with relabeled outcomes */
            grf::Data data_surv(data_vec.data(), n, n_data_cols);
            data_surv.set_outcome_index(y_start);
            data_surv.set_censor_index((size_t)censor_col);

            snprintf(msg, sizeof(msg), "  Training survival forest (failures=%d)...\n", num_failures);
            SF_display(msg);
            grf::Forest forest = trainer.train(data_surv, options);
            SF_display("  Forest trained.\n");

            SF_display("  Computing predictions...\n");
            predictions = predictor.predict_oob(forest, data_surv, false);

            for (int i = 0; i < n; i++) {
                const auto& pred = predictions[i].get_predictions();
                for (int k = 0; k < n_output && (size_t)k < pred.size(); k++) {
                    if (std::isfinite(pred[k])) {
                        SF_vstore(out_col_start + k, obs_map[i], pred[k]);
                    }
                }
                if (!pred.empty()) n_written++;
            }
        }

        snprintf(msg, sizeof(msg), "  Wrote %d survival predictions.\n", n_written);
        SF_display(msg);

    } else if (forest_type == "causal_survival") {
        /* ---- Causal Survival Forest ----
         * Data layout: X Y(time) W(treatment) censor numerator denominator
         * These are pre-computed by the .ado wrapper following grf's pipeline.
         */
        int stabilize = (argc > 20) ? parse_int(argv[20], 1) : 1;

        /* Set special indices for causal survival
         * The .ado sets up: X time W censor cs_numer cs_denom
         * outcome=time, treatment=W, censor, cs_numerator, cs_denominator */
        int cs_numer_col = (argc > 21) ? parse_int(argv[21], -1) : -1;
        int cs_denom_col = (argc > 22) ? parse_int(argv[22], -1) : -1;
        int censor_col_cs = (argc > 23) ? parse_int(argv[23], -1) : -1;

        /* target: 1=RMST (default), 2=survival probability
         * This affects nuisance column interpretation. The .ado computes
         * numerator/denominator differently for each target type. */
        int cs_target = (argc > 24) ? parse_int(argv[24], 1) : 1;
        const char* target_label = (cs_target == 2) ? "survival probability" : "RMST";
        snprintf(msg, sizeof(msg), "  Target estimand: %s (target=%d)\n", target_label, cs_target);
        SF_display(msg);

        if (cs_numer_col >= 0) data.set_causal_survival_numerator_index((size_t)cs_numer_col);
        if (cs_denom_col >= 0) data.set_causal_survival_denominator_index((size_t)cs_denom_col);
        if (censor_col_cs >= 0) data.set_censor_index((size_t)censor_col_cs);
        /* The CausalSurvivalSplittingRule uses get_instrument() for treatment.
         * Set instrument index to the treatment column (W). */
        data.set_instrument_index((size_t)w_start);
        grf::ForestTrainer trainer = grf::causal_survival_trainer((stabilize != 0));
        grf::ForestPredictor predictor = grf::causal_survival_predictor(resolved_threads);
        bool est_var = (estimate_variance != 0);
        int out_col_pred = nvar - n_output + 1;
        int out_col_var = (n_output >= 2 && est_var) ? out_col_pred + 1 : 0;
        int n_written = 0;

        if (predict_mode) {
            int n_test = n - n_train;
            std::vector<double> train_vec, test_vec;
            split_train_test(data_vec, n, n_data_cols, n_train, n_x, train_vec, test_vec);
            grf::Data train_data(train_vec.data(), (size_t)n_train, (size_t)n_data_cols);
            set_data_indices(train_data, y_start, n_y, w_start, n_w, z_start, n_z);
            if (cs_numer_col >= 0) train_data.set_causal_survival_numerator_index((size_t)cs_numer_col);
            if (cs_denom_col >= 0) train_data.set_causal_survival_denominator_index((size_t)cs_denom_col);
            if (censor_col_cs >= 0) train_data.set_censor_index((size_t)censor_col_cs);
            train_data.set_instrument_index((size_t)w_start);
            grf::Data test_data(test_vec.data(), (size_t)n_test, (size_t)n_x);

            SF_display("  Training causal survival forest...\n");
            grf::Forest forest = trainer.train(train_data, options);
            SF_display("  Forest trained. Predicting on new data...\n");
            predictions = predictor.predict(forest, train_data, test_data, est_var);

            for (int i = 0; i < n_test; i++) {
                const auto& pred = predictions[i].get_predictions();
                if (!pred.empty() && std::isfinite(pred[0])) {
                    SF_vstore(out_col_pred, obs_map[n_train + i], pred[0]);
                    n_written++;
                }
                if (out_col_var > 0 && predictions[i].contains_variance_estimates()) {
                    const auto& var_est = predictions[i].get_variance_estimates();
                    if (!var_est.empty()) {
                        SF_vstore(out_col_var, obs_map[n_train + i], var_est[0]);
                    }
                }
            }
        } else {
            SF_display("  Training causal survival forest...\n");
            grf::Forest forest = trainer.train(data, options);
            SF_display("  Forest trained.\n");

            SF_display("  Computing predictions...\n");
            predictions = predictor.predict_oob(forest, data, est_var);

            for (int i = 0; i < n; i++) {
                const auto& pred = predictions[i].get_predictions();
                if (!pred.empty() && std::isfinite(pred[0])) {
                    SF_vstore(out_col_pred, obs_map[i], pred[0]);
                    n_written++;
                }
                if (out_col_var > 0 && predictions[i].contains_variance_estimates()) {
                    const auto& var_est = predictions[i].get_variance_estimates();
                    if (!var_est.empty()) {
                        SF_vstore(out_col_var, obs_map[i], var_est[0]);
                    }
                }
            }
        }

        snprintf(msg, sizeof(msg), "  Wrote %d causal survival predictions.\n", n_written);
        SF_display(msg);

    } else if (forest_type == "multi_arm_causal") {
        /* ---- Multi-arm Causal Forest ----
         * multi_causal_trainer with num_treatments > 1
         */
        int stabilize = (argc > 20) ? parse_int(argv[20], 1) : 1;

        grf::ForestTrainer trainer = grf::multi_causal_trainer(
            (size_t)n_w, (size_t)n_y, (stabilize != 0));
        grf::ForestPredictor predictor = grf::multi_causal_predictor(
            resolved_threads, (size_t)n_w, (size_t)n_y);
        bool est_var = (estimate_variance != 0);
        int out_col_start = nvar - n_output + 1;
        int n_written = 0;
        int n_effects = n_w * n_y;

        if (predict_mode) {
            int n_test = n - n_train;
            std::vector<double> train_vec, test_vec;
            split_train_test(data_vec, n, n_data_cols, n_train, n_x, train_vec, test_vec);
            grf::Data train_data(train_vec.data(), (size_t)n_train, (size_t)n_data_cols);
            set_data_indices(train_data, y_start, n_y, w_start, n_w, z_start, n_z);
            grf::Data test_data(test_vec.data(), (size_t)n_test, (size_t)n_x);

            SF_display("  Training multi-arm causal forest...\n");
            grf::Forest forest = trainer.train(train_data, options);
            SF_display("  Forest trained. Predicting on new data...\n");
            predictions = predictor.predict(forest, train_data, test_data, est_var);

            for (int i = 0; i < n_test; i++) {
                const auto& pred = predictions[i].get_predictions();
                for (int k = 0; k < n_effects && k < n_output; k++) {
                    if ((size_t)k < pred.size() && std::isfinite(pred[k])) {
                        SF_vstore(out_col_start + k, obs_map[n_train + i], pred[k]);
                    }
                }
                if (est_var && predictions[i].contains_variance_estimates()) {
                    const auto& var_est = predictions[i].get_variance_estimates();
                    for (int k = 0; k < n_effects && (n_effects + k) < n_output; k++) {
                        if ((size_t)k < var_est.size()) {
                            SF_vstore(out_col_start + n_effects + k, obs_map[n_train + i], var_est[k]);
                        }
                    }
                }
                if (!pred.empty()) n_written++;
            }
        } else {
            SF_display("  Training multi-arm causal forest...\n");
            grf::Forest forest = trainer.train(data, options);
            SF_display("  Forest trained.\n");

            SF_display("  Computing predictions...\n");
            predictions = predictor.predict_oob(forest, data, est_var);

            for (int i = 0; i < n; i++) {
                const auto& pred = predictions[i].get_predictions();
                for (int k = 0; k < n_effects && k < n_output; k++) {
                    if ((size_t)k < pred.size() && std::isfinite(pred[k])) {
                        SF_vstore(out_col_start + k, obs_map[i], pred[k]);
                    }
                }
                /* Variance estimates after predictions */
                if (est_var && predictions[i].contains_variance_estimates()) {
                    const auto& var_est = predictions[i].get_variance_estimates();
                    for (int k = 0; k < n_effects && (n_effects + k) < n_output; k++) {
                        if ((size_t)k < var_est.size()) {
                            SF_vstore(out_col_start + n_effects + k, obs_map[i], var_est[k]);
                        }
                    }
                }
                if (!pred.empty()) n_written++;
            }
        }

        snprintf(msg, sizeof(msg), "  Wrote %d multi-arm predictions (%d effects).\n",
                 n_written, n_effects);
        SF_display(msg);

    } else if (forest_type == "multi_regression") {
        /* ---- Multi-outcome Regression Forest ---- */
        grf::ForestTrainer trainer = grf::multi_regression_trainer((size_t)n_y);
        grf::ForestPredictor predictor = grf::multi_regression_predictor(
            resolved_threads, (size_t)n_y);
        int out_col_start = nvar - n_output + 1;
        int n_written = 0;

        if (predict_mode) {
            int n_test = n - n_train;
            std::vector<double> train_vec, test_vec;
            split_train_test(data_vec, n, n_data_cols, n_train, n_x, train_vec, test_vec);
            grf::Data train_data(train_vec.data(), (size_t)n_train, (size_t)n_data_cols);
            set_data_indices(train_data, y_start, n_y, w_start, n_w, z_start, n_z);
            grf::Data test_data(test_vec.data(), (size_t)n_test, (size_t)n_x);

            snprintf(msg, sizeof(msg), "  Training multi-regression forest (%d outcomes)...\n", n_y);
            SF_display(msg);
            grf::Forest forest = trainer.train(train_data, options);
            SF_display("  Forest trained. Predicting on new data...\n");
            predictions = predictor.predict(forest, train_data, test_data, false);

            for (int i = 0; i < n_test; i++) {
                const auto& pred = predictions[i].get_predictions();
                for (int k = 0; k < n_y && k < n_output; k++) {
                    if ((size_t)k < pred.size() && std::isfinite(pred[k])) {
                        SF_vstore(out_col_start + k, obs_map[n_train + i], pred[k]);
                    }
                }
                if (!pred.empty()) n_written++;
            }
        } else {
            snprintf(msg, sizeof(msg), "  Training multi-regression forest (%d outcomes)...\n", n_y);
            SF_display(msg);
            grf::Forest forest = trainer.train(data, options);
            SF_display("  Forest trained.\n");

            SF_display("  Computing predictions...\n");
            predictions = predictor.predict_oob(forest, data, false);

            for (int i = 0; i < n; i++) {
                const auto& pred = predictions[i].get_predictions();
                for (int k = 0; k < n_y && k < n_output; k++) {
                    if ((size_t)k < pred.size() && std::isfinite(pred[k])) {
                        SF_vstore(out_col_start + k, obs_map[i], pred[k]);
                    }
                }
                if (!pred.empty()) n_written++;
            }
        }

        snprintf(msg, sizeof(msg), "  Wrote %d multi-regression predictions.\n", n_written);
        SF_display(msg);

    } else if (forest_type == "variable_importance") {
        /* ---- Variable Importance ----
         * Trains a regression forest and computes split-frequency-based
         * variable importance. Writes to Stata matrix.
         */
        int max_depth = (argc > 20) ? parse_int(argv[20], 4) : 4;

        SF_display("  Training forest for variable importance...\n");
        grf::ForestTrainer trainer = grf::regression_trainer();
        grf::Forest forest = trainer.train(data, options);

        grf::SplitFrequencyComputer sfc;
        std::vector<std::vector<size_t>> freqs = sfc.compute(forest, (size_t)max_depth);

        /* Compute variable importance as weighted split frequencies
         * (matches R's variable_importance which is a depth-weighted sum) */
        std::vector<double> vi(n_x, 0.0);
        double total = 0.0;
        for (size_t d = 0; d < freqs.size(); d++) {
            double depth_weight = 1.0 / (1.0 + (double)d);
            for (size_t v = 0; v < (size_t)n_x && v < freqs[d].size(); v++) {
                vi[v] += depth_weight * (double)freqs[d][v];
                total += depth_weight * (double)freqs[d][v];
            }
        }
        if (total > 0.0) {
            for (int v = 0; v < n_x; v++) {
                vi[v] /= total;
            }
        }

        /* Store as Stata scalars: _grf_vi_1, _grf_vi_2, ... */
        for (int v = 0; v < n_x; v++) {
            char scname[64];
            snprintf(scname, sizeof(scname), "_grf_vi_%d", v + 1);
            SF_scal_save(scname, vi[v]);
        }

        /* Also store count */
        SF_scal_save("_grf_vi_n", (double)n_x);

        snprintf(msg, sizeof(msg), "  Computed variable importance for %d variables.\n", n_x);
        SF_display(msg);

    } else if (forest_type == "ll_regression") {
        /* ---- Local Linear Regression Forest ---- */
        int enable_ll_split = (argc > 20) ? parse_int(argv[20], 0) : 0;
        double ll_lambda = (argc > 21) ? parse_double(argv[21], 0.1) : 0.1;
        int ll_weight_penalty_flag = (argc > 22) ? parse_int(argv[22], 0) : 0;
        int ll_split_cutoff = (argc > 23) ? parse_int(argv[23], 0) : 0;

        if (ll_split_cutoff <= 0) {
            ll_split_cutoff = (int)std::ceil(std::sqrt((double)n));
        }

        // All X variables used for local linear correction (0-indexed)
        std::vector<size_t> ll_split_vars;
        for (int j = 0; j < n_x; j++) {
            ll_split_vars.push_back((size_t)j);
        }

        // Compute overall_beta for regularization fallback (only if enable_ll_split && cutoff > 0)
        // overall_beta = (D'D + lambda*J)^{-1} D'Y
        // where D = [1, X], J = diag(0, 1, ..., 1)
        std::vector<double> overall_beta;
        if (enable_ll_split && ll_split_cutoff > 0) {
            int dim = n_x + 1; // intercept + p predictors
            // D'D matrix (dim x dim)
            std::vector<double> DtD(dim * dim, 0.0);
            std::vector<double> DtY(dim, 0.0);

            for (int i = 0; i < n; i++) {
                double yi = data_vec[(size_t)y_start * n + i];
                // D[i, 0] = 1 (intercept)
                DtD[0] += 1.0;
                DtY[0] += yi;
                for (int j = 0; j < n_x; j++) {
                    double xij = data_vec[(size_t)j * n + i];
                    DtD[(j+1) * dim + 0] += xij;       // D'D[j+1, 0]
                    DtD[0 * dim + (j+1)] += xij;       // D'D[0, j+1]
                    DtY[j+1] += xij * yi;
                    for (int k = 0; k <= j; k++) {
                        double xik = data_vec[(size_t)k * n + i];
                        DtD[(j+1) * dim + (k+1)] += xij * xik;
                        if (k != j) DtD[(k+1) * dim + (j+1)] += xij * xik;
                    }
                }
            }

            // Add ridge penalty: lambda to diagonal (skip intercept)
            for (int j = 1; j < dim; j++) {
                DtD[j * dim + j] += ll_lambda;
            }

            // Solve using Eigen LDLT decomposition
            Eigen::MatrixXd A(dim, dim);
            Eigen::VectorXd b(dim);
            for (int i = 0; i < dim; i++) {
                b(i) = DtY[i];
                for (int j = 0; j < dim; j++) {
                    A(i, j) = DtD[i * dim + j];
                }
            }
            Eigen::VectorXd beta = A.ldlt().solve(b);
            overall_beta.resize(dim);
            for (int i = 0; i < dim; i++) overall_beta[i] = beta(i);
        }

        // Create trainer
        grf::ForestTrainer trainer = enable_ll_split
            ? grf::ll_regression_trainer(
                ll_lambda, (ll_weight_penalty_flag != 0),
                overall_beta, (size_t)ll_split_cutoff, ll_split_vars)
            : grf::regression_trainer();

        // Create LL predictor (single lambda)
        std::vector<double> lambdas = {ll_lambda};
        grf::ForestPredictor predictor = grf::ll_regression_predictor(
            resolved_threads, lambdas, (ll_weight_penalty_flag != 0), ll_split_vars);

        bool est_var_ll = (estimate_variance != 0);
        int out_col_pred = nvar - n_output + 1;
        int out_col_var = (n_output >= 2 && est_var_ll) ? out_col_pred + 1 : 0;
        int n_written = 0;

        if (predict_mode) {
            int n_test = n - n_train;
            std::vector<double> train_vec, test_vec;
            split_train_test(data_vec, n, n_data_cols, n_train, n_x, train_vec, test_vec);
            grf::Data train_data(train_vec.data(), (size_t)n_train, (size_t)n_data_cols);
            set_data_indices(train_data, y_start, n_y, w_start, n_w, z_start, n_z);
            grf::Data test_data(test_vec.data(), (size_t)n_test, (size_t)n_x);

            SF_display("  Training LL regression forest...\n");
            grf::Forest forest = trainer.train(train_data, options);
            SF_display("  Forest trained. Predicting on new data...\n");
            predictions = predictor.predict(forest, train_data, test_data, est_var_ll);

            for (int i = 0; i < n_test; i++) {
                const auto& pred = predictions[i].get_predictions();
                if (!pred.empty() && std::isfinite(pred[0])) {
                    SF_vstore(out_col_pred, obs_map[n_train + i], pred[0]);
                    n_written++;
                }
                if (out_col_var > 0 && predictions[i].contains_variance_estimates()) {
                    const auto& var_est = predictions[i].get_variance_estimates();
                    if (!var_est.empty()) {
                        SF_vstore(out_col_var, obs_map[n_train + i], var_est[0]);
                    }
                }
            }
        } else {
            SF_display("  Training LL regression forest...\n");
            grf::Forest forest = trainer.train(data, options);
            SF_display("  Forest trained.\n");

            SF_display("  Computing LL predictions...\n");
            predictions = predictor.predict_oob(forest, data, est_var_ll);

            for (int i = 0; i < n; i++) {
                const auto& pred = predictions[i].get_predictions();
                if (!pred.empty() && std::isfinite(pred[0])) {
                    SF_vstore(out_col_pred, obs_map[i], pred[0]);
                    n_written++;
                }
                if (out_col_var > 0 && predictions[i].contains_variance_estimates()) {
                    const auto& var_est = predictions[i].get_variance_estimates();
                    if (!var_est.empty()) {
                        SF_vstore(out_col_var, obs_map[i], var_est[0]);
                    }
                }
            }
        }

        snprintf(msg, sizeof(msg), "  Wrote %d LL regression predictions.\n", n_written);
        SF_display(msg);

    } else if (forest_type == "boosted_regression") {
        /* ---- Boosted Regression Forest ----
         * Iteratively fit regression forests to residuals.
         * If boost_steps=0, auto-tune via cross-validation.
         */
        int boost_steps = (argc > 20) ? parse_int(argv[20], 0) : 0;
        double boost_error_reduction = (argc > 21) ? parse_double(argv[21], 0.97) : 0.97;
        int boost_max_steps = (argc > 22) ? parse_int(argv[22], 5) : 5;
        int boost_trees_tune = (argc > 23) ? parse_int(argv[23], 10) : 10;
        int stabilize = (argc > 24) ? parse_int(argv[24], 1) : 1;

        /* regression_trainer() does not use stabilize_splits (no treatment to stabilize).
         * The flag is accepted for interface consistency and reserved for future use. */
        (void)stabilize;
        grf::ForestTrainer trainer = grf::regression_trainer();
        grf::ForestPredictor predictor = grf::regression_predictor(resolved_threads);

        int out_col_pred = nvar - n_output + 1;

        // Accumulate predictions
        std::vector<double> y_hat(n, 0.0);
        double prev_mean_error = 1e30;
        int actual_steps = 0;

        // Working copy of data for residuals
        std::vector<double> resid_data = data_vec;

        for (int step = 0; step < (boost_steps > 0 ? boost_steps : boost_max_steps); step++) {

            // Auto-tune: check if another step improves enough
            if (boost_steps == 0 && step > 0) {
                // Fit a small forest on residuals to estimate error
                // Use ci_group_size=2 to get variance estimates for debiased error
                grf::ForestOptions tune_options(
                    (grf::uint)boost_trees_tune,
                    (size_t)2, // ci_group_size=2 for variance
                    sample_fraction, (grf::uint)mtry, (grf::uint)min_node_size,
                    (honesty != 0), honesty_fraction, (honesty_prune != 0),
                    alpha, imbalance_pen, (grf::uint)num_threads, (grf::uint)(seed + step),
                    legacy_seed, clusters, samples_per_cluster);

                grf::Data tune_data(resid_data.data(), (size_t)n, (size_t)n_data_cols);
                set_data_indices(tune_data, y_start, n_y, w_start, n_w, z_start, n_z);

                grf::Forest tune_forest = trainer.train(tune_data, tune_options);
                auto tune_preds = predictor.predict_oob(tune_forest, tune_data, true);

                // Compute mean debiased error: mean((Y_resid - Y_hat)^2 - var_hat)
                double sum_debiased = 0.0;
                int cnt = 0;
                for (int i = 0; i < n; i++) {
                    const auto& p = tune_preds[i].get_predictions();
                    if (p.empty() || !std::isfinite(p[0])) continue;
                    double resid_i = resid_data[(size_t)y_start * n + i];
                    double err = (resid_i - p[0]) * (resid_i - p[0]);
                    if (tune_preds[i].contains_variance_estimates()) {
                        const auto& v = tune_preds[i].get_variance_estimates();
                        if (!v.empty() && std::isfinite(v[0])) {
                            err -= v[0];
                        }
                    }
                    sum_debiased += err;
                    cnt++;
                }
                double mean_error = (cnt > 0) ? sum_debiased / cnt : 1e30;

                if (mean_error > boost_error_reduction * prev_mean_error) {
                    snprintf(msg, sizeof(msg),
                        "  Boosting stopped at step %d (error not improving).\n", step);
                    SF_display(msg);
                    break;
                }
                prev_mean_error = mean_error;
            }

            // Fit full forest on current residuals
            snprintf(msg, sizeof(msg), "  Boosting step %d/%d...\n",
                     step + 1, boost_steps > 0 ? boost_steps : boost_max_steps);
            SF_display(msg);

            // Use seed + step for different randomization each step
            // When auto-tuning (boost_steps==0), need ci_group_size>=2 for variance
            int step_ci = (boost_steps == 0 && ci_group_size < 2) ? 2 : ci_group_size;
            grf::ForestOptions step_options(
                (grf::uint)num_trees,
                (size_t)step_ci,
                sample_fraction, (grf::uint)mtry, (grf::uint)min_node_size,
                (honesty != 0), honesty_fraction, (honesty_prune != 0),
                alpha, imbalance_pen, (grf::uint)num_threads, (grf::uint)(seed + step),
                legacy_seed, clusters, samples_per_cluster);

            grf::Data step_data(resid_data.data(), (size_t)n, (size_t)n_data_cols);
            set_data_indices(step_data, y_start, n_y, w_start, n_w, z_start, n_z);

            grf::Forest forest = trainer.train(step_data, step_options);
            auto step_preds = predictor.predict_oob(forest, step_data, (boost_steps == 0));

            // Accumulate predictions and compute residuals for next step
            for (int i = 0; i < n; i++) {
                const auto& p = step_preds[i].get_predictions();
                if (!p.empty() && std::isfinite(p[0])) {
                    y_hat[i] += p[0];
                }
            }

            // Update residual data: Y_resid = Y_orig - Y_hat
            for (int i = 0; i < n; i++) {
                double y_orig = data_vec[(size_t)y_start * n + i];
                resid_data[(size_t)y_start * n + i] = y_orig - y_hat[i];
            }

            // Compute debiased error for this step (for auto-tune comparison)
            if (boost_steps == 0) {
                double sum_debiased = 0.0;
                int cnt = 0;
                for (int i = 0; i < n; i++) {
                    const auto& p = step_preds[i].get_predictions();
                    if (p.empty() || !std::isfinite(p[0])) continue;
                    double resid_i = resid_data[(size_t)y_start * n + i]; // this is already updated to new residual
                    // Use the step's residual before update for error: orig Y_resid - pred
                    double actual_resid_before = data_vec[(size_t)y_start * n + i] - (y_hat[i] - p[0]);
                    double err = (actual_resid_before - p[0]) * (actual_resid_before - p[0]);
                    if (step_preds[i].contains_variance_estimates()) {
                        const auto& v = step_preds[i].get_variance_estimates();
                        if (!v.empty() && std::isfinite(v[0])) {
                            err -= v[0];
                        }
                    }
                    sum_debiased += err;
                    cnt++;
                }
                if (cnt > 0) prev_mean_error = sum_debiased / cnt;
            }

            actual_steps++;
        }

        // Write final accumulated predictions
        int n_written = 0;
        for (int i = 0; i < n; i++) {
            if (std::isfinite(y_hat[i])) {
                SF_vstore(out_col_pred, obs_map[i], y_hat[i]);
                n_written++;
            }
        }

        // Store number of boosting steps as scalar
        SF_scal_save("_grf_boost_steps", (double)actual_steps);

        snprintf(msg, sizeof(msg), "  Wrote %d boosted predictions (%d steps).\n",
                 n_written, actual_steps);
        SF_display(msg);

    } else if (forest_type == "lm_forest") {
        /* ---- LM Forest ----
         * Conditional linear model: Y = c(x) + h_1(x)*W_1 + ... + h_K(x)*W_K
         * The .ado handles centering (Y - Y.hat, W - W.hat).
         * C++ training: same as multi_arm_causal (multi_causal_trainer).
         * The treatment columns are the "regressors" W.
         */
        int stabilize = (argc > 20) ? parse_int(argv[20], 0) : 0;

        grf::ForestTrainer trainer = grf::multi_causal_trainer(
            (size_t)n_w, (size_t)n_y, (stabilize != 0));
        grf::ForestPredictor predictor = grf::multi_causal_predictor(
            resolved_threads, (size_t)n_w, (size_t)n_y);
        bool est_var = (estimate_variance != 0);
        int out_col_start = nvar - n_output + 1;
        int n_written = 0;

        // Output: n_w * n_y coefficients per observation
        int n_coefs = n_w * n_y;

        if (predict_mode) {
            int n_test = n - n_train;
            std::vector<double> train_vec, test_vec;
            split_train_test(data_vec, n, n_data_cols, n_train, n_x, train_vec, test_vec);
            grf::Data train_data(train_vec.data(), (size_t)n_train, (size_t)n_data_cols);
            set_data_indices(train_data, y_start, n_y, w_start, n_w, z_start, n_z);
            grf::Data test_data(test_vec.data(), (size_t)n_test, (size_t)n_x);

            SF_display("  Training LM forest...\n");
            grf::Forest forest = trainer.train(train_data, options);
            SF_display("  Forest trained. Predicting on new data...\n");
            predictions = predictor.predict(forest, train_data, test_data, est_var);

            for (int i = 0; i < n_test; i++) {
                const auto& pred = predictions[i].get_predictions();
                for (int k = 0; k < n_coefs && k < (int)pred.size(); k++) {
                    if (std::isfinite(pred[k])) {
                        SF_vstore(out_col_start + k, obs_map[n_train + i], pred[k]);
                    }
                }
                if (!pred.empty()) n_written++;
                // Variance estimates
                if (est_var && predictions[i].contains_variance_estimates()) {
                    const auto& var_est = predictions[i].get_variance_estimates();
                    for (int k = 0; k < (int)var_est.size(); k++) {
                        int var_col = out_col_start + n_coefs + k;
                        if (var_col <= nvar && std::isfinite(var_est[k])) {
                            SF_vstore(var_col, obs_map[n_train + i], var_est[k]);
                        }
                    }
                }
            }
        } else {
            SF_display("  Training LM forest...\n");
            grf::Forest forest = trainer.train(data, options);
            SF_display("  Forest trained.\n");

            SF_display("  Computing LM predictions...\n");
            predictions = predictor.predict_oob(forest, data, est_var);

            for (int i = 0; i < n; i++) {
                const auto& pred = predictions[i].get_predictions();
                for (int k = 0; k < n_coefs && k < (int)pred.size(); k++) {
                    if (std::isfinite(pred[k])) {
                        SF_vstore(out_col_start + k, obs_map[i], pred[k]);
                    }
                }
                if (!pred.empty()) n_written++;
                if (est_var && predictions[i].contains_variance_estimates()) {
                    const auto& var_est = predictions[i].get_variance_estimates();
                    for (int k = 0; k < (int)var_est.size(); k++) {
                        int var_col = out_col_start + n_coefs + k;
                        if (var_col <= nvar && std::isfinite(var_est[k])) {
                            SF_vstore(var_col, obs_map[i], var_est[k]);
                        }
                    }
                }
            }
        }

        snprintf(msg, sizeof(msg), "  Wrote %d LM forest predictions (%d coefficients each).\n",
                 n_written, n_coefs);
        SF_display(msg);

    } else {
        snprintf(msg, sizeof(msg), "GRF error: unknown forest type '%s'\n",
                 forest_type.c_str());
        SF_error(msg);
        return 198;
    }

    } catch (const std::exception& e) {
        snprintf(msg, sizeof(msg), "GRF C++ exception: %s\n", e.what());
        SF_error(msg);
        return 198;
    } catch (...) {
        SF_error("GRF unknown C++ exception.\n");
        return 198;
    }

    return 0;
}
