# Introspection Discrepancies vs R `grf` (Current Baseline `v2.5.0`)

## What is implemented (partial introspection)

- `grf_forest_summary`: model-level metadata summary; with `, all` lists all stored `e()` scalars/macros.
- `grf_tree_summary` / `grf_get_tree`: tree-index metadata checks against fitted model state.
- `grf_get_leaf_node`: prediction-based leaf-group proxy (`xtile` grouping).
- `grf_get_forest_weights`: proxy forest weights using prediction distance and optional feature-space distance (`xvars()`, `predweight()`).
- `grf_split_frequencies`: depth-aggregated split-frequency proxy via variable-importance backend.
- `grf_plot_tree`: plot wrapper over split-frequency proxy.

## Remaining discrepancies (not substantively equivalent to full R object introspection)

1. Exact tree node structure is not exposed.
   - R-style raw node-by-node extraction with split thresholds, child pointers, and per-node sample membership is not available.

2. Exact forest kernel weights are not exposed.
   - R-style exact weights from internal forest traversal are not returned; Stata command returns a deterministic proxy.

3. Exact leaf assignment APIs are not exposed.
   - Stata leaf IDs are proxy groups from predictions, not exact trained-tree terminal node IDs.

4. Plot methods are utility wrappers, not full S3 method parity.
   - Plots are built from proxy diagnostics and generated variables, not full in-memory forest objects.

5. Full object-level mutation/composition workflows are constrained.
   - R-style workflows that depend on complete in-memory forest object state and internals are only partially representable through plugin boundaries.

## Practical impact

- Standard applied workflows (fit, predict, ATE/BLP/RATE/scores, diagnostics) remain supported.
- Method-development and deep debugging workflows requiring exact internals may need R for full fidelity.
