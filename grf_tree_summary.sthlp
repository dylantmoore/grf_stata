{smcl}
{* *! version 0.1.0}{...}
{title:grf_tree_summary}

{pstd}
{cmd:grf_tree_summary} reports tree-level metadata for the most recent fitted
GRF model.

{pstd}
Syntax:
{cmd:grf_tree_summary} [{cmd:,} {opt tree(#)}]

{pstd}
Stores {cmd:r(tree)}, {cmd:r(n_trees)}, and {cmd:r(forest_type)}.

{pstd}
This is a metadata summary. Exact node-level tree internals are not currently
persisted through the plugin boundary.
