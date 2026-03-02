{smcl}
{* *! version 0.1.0}{...}
{title:grf_get_leaf_node}

{pstd}
{cmd:grf_get_leaf_node} generates a leaf-node proxy assignment variable from
the current model's prediction column (quantile grouping).

{pstd}
Syntax:
{cmd:grf_get_leaf_node,} {opt generate(newvar)}
[{opt groups(#)} {opt replace}]

{pstd}
Stores {cmd:r(N)}, {cmd:r(groups)}, and {cmd:r(generate)}.
