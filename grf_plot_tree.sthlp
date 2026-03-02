{smcl}
{* *! version 0.1.0}{...}
{title:grf_plot_tree}

{pstd}
{cmd:grf_plot_tree} draws a bar chart of split-frequency proxy importance from
the variable-importance backend.

{pstd}
Syntax:
{cmd:grf_plot_tree} {it:depvar indepvars} [{ifin}]
[{cmd:,} {opt ntrees(#)} {opt seed(#)} {opt maxdepth(#)} {opt name(name)} {opt replace}]

{pstd}
Stores matrix {cmd:r(importance)} and graph name in {cmd:r(graph_name)}.
