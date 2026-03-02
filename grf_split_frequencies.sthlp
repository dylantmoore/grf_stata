{smcl}
{* *! version 0.1.0}{...}
{title:grf_split_frequencies}

{pstd}
{cmd:grf_split_frequencies} returns a depth-aggregated split-frequency proxy
matrix by delegating to the variable-importance analysis backend.

{pstd}
Syntax:
{cmd:grf_split_frequencies} {it:depvar indepvars} [{ifin}]
[{cmd:,} {opt ntrees(#)} {opt seed(#)} {opt maxdepth(#)} {opt decayexponent(#)}]

{pstd}
Stores matrix {cmd:r(split_frequencies)}.
