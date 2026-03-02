{smcl}
{* *! version 0.1.0}{...}
{title:grf_get_forest_weights}

{pstd}
{cmd:grf_get_forest_weights} generates a proxy kernel-weight vector anchored at
a chosen observation, based on the current prediction variable.

{pstd}
Syntax:
{cmd:grf_get_forest_weights,} {opt obs(#)} {opt generate(newvar)}
[{opt scale(#)} {opt replace}]

{pstd}
The generated weights sum to one over non-missing predictions.
