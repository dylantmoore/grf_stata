{smcl}
{* *! version 0.2.0}{...}
{title:grf_get_forest_weights}

{pstd}
{cmd:grf_get_forest_weights} generates a proxy kernel-weight vector anchored at
a chosen observation, based on the current prediction variable and optionally
feature-space distances.

{pstd}
Syntax:
{cmd:grf_get_forest_weights,} {opt obs(#)} {opt generate(newvar)}
[{opt scale(#)} {opt xscale(#)} {opt predweight(#)} {opt xvars(varlist)} {opt replace}]

{pstd}
The generated weights sum to one over complete observations.

{pstd}
Behavior:

{phang}
If {cmd:xvars()} is omitted, weights use prediction distance only.

{phang}
If {cmd:xvars()} is provided, weights combine prediction and standardized
feature distances as:
{cmd:exp(-[predweight * pred_distance + (1-predweight) * x_distance])}.

{pstd}
{cmd:predweight()} must be in [0,1]. If {cmd:predweight() < 1}, {cmd:xvars()}
is required.
