{smcl}
{* *! version 0.1.0}{...}
{viewerjumpto "Syntax" "grf_multi_regression_forest##syntax"}{...}
{viewerjumpto "Description" "grf_multi_regression_forest##description"}{...}
{viewerjumpto "Options" "grf_multi_regression_forest##options"}{...}
{viewerjumpto "Examples" "grf_multi_regression_forest##examples"}{...}
{viewerjumpto "Stored results" "grf_multi_regression_forest##results"}{...}
{viewerjumpto "References" "grf_multi_regression_forest##references"}{...}

{title:Title}

{phang}
{bf:grf_multi_regression_forest} {hline 2} Multi-output regression forest

{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmd:grf_multi_regression_forest}
{it:depvar1 depvar2 ...} {indepvars}
{ifin}{cmd:,}
{opt gen:erate(stub)}
{opt ndep(#)}
[{it:options}]

{synoptset 28 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Main}
{synopt:{opt gen:erate(stub)}}stub for prediction variables{p_end}
{synopt:{opt ndep(#)}}number of outcome (dependent) variables{p_end}
{synopt:{opt replace}}overwrite existing variables{p_end}

{syntab:Forest}
{synopt:{opt ntrees(#)}}number of trees; default {bf:2000}{p_end}
{synopt:{opt seed(#)}}random seed; default {bf:42}{p_end}
{synopt:{opt mtry(#)}}variables per split; default {bf:ceil(sqrt(p))}{p_end}
{synopt:{opt minnodesize(#)}}minimum node size; default {bf:5}{p_end}
{synopt:{opt samplefrac(#)}}sample fraction; default {bf:0.5}{p_end}
{synopt:{opt numthreads(#)}}threads; default {bf:0} (auto){p_end}

{syntab:Honesty}
{synopt:{opt honesty}/{opt nohonesty}}use honesty splitting; default {bf:honesty}{p_end}
{synopt:{opt honestyfrac(#)}}honesty fraction; default {bf:0.5}{p_end}

{syntab:Regularization}
{synopt:{opt alpha(#)}}balance parameter; default {bf:0.05}{p_end}
{synopt:{opt imbalancepenalty(#)}}imbalance penalty; default {bf:0.0}{p_end}
{synoptline}

{marker description}{...}
{title:Description}

{pstd}
{cmd:grf_multi_regression_forest} jointly predicts multiple outcome variables
using a single forest that accounts for cross-outcome structure. It
implements {cmd:multi_regression_forest()} from R's {bf:grf} package.

{pstd}
The first {opt ndep()} variables in the varlist are the outcome variables;
the remaining are predictors. Output variables are named {it:stub}{bf:_y}{it:k}
for k = 1, ..., K, one per outcome.

{pstd}
Use {cmd:grf_regression_forest} for single-outcome regression. Multi-output
regression is useful when outcomes are correlated and joint estimation
improves predictions.

{marker options}{...}
{title:Options}

{dlgtab:Main}

{phang}
{opt generate(stub)} specifies the stub for output prediction variables.

{phang}
{opt ndep(#)} specifies how many of the leading variables in the varlist
are outcomes (must be at least 2).

{marker examples}{...}
{title:Examples}

{pstd}Joint prediction of two outcomes:{p_end}
{phang}{cmd:. grf_multi_regression_forest y1 y2 x1 x2 x3, gen(pred) ndep(2) ntrees(2000)}{p_end}

{marker results}{...}
{title:Stored results}

{pstd}
{cmd:grf_multi_regression_forest} stores the following in {cmd:e()}:

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(n_trees)}}number of trees{p_end}
{synopt:{cmd:e(n_outcomes)}}number of outcome variables{p_end}
{synopt:{cmd:e(seed)}}random seed{p_end}

{p2col 5 24 28 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}grf_multi_regression_forest{p_end}
{synopt:{cmd:e(forest_type)}}multi_regression{p_end}

{marker references}{...}
{title:References}

{pstd}
Athey, S., J. Tibshirani, and S. Wager. 2019.
Generalized Random Forests. {it:Annals of Statistics} 47(2): 1148{c -}1178.

{title:Author}

{pstd}
GRF Stata plugin. Wraps the grf C++ library (grf-labs/grf, v2.5.0, GPL-3.0).
