{smcl}
{* *! version 0.1.0}{...}
{viewerjumpto "Syntax" "grf_multi_arm_causal_forest##syntax"}{...}
{viewerjumpto "Description" "grf_multi_arm_causal_forest##description"}{...}
{viewerjumpto "Options" "grf_multi_arm_causal_forest##options"}{...}
{viewerjumpto "Examples" "grf_multi_arm_causal_forest##examples"}{...}
{viewerjumpto "Stored results" "grf_multi_arm_causal_forest##results"}{...}
{viewerjumpto "References" "grf_multi_arm_causal_forest##references"}{...}

{title:Title}

{phang}
{bf:grf_multi_arm_causal_forest} {hline 2} Multi-arm causal forest for
heterogeneous treatment effects with multiple discrete treatment arms

{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmd:grf_multi_arm_causal_forest}
{it:depvar} {it:treatvar} {indepvars}
{ifin}{cmd:,}
{opt gen:erate(stub)}
{opt ntreat(#)}
[{it:options}]

{synoptset 28 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Main}
{synopt:{opt gen:erate(stub)}}stub for CATE prediction variables{p_end}
{synopt:{opt ntreat(#)}}number of treatment indicator variables{p_end}
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

{syntab:Multi-arm}
{synopt:{opt stab:ilizesplits}/{opt nostab:ilizesplits}}stabilize splits; default {bf:stabilize}{p_end}
{synoptline}

{marker description}{...}
{title:Description}

{pstd}
{cmd:grf_multi_arm_causal_forest} estimates heterogeneous treatment effects
with multiple discrete treatment arms. It implements {cmd:multi_arm_causal_forest()}
from R's {bf:grf} package.

{pstd}
The {it:treatvar} should be a single categorical treatment variable taking
values 0, 1, ..., K-1, where 0 is the control arm. The {opt ntreat(1)}
option specifies 1 treatment variable.

{pstd}
Output variables are named {it:stub}{bf:_t}{it:k} for k = 1, ..., K-1,
each containing CATE estimates for treatment arm k vs. control.

{pstd}
The command fits nuisance models Y ~ X and W ~ X automatically, then
trains the multi-arm causal forest on centered data.

{marker options}{...}
{title:Options}

{dlgtab:Main}

{phang}
{opt generate(stub)} specifies the stub for output variables.

{phang}
{opt ntreat(#)} specifies the number of treatment indicator variables
in the varlist (placed between {it:depvar} and {it:indepvars}).

{marker examples}{...}
{title:Examples}

{pstd}Multi-arm treatment with 3 arms (0=control, 1, 2):{p_end}
{phang}{cmd:. grf_multi_arm_causal_forest y treat x1 x2 x3, gen(cate) ntreat(1) ntrees(2000)}{p_end}

{marker results}{...}
{title:Stored results}

{pstd}
{cmd:grf_multi_arm_causal_forest} stores the following in {cmd:e()}:

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(n_trees)}}number of trees{p_end}
{synopt:{cmd:e(seed)}}random seed{p_end}

{p2col 5 24 28 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}grf_multi_arm_causal_forest{p_end}
{synopt:{cmd:e(forest_type)}}multi_arm_causal{p_end}

{marker references}{...}
{title:References}

{pstd}
Athey, S., J. Tibshirani, and S. Wager. 2019.
Generalized Random Forests. {it:Annals of Statistics} 47(2): 1148{c -}1178.

{title:Author}

{pstd}
GRF Stata plugin. Wraps the grf C++ library (grf-labs/grf, v2.5.0, GPL-3.0).
