{smcl}
{* *! version 0.1.0}{...}
{viewerjumpto "Syntax" "grf_probability_forest##syntax"}{...}
{viewerjumpto "Description" "grf_probability_forest##description"}{...}
{viewerjumpto "Options" "grf_probability_forest##options"}{...}
{viewerjumpto "Examples" "grf_probability_forest##examples"}{...}
{viewerjumpto "Stored results" "grf_probability_forest##results"}{...}
{viewerjumpto "References" "grf_probability_forest##references"}{...}

{title:Title}

{phang}
{bf:grf_probability_forest} {hline 2} Probability forest for multi-class classification

{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmd:grf_probability_forest}
{it:depvar}
{it:indepvars}
{ifin}{cmd:,}
{opt gen:erate(stub)}
[{it:options}]

{synoptset 28 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Required}
{synopt:{opt gen:erate(stub)}}stub for output variables; creates {it:stub}{cmd:_c0}, {it:stub}{cmd:_c1}, ...{p_end}

{syntab:Classification}
{synopt:{opt ncl:asses(#)}}number of classes; default {cmd:0} (auto-detect from data){p_end}

{syntab:Forest tuning}
{synopt:{opt nt:rees(#)}}number of trees; default is {cmd:2000}{p_end}
{synopt:{opt seed(#)}}random seed; default is {cmd:42}{p_end}
{synopt:{opt mtry(#)}}variables to consider at each split; default {cmd:0} (auto){p_end}
{synopt:{opt minn:odesize(#)}}minimum leaf size; default is {cmd:5}{p_end}
{synopt:{opt sample:frac(#)}}fraction of observations per tree; default is {cmd:0.5}{p_end}
{synopt:{opt numt:hreads(#)}}number of threads; default {cmd:0} (all cores){p_end}

{syntab:Honesty}
{synopt:{opt ho:nesty}}use honest splitting (the default){p_end}
{synopt:{opt noho:nesty}}disable honest splitting{p_end}
{synopt:{opt honesty:frac(#)}}fraction reserved for honesty; default is {cmd:0.5}{p_end}
{synopt:{opt honestyp:rune}}prune empty honest leaves (the default){p_end}
{synopt:{opt nohonestyp:rune}}do not prune empty honest leaves{p_end}

{syntab:Regularization}
{synopt:{opt al:pha(#)}}minimum fraction of data in each child; default is {cmd:0.05}{p_end}
{synopt:{opt imb:alancepenalty(#)}}penalty on imbalanced splits; default is {cmd:0.0}{p_end}

{syntab:Output}
{synopt:{opt re:place}}overwrite existing output variables{p_end}
{synoptline}

{marker description}{...}
{title:Description}

{pstd}
{cmd:grf_probability_forest} estimates class membership probabilities P(Y=k|X)
for a multi-class outcome using the probability forest method from the
generalized random forests framework. The dependent variable {it:depvar} must
contain non-negative integers (0, 1, 2, ..., K-1) representing class labels.

{pstd}
The command generates K output variables named {it:stub}{cmd:_c0},
{it:stub}{cmd:_c1}, ..., {it:stub}{cmd:_c}{it:K-1}, where each variable
contains the estimated probability that an observation belongs to that class.
The probabilities for each observation sum to 1.

{marker options}{...}
{title:Options}

{dlgtab:Required}

{phang}
{opt generate(stub)} specifies the stub name for the output variables.
For K classes, the command creates {it:stub}{cmd:_c0} through
{it:stub}{cmd:_c}{it:K-1}, each containing the estimated probability of the
corresponding class.

{dlgtab:Classification}

{phang}
{opt nclasses(#)} sets the number of classes. The default {cmd:0} auto-detects
from the data as max({it:depvar}) + 1. If specified, must exceed
max({it:depvar}).

{dlgtab:Forest tuning}

{phang}
{opt ntrees(#)} sets the number of trees grown in the forest.

{phang}
{opt seed(#)} sets the random number seed for reproducibility.

{phang}
{opt mtry(#)} sets the number of variables randomly sampled as candidates at
each split. The default {cmd:0} lets the library choose automatically.

{phang}
{opt minnodesize(#)} sets the minimum number of observations in each leaf node.

{phang}
{opt samplefrac(#)} sets the fraction of observations sampled without
replacement to build each tree.

{phang}
{opt numthreads(#)} sets the number of threads. The default {cmd:0} uses all
available cores.

{dlgtab:Honesty}

{phang}
{opt honesty} and {opt nohonesty} control whether trees use separate
subsamples for splitting and estimation. Honesty is enabled by default.

{phang}
{opt honestyfrac(#)} sets the fraction of the subsample reserved for honest
estimation.

{phang}
{opt honestyprune} and {opt nohonestyprune} control pruning of honest leaves
that contain no estimation-sample observations. Pruning is enabled by default.

{marker examples}{...}
{title:Examples}

{pstd}Setup with iris-like data:{p_end}
{phang2}{cmd:. clear}{p_end}
{phang2}{cmd:. set obs 500}{p_end}
{phang2}{cmd:. set seed 12345}{p_end}
{phang2}{cmd:. gen x1 = rnormal()}{p_end}
{phang2}{cmd:. gen x2 = rnormal()}{p_end}
{phang2}{cmd:. gen x3 = rnormal()}{p_end}
{phang2}{cmd:. gen byte class = cond(x1 + x2 > 1, 2, cond(x1 > 0, 1, 0))}{p_end}

{pstd}Basic probability forest (3 classes):{p_end}
{phang2}{cmd:. grf_probability_forest class x1 x2 x3, gen(prob)}{p_end}
{phang2}{cmd:. list prob_c0 prob_c1 prob_c2 in 1/5}{p_end}

{pstd}Predicted class with highest probability:{p_end}
{phang2}{cmd:. gen predicted = 0}{p_end}
{phang2}{cmd:. replace predicted = 1 if prob_c1 > prob_c0 & prob_c1 > prob_c2}{p_end}
{phang2}{cmd:. replace predicted = 2 if prob_c2 > prob_c0 & prob_c2 > prob_c1}{p_end}
{phang2}{cmd:. tab class predicted}{p_end}

{pstd}Binary classification (2 classes):{p_end}
{phang2}{cmd:. gen byte treat = (x1 > 0)}{p_end}
{phang2}{cmd:. grf_probability_forest treat x1 x2 x3, gen(phat)}{p_end}
{phang2}{cmd:. sum phat_c0 phat_c1}{p_end}

{marker results}{...}
{title:Stored results}

{pstd}{cmd:grf_probability_forest} stores the following in {cmd:e()}:

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(n_trees)}}number of trees{p_end}
{synopt:{cmd:e(seed)}}random seed{p_end}
{synopt:{cmd:e(mtry)}}number of split candidates{p_end}
{synopt:{cmd:e(min_node)}}minimum node size{p_end}
{synopt:{cmd:e(alpha)}}alpha parameter{p_end}
{synopt:{cmd:e(honesty)}}1 if honesty enabled, 0 otherwise{p_end}
{synopt:{cmd:e(n_classes)}}number of classes{p_end}

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:grf_probability_forest}{p_end}
{synopt:{cmd:e(forest_type)}}{cmd:probability}{p_end}
{synopt:{cmd:e(depvar)}}outcome variable name{p_end}
{synopt:{cmd:e(indepvars)}}predictor variable names{p_end}
{synopt:{cmd:e(predict_vars)}}names of all output variables{p_end}

{marker references}{...}
{title:References}

{pstd}
Athey, S., J. Tibshirani, and S. Wager. 2019.
Generalized Random Forests. {it:Annals of Statistics} 47(2): 1148-1178.

{title:Author}

{pstd}
GRF Stata plugin. Wraps the grf C++ library (grf-labs/grf, v2.5.0, GPL-3.0).
