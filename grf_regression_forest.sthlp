{smcl}
{* *! version 0.1.0}{...}
{viewerjumpto "Syntax" "grf_regression_forest##syntax"}{...}
{viewerjumpto "Description" "grf_regression_forest##description"}{...}
{viewerjumpto "Options" "grf_regression_forest##options"}{...}
{viewerjumpto "Examples" "grf_regression_forest##examples"}{...}
{viewerjumpto "Stored results" "grf_regression_forest##results"}{...}
{viewerjumpto "References" "grf_regression_forest##references"}{...}

{title:Title}

{phang}
{bf:grf_regression_forest} {hline 2} Generalized Random Forest for regression

{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmd:grf_regression_forest}
{depvar}
{indepvars}
{ifin}{cmd:,}
{opt gen:erate(newvar)}
[{it:options}]

{synoptset 30 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Required}
{synopt:{opt gen:erate(newvar)}}name of variable to store out-of-bag predictions{p_end}

{syntab:Forest}
{synopt:{opt ntr:ees(#)}}number of trees; default is {cmd:ntrees(2000)}{p_end}
{synopt:{opt seed(#)}}random-number seed; default is {cmd:seed(42)}{p_end}
{synopt:{opt mtry(#)}}variables tried at each split; default is {cmd:mtry(0)} (= sqrt(p)){p_end}
{synopt:{opt minn:odesize(#)}}minimum leaf size; default is {cmd:minnodesize(5)}{p_end}
{synopt:{opt sample:frac(#)}}fraction of observations per tree; default is {cmd:samplefrac(0.5)}{p_end}
{synopt:{opt numt:hreads(#)}}threads for fitting; default is {cmd:numthreads(0)} (= all cores){p_end}

{syntab:Honesty}
{synopt:{opt hon:esty}}use honest splitting (the default){p_end}
{synopt:{opt nohon:esty}}disable honest splitting{p_end}
{synopt:{opt hon:estyfrac(#)}}fraction held out for honest re-estimation; default is {cmd:honestyfrac(0.5)}{p_end}
{synopt:{opt hon:estyprune}}prune empty honest leaves (the default){p_end}
{synopt:{opt nohon:estyprune}}keep empty honest leaves{p_end}

{syntab:Tuning}
{synopt:{opt alp:ha(#)}}maximum imbalance of a split; default is {cmd:alpha(0.05)}{p_end}
{synopt:{opt imb:alancepenalty(#)}}penalty on imbalanced splits; default is {cmd:imbalancepenalty(0.0)}{p_end}

{syntab:Variance estimation}
{synopt:{opt est:imatevariance}}compute variance of predictions{p_end}
{synopt:{opt varg:enerate(newvar)}}name for variance variable; default is {it:generate}{cmd:_var}{p_end}
{synopt:{opt cig:roupsize(#)}}cluster size for variance; default is {cmd:cigroupsize(1)}; forced to 2+ with {cmd:estimatevariance}{p_end}

{syntab:Output}
{synopt:{opt replace}}overwrite existing output variables{p_end}
{synoptline}

{marker description}{...}
{title:Description}

{pstd}
{cmd:grf_regression_forest} fits a generalized random forest for
non-parametric regression of {depvar} on {indepvars} using the grf C++
library. It produces out-of-bag predictions for every observation in the
estimation sample and optionally estimates prediction variance.

{pstd}
The method is adaptive: each prediction is a weighted average of {depvar}
where the weights are determined by how often observations fall into the
same terminal leaves across trees. Honesty (the default) uses separate
subsamples for splitting and estimation, yielding valid asymptotic
confidence intervals.

{marker options}{...}
{title:Options}

{dlgtab:Required}

{phang}
{opt generate(newvar)} creates {it:newvar} containing out-of-bag
predictions of E[Y|X=x] for each observation in the estimation sample.
Observations outside the {cmd:if}/{cmd:in} sample are set to missing.

{dlgtab:Forest}

{phang}
{opt ntrees(#)} sets the number of trees. More trees improve stability at
the cost of computation. Default is 2000.

{phang}
{opt seed(#)} sets the random-number seed for reproducibility. Default is 42.

{phang}
{opt mtry(#)} sets the number of candidate variables considered at each
split. The default {cmd:mtry(0)} uses floor(sqrt(p)) + 20 where p is the
number of predictors.

{phang}
{opt minnodesize(#)} sets the minimum number of observations in each tree
leaf. Default is 5.

{phang}
{opt samplefrac(#)} sets the fraction of observations drawn (without
replacement) for each tree. Default is 0.5.

{phang}
{opt numthreads(#)} sets the number of threads. Default is 0 (use all
available cores).

{dlgtab:Honesty}

{phang}
{opt honesty} | {opt nohonesty} controls whether the forest uses honest
estimation. With honesty (the default), each tree's subsample is split into
a partition used for determining splits and a partition used for populating
leaf estimates. This reduces bias and enables valid inference.

{phang}
{opt honestyfrac(#)} sets the fraction of each tree's subsample reserved
for the honest estimation partition. Default is 0.5.

{phang}
{opt honestyprune} | {opt nohonestyprune} controls whether leaves that
receive no honest-estimation observations are pruned. Pruning (the default)
prevents predictions from degenerating.

{dlgtab:Tuning}

{phang}
{opt alpha(#)} controls the maximum imbalance of a split. Each child node
must contain at least a fraction {it:alpha} of the parent's observations.
Default is 0.05.

{phang}
{opt imbalancepenalty(#)} penalizes splits that produce children of unequal
size. Default is 0.0 (no penalty).

{dlgtab:Variance estimation}

{phang}
{opt estimatevariance} requests that the forest compute a variance estimate
for each prediction. When specified, {cmd:cigroupsize} is forced to at
least 2.

{phang}
{opt vargenerate(newvar)} names the variable for variance estimates.  If
omitted, the default name is {it:generate}{cmd:_var}.

{phang}
{opt cigroupsize(#)} sets the number of trees in each half-sample group
used for variance estimation. Default is 1; forced to at least 2 when
{cmd:estimatevariance} is specified.

{dlgtab:Output}

{phang}
{opt replace} allows overwriting an existing variable with the same name as
{it:generate} (and {it:vargenerate}, if applicable).

{marker examples}{...}
{title:Examples}

{pstd}Setup{p_end}
{phang2}{cmd:. sysuse auto, clear}{p_end}

{pstd}Basic regression forest{p_end}
{phang2}{cmd:. grf_regression_forest price mpg weight length, gen(yhat)}{p_end}

{pstd}With variance estimation and fewer trees for speed{p_end}
{phang2}{cmd:. grf_regression_forest price mpg weight length, gen(yhat2) estimatevariance ntrees(500) replace}{p_end}

{pstd}Without honesty{p_end}
{phang2}{cmd:. grf_regression_forest price mpg weight, gen(yhat3) nohonesty replace}{p_end}

{marker results}{...}
{title:Stored results}

{pstd}
{cmd:grf_regression_forest} stores the following in {cmd:e()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(n_trees)}}number of trees{p_end}
{synopt:{cmd:e(seed)}}random-number seed{p_end}
{synopt:{cmd:e(mtry)}}number of split candidates{p_end}
{synopt:{cmd:e(min_node)}}minimum node size{p_end}
{synopt:{cmd:e(alpha)}}split imbalance bound{p_end}
{synopt:{cmd:e(honesty)}}1 if honest, 0 otherwise{p_end}

{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:grf_regression_forest}{p_end}
{synopt:{cmd:e(forest_type)}}{cmd:regression}{p_end}
{synopt:{cmd:e(depvar)}}name of dependent variable{p_end}
{synopt:{cmd:e(indepvars)}}names of predictor variables{p_end}
{synopt:{cmd:e(predict_var)}}name of prediction variable{p_end}
{synopt:{cmd:e(variance_var)}}name of variance variable (if {cmd:estimatevariance}){p_end}

{marker references}{...}
{title:References}

{pstd}
Athey, S., J. Tibshirani, and S. Wager. 2019.
Generalized Random Forests.
{it:Annals of Statistics} 47(2): 1148{c -}1178.

{pstd}
See also {cmd:{help grf_causal_forest}}, {cmd:{help grf_quantile_forest}}.

{title:Author}

{pstd}
GRF Stata plugin. Wraps the grf C++ library (grf-labs/grf, v2.5.0, GPL-3.0).
