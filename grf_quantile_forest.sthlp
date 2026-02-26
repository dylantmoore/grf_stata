{smcl}
{* *! version 0.1.0}{...}
{viewerjumpto "Syntax" "grf_quantile_forest##syntax"}{...}
{viewerjumpto "Description" "grf_quantile_forest##description"}{...}
{viewerjumpto "Options" "grf_quantile_forest##options"}{...}
{viewerjumpto "Examples" "grf_quantile_forest##examples"}{...}
{viewerjumpto "Stored results" "grf_quantile_forest##results"}{...}
{viewerjumpto "References" "grf_quantile_forest##references"}{...}

{title:Title}

{phang}
{bf:grf_quantile_forest} {hline 2} Generalized Random Forest for conditional quantiles

{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmd:grf_quantile_forest}
{depvar}
{indepvars}
{ifin}{cmd:,}
{opt gen:erate(stub)}
[{it:options}]

{synoptset 30 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Required}
{synopt:{opt gen:erate(stub)}}stub for output variables; creates {it:stub}{cmd:_q}{it:NN} for each quantile{p_end}

{syntab:Quantiles}
{synopt:{opt q:uantiles(numlist)}}quantiles to estimate; default is {cmd:quantiles(0.1 0.5 0.9)}{p_end}

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

{syntab:Output}
{synopt:{opt replace}}overwrite existing output variables{p_end}
{synoptline}

{marker description}{...}
{title:Description}

{pstd}
{cmd:grf_quantile_forest} estimates conditional quantile functions
Q_tau(Y | X = x) using a generalized random forest. For each requested
quantile tau and each observation, it produces an out-of-bag estimate of
the tau-th conditional quantile of the outcome distribution.

{pstd}
Output variables are named {it:stub}{cmd:_q}{it:NN} where {it:NN} =
tau * 100. For example, with {cmd:generate(pred)} and the default
quantiles {cmd:(0.1 0.5 0.9)}, the command creates {cmd:pred_q10},
{cmd:pred_q50}, and {cmd:pred_q90}.

{pstd}
A single forest is fit once and used to predict all requested quantiles
simultaneously.

{marker options}{...}
{title:Options}

{dlgtab:Required}

{phang}
{opt generate(stub)} specifies the variable-name stub.  One variable is
created per quantile, named {it:stub}{cmd:_q}{it:NN}.

{dlgtab:Quantiles}

{phang}
{opt quantiles(numlist)} specifies the quantiles to estimate. Each value
must be strictly between 0 and 1. Default is {cmd:0.1 0.5 0.9}.

{dlgtab:Forest}

{phang}
{opt ntrees(#)} number of trees. Default is 2000.

{phang}
{opt seed(#)} random-number seed. Default is 42.

{phang}
{opt mtry(#)} candidate variables per split. Default is 0
(= floor(sqrt(p)) + 20).

{phang}
{opt minnodesize(#)} minimum leaf size. Default is 5.

{phang}
{opt samplefrac(#)} subsample fraction per tree. Default is 0.5.

{phang}
{opt numthreads(#)} number of threads. Default is 0 (all cores).

{dlgtab:Honesty}

{phang}
{opt honesty} | {opt nohonesty} controls honest estimation. Default is
honest.

{phang}
{opt honestyfrac(#)} fraction reserved for honest estimation.
Default is 0.5.

{phang}
{opt honestyprune} | {opt nohonestyprune} controls pruning of empty
honest leaves. Default is to prune.

{dlgtab:Tuning}

{phang}
{opt alpha(#)} maximum split imbalance. Default is 0.05.

{phang}
{opt imbalancepenalty(#)} penalty on imbalanced splits. Default is 0.0.

{dlgtab:Output}

{phang}
{opt replace} allows overwriting existing output variables.

{marker examples}{...}
{title:Examples}

{pstd}Setup{p_end}
{phang2}{cmd:. sysuse auto, clear}{p_end}

{pstd}Default quantiles (10th, 50th, 90th){p_end}
{phang2}{cmd:. grf_quantile_forest price mpg weight length, gen(phat)}{p_end}
{phang2}{cmd:. list phat_q10 phat_q50 phat_q90 in 1/5}{p_end}

{pstd}Custom quantiles: quartiles{p_end}
{phang2}{cmd:. grf_quantile_forest price mpg weight length, gen(qf) quantiles(0.25 0.5 0.75) replace}{p_end}

{pstd}Single median{p_end}
{phang2}{cmd:. grf_quantile_forest price mpg weight, gen(med) quantiles(0.5) ntrees(500)}{p_end}

{marker results}{...}
{title:Stored results}

{pstd}
{cmd:grf_quantile_forest} stores the following in {cmd:e()}:

{synoptset 22 tabbed}{...}
{p2col 5 22 26 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(n_trees)}}number of trees{p_end}
{synopt:{cmd:e(seed)}}random-number seed{p_end}
{synopt:{cmd:e(mtry)}}number of split candidates{p_end}
{synopt:{cmd:e(min_node)}}minimum node size{p_end}
{synopt:{cmd:e(alpha)}}split imbalance bound{p_end}
{synopt:{cmd:e(honesty)}}1 if honest, 0 otherwise{p_end}
{synopt:{cmd:e(n_quantiles)}}number of quantiles estimated{p_end}

{p2col 5 22 26 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:grf_quantile_forest}{p_end}
{synopt:{cmd:e(forest_type)}}{cmd:quantile}{p_end}
{synopt:{cmd:e(depvar)}}name of dependent variable{p_end}
{synopt:{cmd:e(indepvars)}}names of predictor variables{p_end}
{synopt:{cmd:e(quantiles)}}quantile values estimated{p_end}
{synopt:{cmd:e(predict_vars)}}names of all output variables{p_end}

{marker references}{...}
{title:References}

{pstd}
Athey, S., J. Tibshirani, and S. Wager. 2019.
Generalized Random Forests.
{it:Annals of Statistics} 47(2): 1148{c -}1178.

{pstd}
Meinshausen, N. 2006.
Quantile Regression Forests.
{it:Journal of Machine Learning Research} 7: 983{c -}999.

{pstd}
See also {cmd:{help grf_regression_forest}}, {cmd:{help grf_causal_forest}}.

{title:Author}

{pstd}
GRF Stata plugin. Wraps the grf C++ library (grf-labs/grf, v2.5.0, GPL-3.0).
