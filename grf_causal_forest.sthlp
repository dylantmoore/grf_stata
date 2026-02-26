{smcl}
{* *! version 0.1.0}{...}
{viewerjumpto "Syntax" "grf_causal_forest##syntax"}{...}
{viewerjumpto "Description" "grf_causal_forest##description"}{...}
{viewerjumpto "Options" "grf_causal_forest##options"}{...}
{viewerjumpto "Examples" "grf_causal_forest##examples"}{...}
{viewerjumpto "Stored results" "grf_causal_forest##results"}{...}
{viewerjumpto "References" "grf_causal_forest##references"}{...}

{title:Title}

{phang}
{bf:grf_causal_forest} {hline 2} Generalized Random Forest for heterogeneous treatment effects

{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmd:grf_causal_forest}
{depvar}
{it:treatvar}
{indepvars}
{ifin}{cmd:,}
{opt gen:erate(newvar)}
[{it:options}]

{synoptset 30 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Required}
{synopt:{opt gen:erate(newvar)}}name of variable to store CATE predictions{p_end}

{syntab:Forest}
{synopt:{opt ntr:ees(#)}}number of trees; default is {cmd:ntrees(2000)}{p_end}
{synopt:{opt seed(#)}}random-number seed; default is {cmd:seed(42)}{p_end}
{synopt:{opt mtry(#)}}variables tried at each split; default is {cmd:mtry(0)} (= sqrt(p)){p_end}
{synopt:{opt minn:odesize(#)}}minimum leaf size; default is {cmd:minnodesize(5)}{p_end}
{synopt:{opt sample:frac(#)}}fraction of observations per tree; default is {cmd:samplefrac(0.5)}{p_end}
{synopt:{opt nuis:ancetrees(#)}}trees for nuisance models (Y.hat, W.hat); default is {cmd:nuisancetrees(500)}{p_end}
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
{synopt:{opt stab:ilizesplits}}stabilize treatment splits (the default){p_end}
{synopt:{opt nostab:ilizesplits}}disable split stabilization{p_end}

{syntab:Variance estimation}
{synopt:{opt est:imatevariance}}compute variance of CATE predictions{p_end}
{synopt:{opt varg:enerate(newvar)}}name for variance variable; default is {it:generate}{cmd:_var}{p_end}
{synopt:{opt cig:roupsize(#)}}cluster size for variance; default is {cmd:cigroupsize(1)}; forced to 2+ with {cmd:estimatevariance}{p_end}

{syntab:Nuisance output}
{synopt:{opt yhatg:enerate(newvar)}}save Y.hat (outcome model predictions){p_end}
{synopt:{opt whatg:enerate(newvar)}}save W.hat (propensity model predictions){p_end}

{syntab:Output}
{synopt:{opt replace}}overwrite existing output variables{p_end}
{synoptline}

{marker description}{...}
{title:Description}

{pstd}
{cmd:grf_causal_forest} estimates conditional average treatment effects
(CATEs) using a generalized random forest. For each observation it
produces tau(x) = E[Y(1) - Y(0) | X = x], the expected effect of
treatment for units with covariates x.

{pstd}
The command runs a three-step nuisance pipeline automatically:

{phang2}1. Fit a regression forest of Y on X to obtain Y.hat.{p_end}
{phang2}2. Fit a regression forest of W on X to obtain W.hat (propensity scores).{p_end}
{phang2}3. Center the outcome and treatment (Robinson, 1988), then fit the causal forest on centered data.{p_end}

{pstd}
This orthogonalization removes confounding bias. The nuisance models
Y.hat and W.hat are always saved as {cmd:_grf_yhat} and {cmd:_grf_what}
for use by post-estimation commands such as {cmd:{help grf_ate}} and
{cmd:{help grf_test_calibration}}.

{pstd}
The treatment variable may be binary or continuous. For binary treatment
the CATE is the unit-level treatment effect; for continuous treatment it
is the partial effect.

{marker options}{...}
{title:Options}

{dlgtab:Required}

{phang}
{opt generate(newvar)} creates {it:newvar} containing out-of-bag CATE
predictions tau.hat(x) for each observation. Observations outside the
estimation sample are set to missing.

{dlgtab:Forest}

{phang}
{opt ntrees(#)} sets the number of trees. Default is 2000. This number is
used for all three forests (two nuisance + one causal).

{phang}
{opt seed(#)} sets the random-number seed. Default is 42.

{phang}
{opt mtry(#)} number of candidate variables at each split. Default is 0
(= floor(sqrt(p)) + 20).

{phang}
{opt minnodesize(#)} minimum leaf size. Default is 5.

{phang}
{opt samplefrac(#)} fraction of observations per tree. Default is 0.5.

{phang}
{opt numthreads(#)} number of threads. Default is 0 (all cores).

{dlgtab:Honesty}

{phang}
{opt honesty} | {opt nohonesty} controls honest estimation. With honesty
(the default), split-selection and leaf-estimation use disjoint subsamples
within each tree.

{phang}
{opt honestyfrac(#)} fraction reserved for honest estimation. Default is 0.5.

{phang}
{opt honestyprune} | {opt nohonestyprune} controls pruning of empty
honest leaves. Default is to prune.

{dlgtab:Tuning}

{phang}
{opt alpha(#)} maximum split imbalance. Default is 0.05.

{phang}
{opt imbalancepenalty(#)} penalty on imbalanced splits. Default is 0.0.

{phang}
{opt stabilizesplits} | {opt nostabilizesplits} controls whether
treatment-split criteria are stabilized by the propensity score. Enabled
by default. Stabilization reduces variance when the propensity score
is close to 0 or 1.

{dlgtab:Variance estimation}

{phang}
{opt estimatevariance} requests variance estimates for each CATE
prediction. When specified, {cmd:cigroupsize} is forced to at least 2.

{phang}
{opt vargenerate(newvar)} names the variance variable.  Default is
{it:generate}{cmd:_var}.

{phang}
{opt cigroupsize(#)} cluster size for variance estimation. Default is 1;
forced to 2+ with {cmd:estimatevariance}.

{dlgtab:Nuisance output}

{phang}
{opt yhatgenerate(newvar)} saves the outcome-model predictions Y.hat in a
user-named variable (in addition to the internal {cmd:_grf_yhat}).

{phang}
{opt whatgenerate(newvar)} saves the propensity-model predictions W.hat in
a user-named variable (in addition to the internal {cmd:_grf_what}).

{dlgtab:Output}

{phang}
{opt replace} allows overwriting existing output variables.

{marker examples}{...}
{title:Examples}

{pstd}Setup: simulate treatment effects{p_end}
{phang2}{cmd:. clear}{p_end}
{phang2}{cmd:. set obs 1000}{p_end}
{phang2}{cmd:. set seed 12345}{p_end}
{phang2}{cmd:. gen x1 = rnormal()}{p_end}
{phang2}{cmd:. gen x2 = rnormal()}{p_end}
{phang2}{cmd:. gen w = rbinomial(1, normal(x1))}{p_end}
{phang2}{cmd:. gen y = 2*x1 + w*(1 + x2) + rnormal()}{p_end}

{pstd}Estimate CATEs{p_end}
{phang2}{cmd:. grf_causal_forest y w x1 x2, gen(tau_hat)}{p_end}

{pstd}With variance estimation and nuisance output{p_end}
{phang2}{cmd:. grf_causal_forest y w x1 x2, gen(tau_hat) estimatevariance yhatgenerate(yhat) whatgenerate(what) replace}{p_end}

{pstd}Compute ATE from stored results{p_end}
{phang2}{cmd:. display "ATE = " %6.3f e(ate) " (se = " %6.3f e(ate_se) ")"}{p_end}

{marker results}{...}
{title:Stored results}

{pstd}
{cmd:grf_causal_forest} stores the following in {cmd:e()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(n_trees)}}number of trees{p_end}
{synopt:{cmd:e(seed)}}random-number seed{p_end}
{synopt:{cmd:e(mtry)}}number of split candidates{p_end}
{synopt:{cmd:e(min_node)}}minimum node size{p_end}
{synopt:{cmd:e(alpha)}}split imbalance bound{p_end}
{synopt:{cmd:e(honesty)}}1 if honest, 0 otherwise{p_end}
{synopt:{cmd:e(stabilize)}}1 if splits stabilized, 0 otherwise{p_end}
{synopt:{cmd:e(ate)}}average treatment effect (mean of tau.hat){p_end}
{synopt:{cmd:e(ate_se)}}standard error of the ATE{p_end}

{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:grf_causal_forest}{p_end}
{synopt:{cmd:e(forest_type)}}{cmd:causal}{p_end}
{synopt:{cmd:e(depvar)}}name of outcome variable{p_end}
{synopt:{cmd:e(treatvar)}}name of treatment variable{p_end}
{synopt:{cmd:e(indepvars)}}names of predictor variables{p_end}
{synopt:{cmd:e(predict_var)}}name of CATE prediction variable{p_end}
{synopt:{cmd:e(variance_var)}}name of variance variable (if {cmd:estimatevariance}){p_end}
{synopt:{cmd:e(yhat_var)}}name of Y.hat variable ({cmd:_grf_yhat}){p_end}
{synopt:{cmd:e(what_var)}}name of W.hat variable ({cmd:_grf_what}){p_end}

{marker references}{...}
{title:References}

{pstd}
Athey, S., J. Tibshirani, and S. Wager. 2019.
Generalized Random Forests.
{it:Annals of Statistics} 47(2): 1148{c -}1178.

{pstd}
Robinson, P. M. 1988.
Root-N-consistent semiparametric regression.
{it:Econometrica} 56(4): 931{c -}954.

{pstd}
See also {cmd:{help grf_regression_forest}}, {cmd:{help grf_quantile_forest}},
{cmd:{help grf_ate}}, {cmd:{help grf_test_calibration}}.

{title:Author}

{pstd}
GRF Stata plugin. Wraps the grf C++ library (grf-labs/grf, v2.5.0, GPL-3.0).
