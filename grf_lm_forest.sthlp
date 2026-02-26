{smcl}
{* *! version 0.1.0}{...}
{viewerjumpto "Syntax" "grf_lm_forest##syntax"}{...}
{viewerjumpto "Description" "grf_lm_forest##description"}{...}
{viewerjumpto "Options" "grf_lm_forest##options"}{...}
{viewerjumpto "Examples" "grf_lm_forest##examples"}{...}
{viewerjumpto "Stored results" "grf_lm_forest##results"}{...}
{viewerjumpto "References" "grf_lm_forest##references"}{...}

{title:Title}

{phang}
{bf:grf_lm_forest} {hline 2} LM Forest for conditional linear models

{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmd:grf_lm_forest}
{depvar}
{it:regvars}
{ifin}{cmd:,}
{opt gen:erate(stub)}
{opt xvars(varlist)}
[{it:options}]

{synoptset 30 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Required}
{synopt:{opt gen:erate(stub)}}stub for output variables; creates {it:stub}{bf:_1}, {it:stub}{bf:_2}, ..., {it:stub}{bf:_K}{p_end}
{synopt:{opt xvars(varlist)}}covariates X used for splitting{p_end}

{syntab:Forest}
{synopt:{opt ntr:ees(#)}}number of trees; default is {cmd:ntrees(2000)}{p_end}
{synopt:{opt seed(#)}}random-number seed; default is {cmd:seed(42)}{p_end}
{synopt:{opt mtry(#)}}variables tried at each split; default is {cmd:mtry(0)} (= sqrt(p)){p_end}
{synopt:{opt minn:odesize(#)}}minimum leaf size; default is {cmd:minnodesize(5)}{p_end}
{synopt:{opt sample:frac(#)}}fraction of observations per tree; default is {cmd:samplefrac(0.5)}{p_end}
{synopt:{opt nuis:ancetrees(#)}}trees for nuisance models; default is {cmd:nuisancetrees(500)}{p_end}
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
{synopt:{opt stab:ilizesplits}}stabilize treatment splits{p_end}
{synopt:{opt nostab:ilizesplits}}disable split stabilization (the default){p_end}

{syntab:Output}
{synopt:{opt replace}}overwrite existing output variables{p_end}
{synoptline}

{marker description}{...}
{title:Description}

{pstd}
{cmd:grf_lm_forest} estimates the locally heterogeneous coefficients in a
conditional linear model:

{p 8 8 2}
Y = c(x) + h_1(x) * W_1 + h_2(x) * W_2 + ... + h_K(x) * W_K

{pstd}
where Y is the outcome ({depvar}), W_1, ..., W_K are the regressors
({it:regvars}), and X ({opt xvars()}) are the covariates that drive
heterogeneity. The functions h_k(x) capture how the association between
each W_k and Y varies across the covariate space.

{pstd}
The command runs a multi-step procedure:

{phang2}1. Fit a regression forest of Y on X to obtain E[Y|X].{p_end}
{phang2}2. Fit a regression forest of each W_k on X to obtain E[W_k|X].{p_end}
{phang2}3. Center the outcome and regressors by subtracting their conditional means.{p_end}
{phang2}4. Fit a multi-arm causal forest on the centered data to estimate h_1(x), ..., h_K(x).{p_end}

{pstd}
This orthogonalization follows Robinson (1988) and removes the
confounding influence of X. The nuisance forests use {opt nuisancetrees(#)}
trees. The final forest uses {opt ntrees(#)} trees.

{pstd}
The output variables {it:stub}{bf:_1}, ..., {it:stub}{bf:_K} contain
out-of-bag estimates of h_1(x), ..., h_K(x), one variable per regressor
in {it:regvars}, in the order they appear.

{marker options}{...}
{title:Options}

{dlgtab:Required}

{phang}
{opt generate(stub)} specifies the stub for output variables. For K
regressors in {it:regvars}, the command creates {it:stub}{bf:_1},
{it:stub}{bf:_2}, ..., {it:stub}{bf:_K}, each containing out-of-bag
estimates of the corresponding coefficient h_k(x). Observations outside
the {cmd:if}/{cmd:in} sample are set to missing.

{phang}
{opt xvars(varlist)} specifies the covariates X used for splitting in the
forest. These variables determine along which dimensions the coefficients
h_k(x) are allowed to vary.

{dlgtab:Forest}

{phang}
{opt ntrees(#)} sets the number of trees for the final multi-arm causal
forest. Default is 2000.

{phang}
{opt seed(#)} sets the random-number seed for reproducibility. Default is 42.

{phang}
{opt mtry(#)} sets the number of candidate variables considered at each
split. The default {cmd:mtry(0)} uses floor(sqrt(p)) + 20 where p is the
number of splitting covariates in {opt xvars()}.

{phang}
{opt minnodesize(#)} sets the minimum number of observations in each tree
leaf. Default is 5.

{phang}
{opt samplefrac(#)} sets the fraction of observations drawn (without
replacement) for each tree. Default is 0.5.

{phang}
{opt nuisancetrees(#)} sets the number of trees for the nuisance forests
that estimate E[Y|X] and E[W_k|X]. Default is 500.

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

{phang}
{opt stabilizesplits} | {opt nostabilizesplits} controls whether
treatment-split criteria are stabilized. Off by default.

{dlgtab:Output}

{phang}
{opt replace} allows overwriting existing variables with the same names as
the output variables.

{marker examples}{...}
{title:Examples}

{pstd}Setup: simulate a conditional linear model{p_end}
{phang2}{cmd:. clear}{p_end}
{phang2}{cmd:. set obs 2000}{p_end}
{phang2}{cmd:. set seed 12345}{p_end}
{phang2}{cmd:. gen x1 = rnormal()}{p_end}
{phang2}{cmd:. gen x2 = rnormal()}{p_end}
{phang2}{cmd:. gen w1 = rnormal()}{p_end}
{phang2}{cmd:. gen w2 = rnormal()}{p_end}
{phang2}{cmd:. gen y = 1 + x1*w1 + (1+x2)*w2 + rnormal()}{p_end}

{pstd}Estimate conditional coefficients{p_end}
{phang2}{cmd:. grf_lm_forest y w1 w2, gen(beta) xvars(x1 x2)}{p_end}

{pstd}
This creates {cmd:beta_1} (estimates of h_1(x), the coefficient on w1)
and {cmd:beta_2} (estimates of h_2(x), the coefficient on w2).

{pstd}With more trees and stabilized splits{p_end}
{phang2}{cmd:. grf_lm_forest y w1 w2, gen(beta) xvars(x1 x2) ntrees(4000) stabilizesplits replace}{p_end}

{pstd}Single regressor{p_end}
{phang2}{cmd:. grf_lm_forest y w1, gen(coef) xvars(x1 x2)}{p_end}

{marker results}{...}
{title:Stored results}

{pstd}
{cmd:grf_lm_forest} stores the following in {cmd:e()}:

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(n_trees)}}number of trees{p_end}
{synopt:{cmd:e(seed)}}random-number seed{p_end}
{synopt:{cmd:e(n_regressors)}}number of regressors (K){p_end}

{p2col 5 24 28 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:grf_lm_forest}{p_end}
{synopt:{cmd:e(forest_type)}}{cmd:lm_forest}{p_end}
{synopt:{cmd:e(depvar)}}name of outcome variable{p_end}
{synopt:{cmd:e(regvars)}}names of regressor variables (W_1, ..., W_K){p_end}
{synopt:{cmd:e(indepvars)}}names of splitting covariates (X){p_end}
{synopt:{cmd:e(predict_var)}}stub used for output variables{p_end}

{marker references}{...}
{title:References}

{pstd}
Athey, S., J. Tibshirani, and S. Wager. 2019.
Generalized Random Forests.
{it:Annals of Statistics} 47(2): 1148{c -}1178.

{pstd}
Zeileis, A., T. Hothorn, and K. Hornik. 2008.
Model-based Recursive Partitioning.
{it:Journal of Computational and Graphical Statistics} 17(2): 492{c -}514.

{pstd}
See also {cmd:{help grf_causal_forest}}, {cmd:{help grf_multi_arm_causal_forest}},
{cmd:{help grf_regression_forest}}.

{title:Author}

{pstd}
GRF Stata plugin. Wraps the grf C++ library (grf-labs/grf, v2.5.0, GPL-3.0).
