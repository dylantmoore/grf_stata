{smcl}
{* *! version 0.1.0}{...}
{viewerjumpto "Syntax" "grf_boosted_regression_forest##syntax"}{...}
{viewerjumpto "Description" "grf_boosted_regression_forest##description"}{...}
{viewerjumpto "Options" "grf_boosted_regression_forest##options"}{...}
{viewerjumpto "Examples" "grf_boosted_regression_forest##examples"}{...}
{viewerjumpto "Stored results" "grf_boosted_regression_forest##results"}{...}
{viewerjumpto "References" "grf_boosted_regression_forest##references"}{...}

{title:Title}

{phang}
{bf:grf_boosted_regression_forest} {hline 2} Boosted regression forest

{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmd:grf_boosted_regression_forest}
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
{synopt:{opt ntr:ees(#)}}number of trees per boosting step; default is {cmd:ntrees(2000)}{p_end}
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

{syntab:Boosting}
{synopt:{opt booststeps(#)}}fixed number of boosting steps; default is {cmd:booststeps(0)} (auto-tune){p_end}
{synopt:{opt boostmaxsteps(#)}}maximum steps when auto-tuning; default is {cmd:boostmaxsteps(5)}{p_end}
{synopt:{opt boosterrorreduction(#)}}error reduction threshold for auto-stopping; default is {cmd:boosterrorreduction(0.97)}{p_end}
{synopt:{opt boosttreestune(#)}}trees for cross-validation check; default is {cmd:boosttreestune(10)}{p_end}

{syntab:Output}
{synopt:{opt replace}}overwrite existing output variables{p_end}
{synoptline}

{marker description}{...}
{title:Description}

{pstd}
{cmd:grf_boosted_regression_forest} fits a boosted regression forest for
non-parametric regression of {depvar} on {indepvars}. Rather than training
a single forest, it trains multiple regression forests sequentially, where
each subsequent forest fits the residuals from the accumulated predictions
of all previous forests. The final prediction is the sum of all forests'
predictions.

{pstd}
When {opt booststeps(0)} is specified (the default), the number of boosting
steps is determined automatically by cross-validation. At each candidate
step a small forest with {opt boosttreestune(#)} trees is fit to the
current residuals and its cross-validated error is compared to the previous
step's error. A new step is taken only if the estimated error is below
{opt boosterrorreduction(#)} times the previous step's error. This
continues until the improvement threshold is not met or {opt boostmaxsteps(#)}
is reached.

{pstd}
Boosting can improve predictive accuracy when a single forest underfits,
particularly in settings with strong signal-to-noise ratios or when the
conditional mean function is complex.

{marker options}{...}
{title:Options}

{dlgtab:Required}

{phang}
{opt generate(newvar)} creates {it:newvar} containing out-of-bag
predictions of E[Y|X=x] for each observation in the estimation sample.
Predictions are the sum across all boosting steps. Observations outside the
{cmd:if}/{cmd:in} sample are set to missing.

{dlgtab:Forest}

{phang}
{opt ntrees(#)} sets the number of trees in each boosting step. Default is
2000.

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

{dlgtab:Boosting}

{phang}
{opt booststeps(#)} sets a fixed number of boosting steps. When set to 0
(the default), the number of steps is determined automatically by
cross-validation.

{phang}
{opt boostmaxsteps(#)} sets the maximum number of boosting steps allowed
when auto-tuning ({cmd:booststeps(0)}). Default is 5.

{phang}
{opt boosterrorreduction(#)} sets the error reduction threshold for
auto-stopping. A new boosting step is taken only if its cross-validated
error is below this fraction of the previous step's error. Default is 0.97,
meaning a step must reduce error by at least 3% to continue.

{phang}
{opt boosttreestune(#)} sets the number of trees used in the small
cross-validation forest that evaluates whether an additional boosting step
improves fit. Default is 10.

{dlgtab:Output}

{phang}
{opt replace} allows overwriting an existing variable with the same name as
{it:generate}.

{marker examples}{...}
{title:Examples}

{pstd}Setup{p_end}
{phang2}{cmd:. sysuse auto, clear}{p_end}

{pstd}Boosted regression forest with auto-tuning (default){p_end}
{phang2}{cmd:. grf_boosted_regression_forest price mpg weight length, gen(yhat)}{p_end}

{pstd}Fixed number of boosting steps{p_end}
{phang2}{cmd:. grf_boosted_regression_forest price mpg weight length, gen(yhat2) booststeps(3) replace}{p_end}

{pstd}Auto-tuning with stricter stopping criterion and more max steps{p_end}
{phang2}{cmd:. grf_boosted_regression_forest price mpg weight, gen(yhat3) boostmaxsteps(10) boosterrorreduction(0.95) replace}{p_end}

{pstd}Fewer trees per step for speed{p_end}
{phang2}{cmd:. grf_boosted_regression_forest price mpg weight length, gen(yhat4) ntrees(500) booststeps(5) replace}{p_end}

{marker results}{...}
{title:Stored results}

{pstd}
{cmd:grf_boosted_regression_forest} stores the following in {cmd:e()}:

{synoptset 28 tabbed}{...}
{p2col 5 28 32 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(n_trees)}}number of trees per boosting step{p_end}
{synopt:{cmd:e(seed)}}random-number seed{p_end}
{synopt:{cmd:e(boost_steps)}}number of boosting steps taken{p_end}
{synopt:{cmd:e(boost_max_steps)}}maximum boosting steps allowed{p_end}
{synopt:{cmd:e(boost_error_reduction)}}error reduction threshold{p_end}

{p2col 5 28 32 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:grf_boosted_regression_forest}{p_end}
{synopt:{cmd:e(forest_type)}}{cmd:boosted_regression}{p_end}
{synopt:{cmd:e(depvar)}}name of dependent variable{p_end}
{synopt:{cmd:e(indepvars)}}names of predictor variables{p_end}
{synopt:{cmd:e(predict_var)}}name of prediction variable{p_end}

{marker references}{...}
{title:References}

{pstd}
Athey, S., J. Tibshirani, and S. Wager. 2019.
Generalized Random Forests.
{it:Annals of Statistics} 47(2): 1148{c -}1178.

{pstd}
See also {cmd:{help grf_regression_forest}}, {cmd:{help grf_ll_regression_forest}}.

{title:Author}

{pstd}
GRF Stata plugin. Wraps the grf C++ library (grf-labs/grf, v2.5.0, GPL-3.0).
