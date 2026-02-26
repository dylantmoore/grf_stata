{smcl}
{* *! version 0.1.0}{...}
{viewerjumpto "Syntax" "grf_predict##syntax"}{...}
{viewerjumpto "Description" "grf_predict##description"}{...}
{viewerjumpto "Options" "grf_predict##options"}{...}
{viewerjumpto "Examples" "grf_predict##examples"}{...}
{viewerjumpto "Stored results" "grf_predict##results"}{...}

{title:Title}

{phang}
{bf:grf_predict} {hline 2} Predict on new data using a previously estimated GRF forest

{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmd:grf_predict}{cmd:,}
{opth gen:erate(newvar)}
[{it:options}]

{synoptset 28 tabbed}{...}
{synopthdr}
{synoptline}
{p2coldent:* {opth gen:erate(newvar)}}name of variable to store predictions{p_end}
{synopt:{opt replace}}replace {it:newvar} if it already exists{p_end}
{synopt:{opt numthreads(#)}}number of threads; default {cmd:0} (all available){p_end}
{synoptline}
{p 4 6 2}* {opt generate()} is required.{p_end}
{p 4 6 2}Forest parameters (trees, mtry, sample fraction, honesty, alpha, etc.)
are read automatically from {cmd:e()} stored by the prior estimation command.{p_end}

{marker description}{...}
{title:Description}

{pstd}
{cmd:grf_predict} generates predictions on new (test) observations using
a forest previously estimated by {helpb grf_regression_forest} or
{helpb grf_causal_forest}.

{pstd}
The workflow is:

{phang2}1. Estimate a forest on training data.{p_end}
{phang2}2. {cmd:append} the test data so that test observations follow the training observations.{p_end}
{phang2}3. Run {cmd:grf_predict} to write predictions for the test rows.{p_end}

{pstd}
The command reads all forest parameters from {cmd:e()} stored by the
prior estimation command.  Training observations are identified as the
first {cmd:e(N)} rows; remaining rows are treated as test data.
Predictions for training rows are set to missing.

{pstd}
For causal forests, {cmd:grf_predict} automatically re-fits nuisance
models (Y~X and W~X) on the training data to center outcomes before
predicting CATEs on test observations.

{marker options}{...}
{title:Options}

{phang}
{opth gen:erate(newvar)} specifies the name of the new variable to hold
predictions.  For regression and causal forests this is a single variable.
Required.

{phang}
{opt replace} allows {cmd:grf_predict} to drop and recreate {it:newvar}
if it already exists.

{phang}
{opt numthreads(#)} number of threads for the C plugin.  Default {cmd:0}
uses all available cores.

{marker examples}{...}
{title:Examples}

{pstd}Predict with a regression forest:{p_end}

{phang2}{cmd:. sysuse auto, clear}{p_end}
{phang2}{cmd:. grf_regression_forest price mpg weight length, gen(yhat) ntrees(500)}{p_end}
{phang2}{cmd:. append using test_cars}{p_end}
{phang2}{cmd:. grf_predict, gen(yhat_new)}{p_end}

{pstd}Predict CATEs with a causal forest:{p_end}

{phang2}{cmd:. grf_causal_forest y w x1 x2 x3, gen(tau) ntrees(2000)}{p_end}
{phang2}{cmd:. append using newdata}{p_end}
{phang2}{cmd:. grf_predict, gen(tau_new)}{p_end}

{marker results}{...}
{title:Stored results}

{pstd}
{cmd:grf_predict} stores the following in {cmd:r()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:r(N_train)}}number of training observations{p_end}
{synopt:{cmd:r(N_test)}}number of test observations predicted{p_end}
{synopt:{cmd:r(mean)}}mean prediction (regression forest){p_end}
{synopt:{cmd:r(sd)}}SD of predictions (regression forest){p_end}
{synopt:{cmd:r(mean_cate)}}mean CATE prediction (causal forest){p_end}
{synopt:{cmd:r(sd_cate)}}SD of CATE predictions (causal forest){p_end}

{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:r(forest_type)}}type of forest ({cmd:regression}, {cmd:causal}){p_end}
{synopt:{cmd:r(predict_var)}}name of the generated prediction variable{p_end}

{title:References}

{pstd}
Athey, S., J. Tibshirani, and S. Wager. 2019.
Generalized Random Forests. {it:Annals of Statistics} 47(2): 1148-1178.

{title:Author}

{pstd}
GRF Stata plugin. Wraps the grf C++ library (grf-labs/grf, v2.5.0, GPL-3.0).
