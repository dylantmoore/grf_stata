{smcl}
{* *! version 0.1.0}{...}
{viewerjumpto "Syntax" "grf_tune##syntax"}{...}
{viewerjumpto "Description" "grf_tune##description"}{...}
{viewerjumpto "Options" "grf_tune##options"}{...}
{viewerjumpto "Examples" "grf_tune##examples"}{...}
{viewerjumpto "Stored results" "grf_tune##results"}{...}

{title:Title}

{phang}
{bf:grf_tune} {hline 2} Cross-validation tuning for GRF forests via random search

{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmd:grf_tune}
{depvar} [{it:treatvar}] {indepvars}
{ifin}{cmd:,}
{opt foresttype(string)}
[{it:options}]

{synoptset 28 tabbed}{...}
{synopthdr}
{synoptline}
{p2coldent:* {opt foresttype(string)}}forest type: {cmd:regression}, {cmd:causal}, or {cmd:quantile}{p_end}
{synopt:{opt numreps(#)}}number of random search candidates; default {cmd:50}{p_end}
{synopt:{opt tunetrees(#)}}trees per candidate forest; default {cmd:200}{p_end}
{synopt:{opt seed(#)}}random number seed; default {cmd:42}{p_end}
{synopt:{opt numthreads(#)}}number of threads; default {cmd:0} (all available){p_end}
{synopt:{opt honesty}}use honesty (default){p_end}
{synopt:{opt nohonesty}}disable honesty{p_end}
{synopt:{opt honestyprune}}prune honest leaves (default){p_end}
{synopt:{opt nohonestyprune}}disable honesty pruning{p_end}
{synopt:{opt stabilizesplits}}stabilize splits for causal forests (default){p_end}
{synopt:{opt nostabilizesplits}}disable split stabilization{p_end}
{synoptline}
{p 4 6 2}* {opt foresttype()} is required.{p_end}
{p 4 6 2}For {cmd:causal} forests, the varlist is {depvar} {it:treatvar} {indepvars}.{p_end}

{marker description}{...}
{title:Description}

{pstd}
{cmd:grf_tune} performs random-search hyperparameter tuning for GRF forests.
It evaluates {opt numreps()} random parameter configurations, each fitted with
a small forest of {opt tunetrees()} trees, and selects the configuration that
minimizes out-of-bag MSE.

{pstd}
The tuned parameters are: {cmd:mtry}, {cmd:min_node_size},
{cmd:sample_fraction}, {cmd:honesty_fraction}, {cmd:alpha}, and
{cmd:imbalance_penalty}.  Results are returned in {cmd:r()} for use
in a subsequent forest estimation call.

{pstd}
For causal forests, nuisance models (Y~X and W~X) are fitted once before
the search loop, and the causal forest MSE is computed on centered outcomes.

{marker options}{...}
{title:Options}

{phang}
{opt foresttype(string)} specifies the forest type.  Must be
{cmd:regression}, {cmd:causal}, or {cmd:quantile}.  Required.

{phang}
{opt numreps(#)} number of random candidate parameter vectors to evaluate.
Default is {cmd:50}.

{phang}
{opt tunetrees(#)} number of trees grown per candidate.  Smaller values
speed up tuning at the cost of noisier OOB estimates.  Default is {cmd:200}.
Minimum is {cmd:10}.

{phang}
{opt seed(#)} random number seed for reproducibility.  Default is {cmd:42}.

{phang}
{opt numthreads(#)} threads for the C plugin.  Default {cmd:0} uses all
available cores.

{phang}
{opt honesty} / {opt nohonesty} enable or disable honesty.  Default is on.

{phang}
{opt honestyprune} / {opt nohonestyprune} enable or disable honesty leaf
pruning.  Default is on.

{phang}
{opt stabilizesplits} / {opt nostabilizesplits} enable or disable split
stabilization (causal forests only).  Default is on.

{marker examples}{...}
{title:Examples}

{pstd}Tune a regression forest:{p_end}

{phang2}{cmd:. grf_tune price mpg weight length, foresttype(regression) numreps(100)}{p_end}
{phang2}{cmd:. grf_regression_forest price mpg weight length, gen(yhat) mtry(`r(best_mtry)') minnodesize(`r(best_min_node_size)')}{p_end}

{pstd}Tune a causal forest:{p_end}

{phang2}{cmd:. grf_tune y w x1 x2 x3, foresttype(causal) numreps(100) seed(123)}{p_end}
{phang2}{cmd:. local sf = r(best_sample_fraction)}{p_end}
{phang2}{cmd:. grf_causal_forest y w x1 x2 x3, gen(tau) samplefrac(`sf')}{p_end}

{marker results}{...}
{title:Stored results}

{pstd}
{cmd:grf_tune} stores the following in {cmd:r()}:

{synoptset 30 tabbed}{...}
{p2col 5 30 34 2: Scalars}{p_end}
{synopt:{cmd:r(best_mtry)}}optimal number of variables tried at each split{p_end}
{synopt:{cmd:r(best_min_node_size)}}optimal minimum node size{p_end}
{synopt:{cmd:r(best_sample_fraction)}}optimal sample fraction{p_end}
{synopt:{cmd:r(best_honesty_fraction)}}optimal honesty fraction{p_end}
{synopt:{cmd:r(best_alpha)}}optimal alpha{p_end}
{synopt:{cmd:r(best_imbalance_penalty)}}optimal imbalance penalty{p_end}
{synopt:{cmd:r(best_mse)}}OOB MSE of the best configuration{p_end}
{synopt:{cmd:r(n_reps)}}number of search candidates evaluated{p_end}
{synopt:{cmd:r(N)}}number of observations used{p_end}

{p2col 5 30 34 2: Macros}{p_end}
{synopt:{cmd:r(forest_type)}}forest type specified{p_end}

{title:References}

{pstd}
Athey, S., J. Tibshirani, and S. Wager. 2019.
Generalized Random Forests. {it:Annals of Statistics} 47(2): 1148-1178.

{title:Author}

{pstd}
GRF Stata plugin. Wraps the grf C++ library (grf-labs/grf, v2.5.0, GPL-3.0).
