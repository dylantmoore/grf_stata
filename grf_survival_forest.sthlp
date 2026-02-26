{smcl}
{* *! version 0.1.0}{...}
{viewerjumpto "Syntax" "grf_survival_forest##syntax"}{...}
{viewerjumpto "Description" "grf_survival_forest##description"}{...}
{viewerjumpto "Options" "grf_survival_forest##options"}{...}
{viewerjumpto "Examples" "grf_survival_forest##examples"}{...}
{viewerjumpto "Stored results" "grf_survival_forest##results"}{...}
{viewerjumpto "References" "grf_survival_forest##references"}{...}

{title:Title}

{phang}
{bf:grf_survival_forest} {hline 2} Survival forest for conditional survival function estimation

{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmd:grf_survival_forest}
{it:timevar}
{it:statusvar}
{it:indepvars}
{ifin}{cmd:,}
{opt gen:erate(stub)}
[{it:options}]

{synoptset 28 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Required}
{synopt:{opt gen:erate(stub)}}stub for output variables; creates {it:stub}{cmd:_s1}, {it:stub}{cmd:_s2}, ...{p_end}

{syntab:Survival}
{synopt:{opt pred:type(#)}}0 = Kaplan-Meier, 1 = Nelson-Aalen; default is {cmd:1}{p_end}
{synopt:{opt nout:put(#)}}number of failure-time columns; default is {cmd:20}{p_end}
{synopt:{opt numf:ailures(#)}}expected number of unique failure times; default {cmd:0} (auto-detect){p_end}

{syntab:Forest tuning}
{synopt:{opt nt:rees(#)}}number of trees; default is {cmd:2000}{p_end}
{synopt:{opt seed(#)}}random seed; default is {cmd:42}{p_end}
{synopt:{opt mtry(#)}}variables to consider at each split; default {cmd:0} (auto){p_end}
{synopt:{opt minn:odesize(#)}}minimum leaf size; default is {cmd:15}{p_end}
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
{synopt:{opt cig:roupsize(#)}}cluster size for splitting; default is {cmd:1}{p_end}

{syntab:Output}
{synopt:{opt re:place}}overwrite existing output variables{p_end}
{synoptline}

{marker description}{...}
{title:Description}

{pstd}
{cmd:grf_survival_forest} estimates conditional survival functions S(t|X) using
the survival forest method of Cui et al. (2023), implemented within the
generalized random forests framework. The command takes a survival time variable
{it:timevar} (must be strictly positive) and a censoring indicator {it:statusvar}
(1 = event, 0 = censored).

{pstd}
The command generates {opt noutput()} variables named {it:stub}{cmd:_s1},
{it:stub}{cmd:_s2}, ..., {it:stub}{cmd:_s}{it:N}, where each variable contains
the estimated survival probability at the corresponding failure time. The
failure times are the unique observed event times from the data (or an
evenly-spaced grid if {opt noutput()} is fewer than the number of unique
failure times).

{marker options}{...}
{title:Options}

{dlgtab:Required}

{phang}
{opt generate(stub)} specifies the stub name for the output variables. The
command creates {it:stub}{cmd:_s1} through {it:stub}{cmd:_s}{it:N}, where N is
set by {opt noutput()}.

{dlgtab:Survival}

{phang}
{opt predtype(#)} selects the survival estimator: {cmd:0} for Kaplan-Meier and
{cmd:1} for Nelson-Aalen (the default). Nelson-Aalen is generally preferred for
its smoothness properties.

{phang}
{opt noutput(#)} sets the number of failure-time grid points at which the
survival curve is evaluated. Default is {cmd:20}.

{phang}
{opt numfailures(#)} provides a hint for the number of unique failure times.
The default {cmd:0} auto-detects from the data.

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
The default is {cmd:15}, larger than other forest types to account for censoring.

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

{pstd}Setup with simulated survival data:{p_end}
{phang2}{cmd:. clear}{p_end}
{phang2}{cmd:. set obs 1000}{p_end}
{phang2}{cmd:. set seed 12345}{p_end}
{phang2}{cmd:. gen x1 = rnormal()}{p_end}
{phang2}{cmd:. gen x2 = rnormal()}{p_end}
{phang2}{cmd:. gen truetime = -ln(runiform()) * exp(-0.5*x1)}{p_end}
{phang2}{cmd:. gen censtime = -ln(runiform()) * 2}{p_end}
{phang2}{cmd:. gen time = min(truetime, censtime)}{p_end}
{phang2}{cmd:. gen byte status = (truetime <= censtime)}{p_end}

{pstd}Basic survival forest (Nelson-Aalen, 20 time points):{p_end}
{phang2}{cmd:. grf_survival_forest time status x1 x2, gen(surv)}{p_end}
{phang2}{cmd:. sum surv_s1 surv_s10 surv_s20}{p_end}

{pstd}Kaplan-Meier estimator with more time points:{p_end}
{phang2}{cmd:. grf_survival_forest time status x1 x2, gen(surv_km) predtype(0) noutput(50)}{p_end}

{pstd}Compare survival curves for two groups:{p_end}
{phang2}{cmd:. grf_survival_forest time status x1 x2, gen(sf) replace}{p_end}
{phang2}{cmd:. sum sf_s10 if x1 > 0}{p_end}
{phang2}{cmd:. sum sf_s10 if x1 <= 0}{p_end}

{marker results}{...}
{title:Stored results}

{pstd}{cmd:grf_survival_forest} stores the following in {cmd:e()}:

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(n_events)}}number of events{p_end}
{synopt:{cmd:e(n_censored)}}number of censored observations{p_end}
{synopt:{cmd:e(n_trees)}}number of trees{p_end}
{synopt:{cmd:e(seed)}}random seed{p_end}
{synopt:{cmd:e(mtry)}}number of split candidates{p_end}
{synopt:{cmd:e(min_node)}}minimum node size{p_end}
{synopt:{cmd:e(alpha)}}alpha parameter{p_end}
{synopt:{cmd:e(honesty)}}1 if honesty enabled, 0 otherwise{p_end}
{synopt:{cmd:e(n_output)}}number of output time columns{p_end}
{synopt:{cmd:e(pred_type)}}prediction type (0 = Nelson-Aalen, 1 = Kaplan-Meier){p_end}

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:grf_survival_forest}{p_end}
{synopt:{cmd:e(forest_type)}}{cmd:survival}{p_end}
{synopt:{cmd:e(timevar)}}survival time variable name{p_end}
{synopt:{cmd:e(statusvar)}}censoring indicator variable name{p_end}
{synopt:{cmd:e(indepvars)}}predictor variable names{p_end}
{synopt:{cmd:e(predict_stub)}}output variable stub name{p_end}

{marker references}{...}
{title:References}

{pstd}
Athey, S., J. Tibshirani, and S. Wager. 2019.
Generalized Random Forests. {it:Annals of Statistics} 47(2): 1148-1178.

{pstd}
Cui, Y., M. Zhu, M. Kosorok, and others. 2023.
Survival Analysis with Generalized Random Forests.
{browse "https://github.com/grf-labs/grf"}

{title:Author}

{pstd}
GRF Stata plugin. Wraps the grf C++ library (grf-labs/grf, v2.5.0, GPL-3.0).
