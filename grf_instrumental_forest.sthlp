{smcl}
{* *! version 0.1.0}{...}
{viewerjumpto "Syntax" "grf_instrumental_forest##syntax"}{...}
{viewerjumpto "Description" "grf_instrumental_forest##description"}{...}
{viewerjumpto "Options" "grf_instrumental_forest##options"}{...}
{viewerjumpto "Examples" "grf_instrumental_forest##examples"}{...}
{viewerjumpto "Stored results" "grf_instrumental_forest##results"}{...}
{viewerjumpto "References" "grf_instrumental_forest##references"}{...}

{title:Title}

{phang}
{bf:grf_instrumental_forest} {hline 2} Instrumental forest for heterogeneous treatment effects with instrumental variables

{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmd:grf_instrumental_forest}
{it:depvar}
{it:treatvar}
{it:instrvar}
{it:indepvars}
{ifin}{cmd:,}
{opt gen:erate(newvar)}
[{it:options}]

{synoptset 28 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Required}
{synopt:{opt gen:erate(newvar)}}variable to store LATE estimates{p_end}

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

{syntab:Instrumental forest}
{synopt:{opt sta:bilizesplits}}stabilize instrument splits{p_end}
{synopt:{opt red:ucedformweight(#)}}weight on reduced-form estimates; default is {cmd:0.0}{p_end}
{synopt:{opt nuis:ancetrees(#)}}trees for each nuisance regression forest; default is {cmd:500}{p_end}

{syntab:Variance estimation}
{synopt:{opt est:imatevariance}}estimate asymptotic variance{p_end}
{synopt:{opt vargen:erate(newvar)}}variable for variance estimates; default is {it:newvar}{cmd:_var}{p_end}
{synopt:{opt cig:roupsize(#)}}cluster size for variance; default is {cmd:1} (auto-set to 2 when {opt estimatevariance} is specified){p_end}

{syntab:Output}
{synopt:{opt re:place}}overwrite existing output variable(s){p_end}
{synoptline}

{marker description}{...}
{title:Description}

{pstd}
{cmd:grf_instrumental_forest} estimates conditional local average treatment effects
(CLATEs) using the instrumental forest method of Athey, Tibshirani, and Wager (2019).
It uses an instrument {it:instrvar} to address endogeneity of the treatment
{it:treatvar} when estimating heterogeneous effects on {it:depvar} as a function of
covariates {it:indepvars}.

{pstd}
The command runs a nuisance pipeline of three regression forests to center the
outcome, treatment, and instrument on their conditional expectations given X,
then fits the instrumental forest on the centered values.

{marker options}{...}
{title:Options}

{dlgtab:Required}

{phang}
{opt generate(newvar)} creates {it:newvar} containing the estimated LATE for
each observation.

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

{dlgtab:Instrumental forest}

{phang}
{opt stabilizesplits} applies stabilization to the instrument-based splits,
which can improve performance when the first stage is weak.

{phang}
{opt reducedformweight(#)} controls the weight placed on the reduced-form
(intent-to-treat) moment relative to the IV moment. A value of 0 (the default)
uses pure IV; a value of 1 uses the reduced form only.

{phang}
{opt nuisancetrees(#)} sets the number of trees in each of the three
nuisance regression forests used to center Y, W, and Z.

{dlgtab:Variance estimation}

{phang}
{opt estimatevariance} computes an estimate of the asymptotic variance of the
LATE predictions using the forest's built-in variance estimator.

{phang}
{opt vargenerate(newvar)} names the variable to store variance estimates.
If omitted, defaults to {it:generate}{cmd:_var}.

{phang}
{opt cigroupsize(#)} sets the cluster group size for the variance estimator.
When {opt estimatevariance} is specified, this is automatically set to at
least 2.

{marker examples}{...}
{title:Examples}

{pstd}Setup with simulated data:{p_end}
{phang2}{cmd:. clear}{p_end}
{phang2}{cmd:. set obs 1000}{p_end}
{phang2}{cmd:. set seed 12345}{p_end}
{phang2}{cmd:. gen x1 = rnormal()}{p_end}
{phang2}{cmd:. gen x2 = rnormal()}{p_end}
{phang2}{cmd:. gen z = (rnormal() > 0)}{p_end}
{phang2}{cmd:. gen w = 0.5*z + 0.3*x1 + rnormal()}{p_end}
{phang2}{cmd:. gen y = 2*w*x1 + x2 + rnormal()}{p_end}

{pstd}Basic instrumental forest:{p_end}
{phang2}{cmd:. grf_instrumental_forest y w z x1 x2, gen(iv_tau)}{p_end}

{pstd}With variance estimation and stabilized splits:{p_end}
{phang2}{cmd:. grf_instrumental_forest y w z x1 x2, gen(iv_tau) estimatevariance stabilizesplits replace}{p_end}

{pstd}Increase reduced-form weight:{p_end}
{phang2}{cmd:. grf_instrumental_forest y w z x1 x2, gen(iv_tau) reducedformweight(0.5) replace}{p_end}

{marker results}{...}
{title:Stored results}

{pstd}{cmd:grf_instrumental_forest} stores the following in {cmd:e()}:

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(n_trees)}}number of trees{p_end}
{synopt:{cmd:e(seed)}}random seed{p_end}
{synopt:{cmd:e(mtry)}}number of split candidates{p_end}
{synopt:{cmd:e(min_node)}}minimum node size{p_end}
{synopt:{cmd:e(alpha)}}alpha parameter{p_end}
{synopt:{cmd:e(honesty)}}1 if honesty enabled, 0 otherwise{p_end}
{synopt:{cmd:e(reduced_form_wt)}}reduced-form weight{p_end}
{synopt:{cmd:e(stabilize_splits)}}1 if splits stabilized, 0 otherwise{p_end}

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:grf_instrumental_forest}{p_end}
{synopt:{cmd:e(forest_type)}}{cmd:instrumental}{p_end}
{synopt:{cmd:e(depvar)}}outcome variable name{p_end}
{synopt:{cmd:e(treatment)}}treatment variable name{p_end}
{synopt:{cmd:e(instrument)}}instrument variable name{p_end}
{synopt:{cmd:e(indepvars)}}predictor variable names{p_end}
{synopt:{cmd:e(predict_var)}}name of LATE prediction variable{p_end}
{synopt:{cmd:e(variance_var)}}name of variance variable (if estimated){p_end}

{marker references}{...}
{title:References}

{pstd}
Athey, S., J. Tibshirani, and S. Wager. 2019.
Generalized Random Forests. {it:Annals of Statistics} 47(2): 1148-1178.

{title:Author}

{pstd}
GRF Stata plugin. Wraps the grf C++ library (grf-labs/grf, v2.5.0, GPL-3.0).
