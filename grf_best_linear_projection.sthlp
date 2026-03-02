{smcl}
{* *! version 0.3.0}{...}
{viewerjumpto "Syntax" "grf_best_linear_projection##syntax"}{...}
{viewerjumpto "Description" "grf_best_linear_projection##description"}{...}
{viewerjumpto "Options" "grf_best_linear_projection##options"}{...}
{viewerjumpto "Examples" "grf_best_linear_projection##examples"}{...}
{viewerjumpto "Stored results" "grf_best_linear_projection##results"}{...}

{title:Title}

{phang}
{bf:grf_best_linear_projection} {hline 2} Best linear projection of CATE onto covariates

{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmd:grf_best_linear_projection}
[{varlist}]
{ifin}
[{cmd:,} {it:options}]

{synoptset 30 tabbed}{...}
{synopthdr}
{synoptline}
{synopt:{opt vcovtype(string)}}variance estimator: {cmd:HC0}, {cmd:HC1}, {cmd:HC2}, or {cmd:HC3} (default){p_end}
{synopt:{opt targetsample(string)}}target sample: {cmd:all} (default), {cmd:treated}, {cmd:control}, or {cmd:overlap}{p_end}
{synopt:{opt debiasingweights(varname)}}optional debiasing weights applied to DR scores prior to projection{p_end}
{synoptline}

{marker description}{...}
{title:Description}

{pstd}
{cmd:grf_best_linear_projection} projects doubly-robust scores onto
covariates via linear regression. This mirrors R's
{cmd:best_linear_projection()} behavior and supports both robust
variance selection ({cmd:vcovtype()}) and target-sample weighting
({cmd:targetsample()}).

{pstd}
If {varlist} is omitted, covariates from the prior
{cmd:grf_causal_forest} fit ({cmd:e(indepvars)}) are used.

{pstd}
{bf:Note:} This is an {cmd:eclass} command. Running it replaces
previous {cmd:e()} results from forest estimation.

{marker options}{...}
{title:Options}

{phang}
{opt vcovtype(string)} sets the variance estimator. Allowed values are
{cmd:HC0}, {cmd:HC1}, {cmd:HC2}, and {cmd:HC3}. Default is {cmd:HC3}.

{phang}
{opt targetsample(string)} chooses the target population for projection
weights. Allowed values are {cmd:all}, {cmd:treated}, {cmd:control}, and
{cmd:overlap}. Default is {cmd:all}.

{phang}
{opt debiasingweights(varname)} multiplies DR scores by a user-supplied
weight prior to projection.

{marker examples}{...}
{title:Examples}

{pstd}Project onto all forest covariates:{p_end}
{phang2}{cmd:. grf_causal_forest y w x1 x2 x3 x4 x5, gen(tau)}{p_end}
{phang2}{cmd:. grf_best_linear_projection}{p_end}

{pstd}Project onto a subset with HC1:{p_end}
{phang2}{cmd:. grf_best_linear_projection x1 x3, vcovtype(HC1)}{p_end}

{pstd}Target treated sample:{p_end}
{phang2}{cmd:. grf_best_linear_projection x1 x2, targetsample(treated)}{p_end}

{marker results}{...}
{title:Stored results}

{pstd}
{cmd:grf_best_linear_projection} is an {cmd:eclass} command and stores
standard regression outputs plus GRF-specific metadata, including:

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(F)}}F-statistic{p_end}

{p2col 5 24 28 2: Macros}{p_end}
{synopt:{cmd:e(target_sample)}}target sample used in weighting{p_end}
{synopt:{cmd:e(vcov_type)}}variance estimator{p_end}

{title:References}

{pstd}
Semenova, V. and V. Chernozhukov. 2021.
Debiased Machine Learning of Conditional Average Treatment Effects and
Other Causal Functions. {it:The Econometrics Journal} 24(2): 264-289.

{pstd}
Athey, S., J. Tibshirani, and S. Wager. 2019.
Generalized Random Forests. {it:Annals of Statistics} 47(2): 1148-1178.
