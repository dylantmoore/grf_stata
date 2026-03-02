{smcl}
{* *! version 0.1.0}{...}
{viewerjumpto "Syntax" "grf_rate##syntax"}{...}
{viewerjumpto "Description" "grf_rate##description"}{...}
{viewerjumpto "Options" "grf_rate##options"}{...}
{viewerjumpto "Examples" "grf_rate##examples"}{...}
{viewerjumpto "Stored results" "grf_rate##results"}{...}

{title:Title}

{phang}
{bf:grf_rate} {hline 2} Rank-Weighted Average Treatment Effect (RATE)

{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmd:grf_rate}
{it:priorities_var}
{ifin}
[{cmd:,} {it:options}]

{synoptset 28 tabbed}{...}
{synopthdr}
{synoptline}
{synopt:{opt target(string)}}target metric: {cmd:AUTOC} (default) or {cmd:QINI}{p_end}
{synopt:{opt quantiles(numlist)}}quantile grid for TOC integration; default {cmd:0.1(0.1)1.0}{p_end}
{synopt:{opt bootstrap(#)}}number of bootstrap replications; default {cmd:200}{p_end}
{synopt:{opt catevar(varname)}}variable containing CATE estimates; reads from {cmd:e()} if omitted{p_end}
{synopt:{opt subset(varname)}}binary/indicator subset mask; 0 excludes observation{p_end}
{synopt:{opt compliancescore(varname)}}optional compliance-score weights for IV-style RATE{p_end}
{synopt:{opt debiasingweights(varname)}}optional debiasing weights applied to score variable{p_end}
{synopt:{opt seed(#)}}random number seed; default {cmd:-1} (no seed){p_end}
{synoptline}

{marker description}{...}
{title:Description}

{pstd}
{cmd:grf_rate} computes the Rank-Weighted Average Treatment Effect, a
metric for evaluating treatment prioritization rules.  It measures
whether units ranked highly by {it:priorities_var} have larger
treatment effects than average.

{pstd}
Two targets are available:

{phang2}{cmd:AUTOC}: Area Under the Targeting Operator Characteristic.
Integral of TOC(q) over q in [0,1].{p_end}

{phang2}{cmd:QINI}: Integral of q*TOC(q) over q in [0,1].  Upweights
the right tail of the TOC curve.{p_end}

{pstd}
Standard errors are computed via Poisson bootstrap.

{pstd}
By default, the CATE variable is read from {cmd:e(predict_var)} left by
{helpb grf_causal_forest}.  Use {opt catevar()} to override.

{marker options}{...}
{title:Options}

{phang}
{opt target(string)} selects {cmd:AUTOC} or {cmd:QINI}.  Default is
{cmd:AUTOC}.

{phang}
{opt quantiles(numlist)} specifies the quantile grid used for trapezoidal
integration of the TOC curve.  Values must be in (0,1).  Default is
{cmd:0.1 0.2 ... 1.0}.

{phang}
{opt bootstrap(#)} number of Poisson bootstrap replications for standard
error estimation.  Minimum is {cmd:2}.  Default is {cmd:200}.

{phang}
{opt catevar(varname)} specifies a variable containing CATE estimates.
If omitted, {cmd:grf_rate} reads {cmd:e(predict_var)} from a prior
{cmd:grf_causal_forest} estimation.

{phang}
{opt subset(varname)} restricts RATE computation to observations where
{it:varname} is nonzero. This mirrors R's {cmd:subset} behavior.

{phang}
{opt compliancescore(varname)} multiplies the score variable by the supplied
compliance score before ranking/integration.

{phang}
{opt debiasingweights(varname)} multiplies the score variable by the supplied
weights before RATE integration.

{phang}
{opt seed(#)} sets the random number seed.  Default {cmd:-1} does not
set a seed.

{marker examples}{...}
{title:Examples}

{pstd}Evaluate prioritization using estimated CATEs:{p_end}

{phang2}{cmd:. grf_causal_forest y w x1 x2 x3, gen(tau)}{p_end}
{phang2}{cmd:. grf_rate tau}{p_end}

{pstd}Use QINI target with a custom priorities variable:{p_end}

{phang2}{cmd:. grf_rate my_score, target(QINI) catevar(tau) bootstrap(500)}{p_end}

{pstd}Use a subset mask equivalent to {cmd:if}: {p_end}

{phang2}{cmd:. gen keep = (x1 > 0)}{p_end}
{phang2}{cmd:. grf_rate tau, subset(keep)}{p_end}

{marker results}{...}
{title:Stored results}

{pstd}
{cmd:grf_rate} stores the following in {cmd:r()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:r(estimate)}}RATE point estimate{p_end}
{synopt:{cmd:r(std_err)}}bootstrap standard error{p_end}
{synopt:{cmd:r(z_stat)}}z-statistic (estimate / std_err){p_end}
{synopt:{cmd:r(p_value)}}two-sided p-value{p_end}
{synopt:{cmd:r(n)}}number of observations{p_end}
{synopt:{cmd:r(n_bootstrap)}}number of valid bootstrap replications{p_end}

{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:r(target)}}{cmd:AUTOC} or {cmd:QINI}{p_end}
{synopt:{cmd:r(priorities)}}name of the priorities variable{p_end}
{synopt:{cmd:r(catevar)}}name of the CATE variable{p_end}
{synopt:{cmd:r(subset_var)}}name of subset variable, if specified{p_end}
{synopt:{cmd:r(compliance_score_var)}}name of compliance score variable, if specified{p_end}

{title:References}

{pstd}
Yadlowsky, S., S. Fleming, N. Shah, E. Brunskill, and S. Wager. 2021.
Evaluating Treatment Prioritization Rules via Rank-Weighted Average
Treatment Effects. {it:arXiv:2111.07966}.

{pstd}
Athey, S., J. Tibshirani, and S. Wager. 2019.
Generalized Random Forests. {it:Annals of Statistics} 47(2): 1148-1178.

{title:Author}

{pstd}
GRF Stata plugin. Wraps the grf C++ library (grf-labs/grf, v2.5.0, GPL-3.0).
