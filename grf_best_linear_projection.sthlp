{smcl}
{* *! version 0.2.0}{...}
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

{pstd}
If {varlist} is omitted, the covariates from the prior
{cmd:grf_causal_forest} estimation ({cmd:e(indepvars)}) are used.

{marker description}{...}
{title:Description}

{pstd}
{cmd:grf_best_linear_projection} projects doubly-robust AIPW scores onto
covariates via OLS with HC3 robust standard errors.  This matches R's
{cmd:grf::best_linear_projection()}.

{pstd}
The doubly-robust score for each observation is:

{pmore}
DR_i = tau_hat_i + (W_i - W_hat_i) / Var(W - W_hat) * (Y_i - Y_hat_i - tau_hat_i * (W_i - W_hat_i))

{pstd}
The BLP regression is: DR_i = alpha + beta * X_i + epsilon_i, with HC3
standard errors.  Significant coefficients indicate which covariates are
associated with treatment effect heterogeneity.

{pstd}
{bf:Note:} This is an {cmd:eclass} command.  Running it replaces {cmd:e()}
results from {cmd:grf_causal_forest} with regression output.  Run
{cmd:grf_ate} and {cmd:grf_test_calibration} first if needed.

{marker options}{...}
{title:Options}

{pstd}
{cmd:grf_best_linear_projection} takes no options beyond the optional
{varlist} and {ifin} qualifiers.  If {varlist} is omitted, the original
forest covariates are used.

{marker examples}{...}
{title:Examples}

{pstd}Project onto all forest covariates:{p_end}

{phang2}{cmd:. grf_causal_forest y w x1 x2 x3, gen(tau)}{p_end}
{phang2}{cmd:. grf_best_linear_projection}{p_end}

{pstd}Project onto a subset of covariates:{p_end}

{phang2}{cmd:. grf_best_linear_projection x1 x3}{p_end}

{marker results}{...}
{title:Stored results}

{pstd}
{cmd:grf_best_linear_projection} is an {cmd:eclass} command.  It stores
standard {cmd:regress, vce(hc3)} results in {cmd:e()}, including:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(r2)}}R-squared{p_end}
{synopt:{cmd:e(F)}}F-statistic{p_end}

{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}coefficient vector{p_end}
{synopt:{cmd:e(V)}}HC3 variance-covariance matrix{p_end}

{title:References}

{pstd}
Semenova, V. and V. Chernozhukov. 2021.
Debiased Machine Learning of Conditional Average Treatment Effects and
Other Causal Functions. {it:The Econometrics Journal} 24(2): 264-289.

{pstd}
Athey, S., J. Tibshirani, and S. Wager. 2019.
Generalized Random Forests. {it:Annals of Statistics} 47(2): 1148-1178.

{title:Author}

{pstd}
GRF Stata plugin. Wraps the grf C++ library (grf-labs/grf, v2.5.0, GPL-3.0).
