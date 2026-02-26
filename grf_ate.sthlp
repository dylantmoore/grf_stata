{smcl}
{* *! version 0.1.0}{...}
{viewerjumpto "Syntax" "grf_ate##syntax"}{...}
{viewerjumpto "Description" "grf_ate##description"}{...}
{viewerjumpto "Examples" "grf_ate##examples"}{...}
{viewerjumpto "Stored results" "grf_ate##results"}{...}

{title:Title}

{phang}
{bf:grf_ate} {hline 2} Average Treatment Effect (AIPW) from a causal forest

{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmd:grf_ate}
{ifin}

{marker description}{...}
{title:Description}

{pstd}
{cmd:grf_ate} computes the doubly-robust Augmented Inverse Probability
Weighted (AIPW) Average Treatment Effect using results stored in {cmd:e()}
by a prior {helpb grf_causal_forest} estimation.

{pstd}
The doubly-robust score for each observation is:

{pmore}
DR_i = tau_hat_i + (W_i - W_hat_i) / Var(W - W_hat) * (Y_i - Y_hat_i - tau_hat_i * (W_i - W_hat_i))

{pstd}
The ATE is the sample mean of these scores.  Standard errors use the
normal approximation: SE = sd(DR) / sqrt(N).

{pstd}
{cmd:grf_ate} requires that the causal forest was estimated with nuisance
variables stored (the default behavior of {cmd:grf_causal_forest}).

{marker examples}{...}
{title:Examples}

{phang2}{cmd:. grf_causal_forest y w x1 x2 x3, gen(tau)}{p_end}
{phang2}{cmd:. grf_ate}{p_end}

{marker results}{...}
{title:Stored results}

{pstd}
{cmd:grf_ate} stores the following in {cmd:r()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:r(ate)}}ATE point estimate{p_end}
{synopt:{cmd:r(se)}}standard error{p_end}
{synopt:{cmd:r(ci_lower)}}lower bound of 95% confidence interval{p_end}
{synopt:{cmd:r(ci_upper)}}upper bound of 95% confidence interval{p_end}
{synopt:{cmd:r(pvalue)}}two-sided p-value{p_end}
{synopt:{cmd:r(N)}}number of observations{p_end}

{title:References}

{pstd}
Athey, S., J. Tibshirani, and S. Wager. 2019.
Generalized Random Forests. {it:Annals of Statistics} 47(2): 1148-1178.

{title:Author}

{pstd}
GRF Stata plugin. Wraps the grf C++ library (grf-labs/grf, v2.5.0, GPL-3.0).
