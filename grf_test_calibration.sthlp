{smcl}
{* *! version 0.1.0}{...}
{viewerjumpto "Syntax" "grf_test_calibration##syntax"}{...}
{viewerjumpto "Description" "grf_test_calibration##description"}{...}
{viewerjumpto "Examples" "grf_test_calibration##examples"}{...}
{viewerjumpto "Stored results" "grf_test_calibration##results"}{...}

{title:Title}

{phang}
{bf:grf_test_calibration} {hline 2} Calibration test for a causal forest

{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmd:grf_test_calibration}
{ifin}

{marker description}{...}
{title:Description}

{pstd}
{cmd:grf_test_calibration} implements the calibration test of
Chernozhukov, Demirer, Duflo, and Fernandez-Val (2020) for a causal
forest previously estimated by {helpb grf_causal_forest}.

{pstd}
The test regresses outcome residuals (Y - Y_hat) on two terms
constructed from treatment residuals (W - W_hat), with no constant and
HC3 robust standard errors:

{phang2}1. {bf:mean.forest.prediction}: (W - W_hat) * 1.  Tests whether the
forest's average treatment effect is non-zero.{p_end}

{phang2}2. {bf:differential.forest.prediction}: (W - W_hat) * (tau_hat -
mean(tau_hat)).  Tests heterogeneity calibration.  A coefficient near 1
indicates well-calibrated heterogeneity.{p_end}

{pstd}
{bf:Note:} This command runs {cmd:regress} internally, which overwrites
{cmd:e()} results from the causal forest.  Run {cmd:grf_ate} and
{cmd:grf_best_linear_projection} before this command if you need them.

{marker examples}{...}
{title:Examples}

{phang2}{cmd:. grf_causal_forest y w x1 x2 x3, gen(tau)}{p_end}
{phang2}{cmd:. grf_test_calibration}{p_end}

{marker results}{...}
{title:Stored results}

{pstd}
{cmd:grf_test_calibration} stores the following in {cmd:r()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:r(b_mean)}}coefficient on mean.forest.prediction{p_end}
{synopt:{cmd:r(se_mean)}}standard error of mean coefficient{p_end}
{synopt:{cmd:r(t_mean)}}t-statistic of mean coefficient{p_end}
{synopt:{cmd:r(p_mean)}}p-value of mean coefficient{p_end}
{synopt:{cmd:r(b_diff)}}coefficient on differential.forest.prediction{p_end}
{synopt:{cmd:r(se_diff)}}standard error of differential coefficient{p_end}
{synopt:{cmd:r(t_diff)}}t-statistic of differential coefficient{p_end}
{synopt:{cmd:r(p_diff)}}p-value of differential coefficient{p_end}
{synopt:{cmd:r(N)}}number of observations{p_end}

{title:References}

{pstd}
Chernozhukov, V., M. Demirer, E. Duflo, and I. Fernandez-Val. 2020.
Generic Machine Learning Inference on Heterogeneous Treatment Effects in
Randomized Experiments. {it:arXiv:1712.04802}.

{pstd}
Athey, S., J. Tibshirani, and S. Wager. 2019.
Generalized Random Forests. {it:Annals of Statistics} 47(2): 1148-1178.

{title:Author}

{pstd}
GRF Stata plugin. Wraps the grf C++ library (grf-labs/grf, v2.5.0, GPL-3.0).
