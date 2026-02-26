{smcl}
{* *! version 0.1.0}{...}
{viewerjumpto "Description" "grf##description"}{...}
{viewerjumpto "Forest types" "grf##forests"}{...}
{viewerjumpto "Post-estimation" "grf##postestimation"}{...}
{viewerjumpto "Quick start" "grf##quickstart"}{...}
{viewerjumpto "Installation" "grf##installation"}{...}
{viewerjumpto "References" "grf##references"}{...}

{title:Title}

{phang}
{bf:GRF} {hline 2} Generalized Random Forests for Stata

{marker description}{...}
{title:Description}

{pstd}
{bf:grf_stata} is a Stata implementation of the
{browse "https://github.com/grf-labs/grf":grf} R package (v2.5.0) by Athey,
Tibshirani, and Wager.  It wraps the same C++ backend used by the R package,
compiled as a native Stata plugin for full performance.  No external
dependencies are required{hline 2}the grf C++ library is statically linked
into the plugin binary.

{pstd}
The package provides 12 forest types for non-parametric estimation and
7 post-estimation commands for inference, diagnostics, and prediction.
All commands support {cmd:if}/{cmd:in} sample restrictions, store results
in {cmd:e()} or {cmd:r()}, and run multithreaded by default.

{marker forests}{...}
{title:Forest types}

{p2colset 5 40 42 2}{...}
{p2col:Command}Description{p_end}
{p2line}
{p2col:{helpb grf_regression_forest}}Regression forest for conditional mean estimation{p_end}
{p2col:{helpb grf_causal_forest}}Causal forest for heterogeneous treatment effects (CATE){p_end}
{p2col:{helpb grf_quantile_forest}}Quantile forest for conditional quantile estimation{p_end}
{p2col:{helpb grf_instrumental_forest}}Instrumental forest for LATE with instruments{p_end}
{p2col:{helpb grf_probability_forest}}Probability forest for class probability estimation{p_end}
{p2col:{helpb grf_survival_forest}}Survival forest for conditional survival functions{p_end}
{p2col:{helpb grf_causal_survival_forest}}Causal survival forest for treatment effects on survival{p_end}
{p2col:{helpb grf_multi_arm_causal_forest}}Multi-arm causal forest for multi-valued treatments{p_end}
{p2col:{helpb grf_multi_regression_forest}}Multi-outcome regression forest{p_end}
{p2col:{helpb grf_ll_regression_forest}}Local linear regression forest{p_end}
{p2col:{helpb grf_boosted_regression_forest}}Boosted regression forest (iterative residual fitting){p_end}
{p2col:{helpb grf_lm_forest}}Linear model forest (conditional linear coefficients){p_end}
{p2line}

{marker postestimation}{...}
{title:Post-estimation commands}

{p2colset 5 40 42 2}{...}
{p2col:Command}Description{p_end}
{p2line}
{p2col:{helpb grf_ate}}Average treatment effect (AIPW doubly-robust estimator){p_end}
{p2col:{helpb grf_best_linear_projection}}Project CATEs onto covariates (BLP){p_end}
{p2col:{helpb grf_test_calibration}}Calibration test for forest predictions{p_end}
{p2col:{helpb grf_variable_importance}}Variable importance via weighted split frequencies{p_end}
{p2col:{helpb grf_rate}}Rank-weighted average treatment effect (AUTOC/QINI){p_end}
{p2col:{helpb grf_predict}}Predict on new (out-of-sample) data{p_end}
{p2col:{helpb grf_tune}}Cross-validation tuning of forest parameters{p_end}
{p2line}

{marker quickstart}{...}
{title:Quick start}

{pstd}Simulate data with heterogeneous treatment effects{p_end}
{phang2}{cmd:. clear}{p_end}
{phang2}{cmd:. set obs 2000}{p_end}
{phang2}{cmd:. set seed 12345}{p_end}
{phang2}{cmd:. gen x1 = rnormal()}{p_end}
{phang2}{cmd:. gen x2 = rnormal()}{p_end}
{phang2}{cmd:. gen w = rbinomial(1, normal(x1))}{p_end}
{phang2}{cmd:. gen y = 2*x1 + w*(1 + x2) + rnormal()}{p_end}

{pstd}Estimate conditional average treatment effects{p_end}
{phang2}{cmd:. grf_causal_forest y w x1 x2, gen(tau_hat)}{p_end}

{pstd}Average treatment effect{p_end}
{phang2}{cmd:. grf_ate}{p_end}

{pstd}Best linear projection of CATEs onto covariates{p_end}
{phang2}{cmd:. grf_best_linear_projection x1 x2}{p_end}

{pstd}Calibration test{p_end}
{phang2}{cmd:. grf_test_calibration}{p_end}

{pstd}Variable importance{p_end}
{phang2}{cmd:. grf_variable_importance x1 x2}{p_end}

{pstd}Rank-weighted average treatment effect{p_end}
{phang2}{cmd:. grf_rate tau_hat}{p_end}

{marker installation}{...}
{title:Installation}

{pstd}Install from GitHub:{p_end}
{phang2}{cmd:. net install grf_stata, from("https://raw.githubusercontent.com/dylantmoore/grf_stata/main") replace}{p_end}

{pstd}
The package includes pre-compiled plugin binaries for macOS (ARM64 and x86_64),
Linux (x86_64), and Windows (x86_64).  Stata 14.0 or later is required.

{marker references}{...}
{title:References}

{pstd}
Athey, S., J. Tibshirani, and S. Wager. 2019.
Generalized Random Forests.
{it:Annals of Statistics} 47(2): 1148{c -}1178.

{pstd}
Wager, S. and S. Athey. 2018.
Estimation and Inference of Heterogeneous Treatment Effects using Random Forests.
{it:Journal of the American Statistical Association} 113(523): 1228{c -}1242.

{pstd}
Source code and documentation:{break}
{browse "https://github.com/dylantmoore/grf_stata"}

{title:Author}

{pstd}
Dylan Moore{break}
GRF Stata plugin. Wraps the grf C++ library
({browse "https://github.com/grf-labs/grf":grf-labs/grf}, v2.5.0, GPL-3.0)
by Susan Athey, Julie Tibshirani, Stefan Wager, and contributors.
