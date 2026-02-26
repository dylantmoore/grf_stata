{smcl}
{* *! version 0.1.0}{...}
{viewerjumpto "Syntax" "grf_variable_importance##syntax"}{...}
{viewerjumpto "Description" "grf_variable_importance##description"}{...}
{viewerjumpto "Options" "grf_variable_importance##options"}{...}
{viewerjumpto "Examples" "grf_variable_importance##examples"}{...}
{viewerjumpto "Stored results" "grf_variable_importance##results"}{...}

{title:Title}

{phang}
{bf:grf_variable_importance} {hline 2} Split-frequency variable importance from a GRF regression forest

{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmd:grf_variable_importance}
{depvar} {indepvars}
{ifin}
[{cmd:,} {it:options}]

{synoptset 20 tabbed}{...}
{synopthdr}
{synoptline}
{synopt:{opt ntrees(#)}}number of trees; default {cmd:2000}{p_end}
{synopt:{opt seed(#)}}random number seed; default {cmd:42}{p_end}
{synopt:{opt maxdepth(#)}}maximum tree depth; default {cmd:4}{p_end}
{synoptline}

{marker description}{...}
{title:Description}

{pstd}
{cmd:grf_variable_importance} computes split-frequency weighted variable
importance scores by fitting a regression forest via the GRF C plugin.
Each score measures how often a variable is used for splitting, weighted
by the depth of the split.  Scores sum to 1.

{pstd}
This is a standalone command that fits its own forest; it does not require
a prior estimation step.

{marker options}{...}
{title:Options}

{phang}
{opt ntrees(#)} number of trees to grow.  Default is {cmd:2000}.

{phang}
{opt seed(#)} random number seed.  Default is {cmd:42}.

{phang}
{opt maxdepth(#)} maximum depth of each tree.  Default is {cmd:4}.

{marker examples}{...}
{title:Examples}

{phang2}{cmd:. sysuse auto, clear}{p_end}
{phang2}{cmd:. grf_variable_importance price mpg weight length turn}{p_end}

{phang2}{cmd:. grf_variable_importance y x1 x2 x3 x4, ntrees(5000) maxdepth(6)}{p_end}

{marker results}{...}
{title:Stored results}

{pstd}
{cmd:grf_variable_importance} stores the following in {cmd:r()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:r(N)}}number of observations{p_end}

{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:r(importance)}}1 x p matrix of importance scores (columns named by variable){p_end}

{title:References}

{pstd}
Athey, S., J. Tibshirani, and S. Wager. 2019.
Generalized Random Forests. {it:Annals of Statistics} 47(2): 1148-1178.

{title:Author}

{pstd}
GRF Stata plugin. Wraps the grf C++ library (grf-labs/grf, v2.5.0, GPL-3.0).
