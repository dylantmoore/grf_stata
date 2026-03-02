{smcl}
{* *! version 0.2.0}{...}
{viewerjumpto "Syntax" "grf_average_partial_effect##syntax"}{...}
{viewerjumpto "Description" "grf_average_partial_effect##description"}{...}

{title:Title}

{phang}
{bf:grf_average_partial_effect} {hline 2} Average partial effect from a fitted causal forest

{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmd:grf_average_partial_effect}
{ifin}
[{cmd:,}
{opt debiasingweights(varname)}
{opt numtreesvariance(#)}
{opt nocalibrate}]

{marker description}{...}
{title:Description}

{pstd}
{cmd:grf_average_partial_effect} computes an AIPW-style average partial effect
using nuisance quantities from the most recent {cmd:grf_causal_forest} fit.

{pstd}
{bf:Deprecation parity note:} R {bf:grf} has deprecated
{cmd:average_partial_effect()} in favor of
{cmd:average_treatment_effect()} with continuous treatment.
Stata retains {cmd:grf_average_partial_effect} for backward compatibility.

{pstd}
Returns include estimate, standard error, confidence interval, and p-value in
{cmd:r()}.
