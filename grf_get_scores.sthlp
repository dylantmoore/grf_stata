{smcl}
{* *! version 0.3.0}{...}
{viewerjumpto "Syntax" "grf_get_scores##syntax"}{...}
{viewerjumpto "Description" "grf_get_scores##description"}{...}
{viewerjumpto "Stored results" "grf_get_scores##results"}{...}

{title:Title}

{phang}
{bf:grf_get_scores} {hline 2} Extract doubly-robust scores from fitted GRF causal models

{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmd:grf_get_scores}{cmd:,} {opt gen:erate(newvar)} [{opt replace}]

{marker description}{...}
{title:Description}

{pstd}
{cmd:grf_get_scores} computes score variables from the most recent fitted model.
Supported forest types are:

{phang}{cmd:grf_causal_forest} (AIPW DR score){p_end}
{phang}{cmd:grf_instrumental_forest} (AIPW-style DR score){p_end}
{phang}{cmd:grf_multi_arm_causal_forest} (one DR score variable per arm){p_end}
{phang}{cmd:grf_causal_survival_forest} (DR score from stored nuisance moments){p_end}

{pstd}
For causal-survival forests, scores follow the upstream construction:
{cmd:score_i = tau_i + psi_i / Vhat_i},
where {cmd:psi_i = numer_i - denom_i * tau_i}.
For binary treatment, {cmd:Vhat_i = W.hat_i * (1 - W.hat_i)};
otherwise a denominator-based fallback is used.

{marker results}{...}
{title:Stored results}

{pstd}
{cmd:grf_get_scores} stores in {cmd:r()}:

{synoptset 28 tabbed}{...}
{p2col 5 28 32 2: Scalars}{p_end}
{synopt:{cmd:r(N)}}observations used{p_end}
{synopt:{cmd:r(mean)}}mean score{p_end}
{synopt:{cmd:r(sd)}}score SD{p_end}
{synopt:{cmd:r(se)}}SE of score mean{p_end}
{synopt:{cmd:r(min)}}minimum score{p_end}
{synopt:{cmd:r(max)}}maximum score{p_end}
{synopt:{cmd:r(denom_mean)}}mean denominator nuisance (causal-survival only){p_end}
{synopt:{cmd:r(vhat_mean)}}mean effective Vhat used in score correction (causal-survival only){p_end}

{p2col 5 28 32 2: Macros}{p_end}
{synopt:{cmd:r(generate)}}generated score variable(s){p_end}
{synopt:{cmd:r(forest_type)}}forest type used{p_end}
{synopt:{cmd:r(vhat_source)}}Vhat source used by causal-survival scoring{p_end}
