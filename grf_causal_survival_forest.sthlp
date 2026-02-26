{smcl}
{* *! version 0.1.0}{...}
{viewerjumpto "Syntax" "grf_causal_survival_forest##syntax"}{...}
{viewerjumpto "Description" "grf_causal_survival_forest##description"}{...}
{viewerjumpto "Options" "grf_causal_survival_forest##options"}{...}
{viewerjumpto "Examples" "grf_causal_survival_forest##examples"}{...}
{viewerjumpto "Stored results" "grf_causal_survival_forest##results"}{...}
{viewerjumpto "References" "grf_causal_survival_forest##references"}{...}

{title:Title}

{phang}
{bf:grf_causal_survival_forest} {hline 2} Causal survival forest for estimating
heterogeneous treatment effects with survival outcomes

{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmd:grf_causal_survival_forest}
{it:timevar} {it:statusvar} {it:treatvar} {indepvars}
{ifin}{cmd:,}
{opt gen:erate(newvar)}
[{it:options}]

{synoptset 30 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Main}
{synopt:{opt gen:erate(newvar)}}name for CATE predictions{p_end}
{synopt:{opt replace}}overwrite existing variable{p_end}

{syntab:Forest}
{synopt:{opt ntrees(#)}}number of trees; default {bf:2000}{p_end}
{synopt:{opt seed(#)}}random seed; default {bf:42}{p_end}
{synopt:{opt mtry(#)}}variables per split; default {bf:ceil(sqrt(p))}{p_end}
{synopt:{opt minnodesize(#)}}minimum node size; default {bf:15}{p_end}
{synopt:{opt samplefrac(#)}}sample fraction; default {bf:0.5}{p_end}
{synopt:{opt honesty}/{opt nohonesty}}use honesty splitting; default {bf:honesty}{p_end}
{synopt:{opt honestyfrac(#)}}honesty fraction; default {bf:0.5}{p_end}
{synopt:{opt alpha(#)}}balance parameter; default {bf:0.05}{p_end}
{synopt:{opt imbalancepenalty(#)}}imbalance penalty; default {bf:0.0}{p_end}
{synopt:{opt numthreads(#)}}threads; default {bf:0} (auto){p_end}

{syntab:Causal survival}
{synopt:{opt stab:ilizesplits}/{opt nostab:ilizesplits}}stabilize splits; default {bf:stabilize}{p_end}
{synopt:{opt numer(varname)}}pre-computed numerator (advanced){p_end}
{synopt:{opt denom(varname)}}pre-computed denominator (advanced){p_end}
{synopt:{opt hor:izon(#)}}time horizon for restricted mean survival{p_end}
{synopt:{opt tar:get(#)}}target type: 1=RMST (default), 2=survival probability{p_end}

{syntab:Variance}
{synopt:{opt est:imatevariance}}compute variance estimates{p_end}
{synopt:{opt var:generate(newvar)}}store variance estimates{p_end}
{synoptline}

{marker description}{...}
{title:Description}

{pstd}
{cmd:grf_causal_survival_forest} estimates heterogeneous treatment effects
for survival outcomes using the causal survival forest of Cui et al. (2023).
It implements the {cmd:causal_survival_forest()} function from R's {bf:grf}
package.

{pstd}
The command implements a multi-step nuisance estimation pipeline:
(1) survival forest for censoring model,
(2) regression forests for nuisance parameters,
(3) causal survival forest on the centered data.

{pstd}
{it:timevar} is the survival/failure time variable,
{it:statusvar} is the censoring indicator (1=event, 0=censored),
and {it:treatvar} is the binary treatment variable.

{marker options}{...}
{title:Options}

{dlgtab:Main}

{phang}
{opt generate(newvar)} specifies the name for the CATE prediction variable.

{phang}
{opt numer(varname)} and {opt denom(varname)} provide pre-computed nuisance
estimates for advanced users who want to supply their own censoring model.
If omitted, the command estimates them internally.

{phang}
{opt horizon(#)} sets the time horizon for restricted mean survival time.
If 0 (default), uses the maximum observed failure time.

{marker examples}{...}
{title:Examples}

{pstd}Estimate causal survival forest:{p_end}
{phang}{cmd:. grf_causal_survival_forest time status treatment x1 x2 x3, gen(cate) ntrees(2000) seed(42)}{p_end}

{marker results}{...}
{title:Stored results}

{pstd}
{cmd:grf_causal_survival_forest} stores the following in {cmd:e()}:

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(n_trees)}}number of trees{p_end}
{synopt:{cmd:e(seed)}}random seed{p_end}

{p2col 5 24 28 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}grf_causal_survival_forest{p_end}
{synopt:{cmd:e(forest_type)}}causal_survival{p_end}
{synopt:{cmd:e(predict_var)}}name of prediction variable{p_end}

{marker references}{...}
{title:References}

{pstd}
Cui, Y., M. R. Kosorok, E. Sverdrup, S. Wager, and R. Zhu. 2023.
Estimating Heterogeneous Treatment Effects with Right-Censored Data via
Causal Survival Forests. {it:Journal of the Royal Statistical Society, Series B} 85(2): 179-211.

{pstd}
Athey, S., J. Tibshirani, and S. Wager. 2019.
Generalized Random Forests. {it:Annals of Statistics} 47(2): 1148-1178.

{title:Author}

{pstd}
GRF Stata plugin. Wraps the grf C++ library (grf-labs/grf, v2.5.0, GPL-3.0).
