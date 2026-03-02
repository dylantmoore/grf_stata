{smcl}
{* *! version 0.2.0}{...}
{viewerjumpto "Syntax" "grf_causal_survival_forest##syntax"}{...}
{viewerjumpto "Description" "grf_causal_survival_forest##description"}{...}
{viewerjumpto "Options" "grf_causal_survival_forest##options"}{...}
{viewerjumpto "Examples" "grf_causal_survival_forest##examples"}{...}
{viewerjumpto "Stored results" "grf_causal_survival_forest##results"}{...}
{viewerjumpto "References" "grf_causal_survival_forest##references"}{...}

{title:Title}

{phang}
{bf:grf_causal_survival_forest} {hline 2} Causal survival forest for heterogeneous treatment effects with right-censored outcomes

{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmd:grf_causal_survival_forest}
{it:timevar} {it:statusvar} {it:treatvar} {it:indepvars}
{ifin}{cmd:,}
{opt gen:erate(newvar)}
[{it:options}]

{synoptset 32 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Main}
{synopt:{opt gen:erate(newvar)}}name for CATE predictions{p_end}
{synopt:{opt replace}}overwrite existing output variable{p_end}

{syntab:Forest}
{synopt:{opt ntrees(#)}}number of trees; default {cmd:2000}{p_end}
{synopt:{opt seed(#)}}random seed; default {cmd:42}{p_end}
{synopt:{opt mtry(#)}}variables considered per split; default {cmd:mtry(0)} (auto){p_end}
{synopt:{opt minnodesize(#)}}minimum node size; default {cmd:15}{p_end}
{synopt:{opt samplefrac(#)}}sample fraction per tree; default {cmd:0.5}{p_end}
{synopt:{opt honesty}/{opt nohonesty}}honest splitting; default {cmd:honesty}{p_end}
{synopt:{opt honestyfrac(#)}}honesty split fraction; default {cmd:0.5}{p_end}
{synopt:{opt nohonestyprune}}disable honest-leaf pruning{p_end}
{synopt:{opt alpha(#)}}imbalance bound; default {cmd:0.05}{p_end}
{synopt:{opt imbalancepenalty(#)}}split imbalance penalty; default {cmd:0.0}{p_end}
{synopt:{opt numthreads(#)}}threads; default {cmd:0} (auto){p_end}

{syntab:Causal survival}
{synopt:{opt horizon(#)}}horizon for RMST/survival-probability estimand (0 = median event time){p_end}
{synopt:{opt target(#)}}estimand target: {cmd:1}=RMST (default), {cmd:2}=survival probability{p_end}
{synopt:{opt nostabilizesplits}}disable split stabilization{p_end}
{synopt:{opt numer(varname)}}precomputed IPCW numerator column{p_end}
{synopt:{opt denom(varname)}}precomputed IPCW denominator column{p_end}
{synopt:{opt whatinput(varname)}}user-supplied propensity nuisance {it:W.hat}{p_end}
{synopt:{opt yhatinput(varname)}}user-supplied outcome nuisance {it:Y.hat}{p_end}
{synopt:{opt chatinput(varname)}}user-supplied censoring nuisance {it:C.hat}{p_end}

{syntab:Variance}
{synopt:{opt estimatevariance}}compute variance estimates{p_end}
{synopt:{opt vargenerate(newvar)}}name for variance output (default: {it:generate}{cmd:_var}){p_end}
{synopt:{opt cigroupsize(#)}}ci-group size; forced to >=2 with {cmd:estimatevariance}{p_end}

{syntab:Missingness and weighting}
{synopt:{opt nomia}}disable MIA and use casewise deletion of missing X values{p_end}
{synopt:{opt cluster(varname)}}cluster ID variable{p_end}
{synopt:{opt weights(varname)}}sample weights{p_end}
{synopt:{opt equalizeclusterweights}}reweight to equalize cluster contribution; requires {cmd:cluster()}{p_end}

{syntab:Inline tuning}
{synopt:{opt tuneparameters(string)}}parameters to tune inline{p_end}
{synopt:{opt tunenumtrees(#)}}trees used in tuning stage; default {cmd:200}{p_end}
{synopt:{opt tunenumreps(#)}}tuning repetitions; default {cmd:50}{p_end}
{synoptline}

{marker description}{...}
{title:Description}

{pstd}
{cmd:grf_causal_survival_forest} estimates heterogeneous treatment effects for
survival outcomes using the causal survival forest of Cui et al. (2023).

{pstd}
The command supports two nuisance modes:

{phang}
{bf:Estimated IPCW mode} (default): nuisance quantities are estimated internally
with regression/survival proxies, producing IPCW numerator/denominator columns
used by the causal-survival forest objective.{p_end}

{phang}
{bf:Precomputed mode}: user supplies {cmd:numer()} and {cmd:denom()} directly.
These must be supplied together.{p_end}

{pstd}
Optional nuisance hooks {cmd:whatinput()}, {cmd:yhatinput()}, and
{cmd:chatinput()} let advanced users override parts of the default nuisance
pipeline while keeping the forest fit in Stata.

{marker options}{...}
{title:Options}

{phang}
{opt numer(varname)} and {opt denom(varname)} specify precomputed nuisance
moments. Use either both or neither.

{phang}
{opt whatinput(varname)} supplies propensity nuisance estimates
{it:W.hat = E[W|X]}.

{phang}
{opt yhatinput(varname)} supplies the outcome nuisance used in the IPCW
moment construction.

{phang}
{opt chatinput(varname)} supplies a censoring-survival proxy used in IPCW
weighting; values are clipped internally to a numerically stable range.

{phang}
{opt target(#)} chooses the estimand used when constructing nuisance moments:
{cmd:1} for RMST, {cmd:2} for survival probability at {cmd:horizon()}.

{marker examples}{...}
{title:Examples}

{pstd}Default nuisance estimation{p_end}
{phang2}{cmd:. grf_causal_survival_forest time status w x1 x2 x3, gen(cs_tau)}{p_end}

{pstd}User-supplied nuisance hooks{p_end}
{phang2}{cmd:. grf_causal_survival_forest time status w x1 x2 x3, gen(cs_tau2) horizon(4) whatinput(wh) yhatinput(yh) chatinput(ch)}{p_end}

{pstd}Fully precomputed nuisance moments{p_end}
{phang2}{cmd:. grf_causal_survival_forest time status w x1 x2 x3, gen(cs_tau3) numer(num) denom(den)}{p_end}

{marker results}{...}
{title:Stored results}

{pstd}
{cmd:grf_causal_survival_forest} stores the following in {cmd:e()}:

{synoptset 28 tabbed}{...}
{p2col 5 28 32 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}observations used{p_end}
{synopt:{cmd:e(n_events)}}event count{p_end}
{synopt:{cmd:e(n_censored)}}censored count{p_end}
{synopt:{cmd:e(n_trees)}}number of trees{p_end}
{synopt:{cmd:e(horizon)}}horizon used{p_end}
{synopt:{cmd:e(target)}}target type (1=RMST, 2=survival probability){p_end}
{synopt:{cmd:e(ate)}}mean CATE prediction{p_end}
{synopt:{cmd:e(ate_se)}}SE of mean CATE prediction{p_end}

{p2col 5 28 32 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:grf_causal_survival_forest}{p_end}
{synopt:{cmd:e(forest_type)}}{cmd:causal_survival}{p_end}
{synopt:{cmd:e(predict_var)}}CATE output variable{p_end}
{synopt:{cmd:e(variance_var)}}variance output variable (if requested){p_end}
{synopt:{cmd:e(what_var)}}stored propensity nuisance variable ({cmd:_grf_cs_what}){p_end}
{synopt:{cmd:e(numer_var)}}stored IPCW numerator variable ({cmd:_grf_cs_numer}){p_end}
{synopt:{cmd:e(denom_var)}}stored IPCW denominator variable ({cmd:_grf_cs_denom}){p_end}
{synopt:{cmd:e(yhat_var)}}stored outcome nuisance variable when available ({cmd:_grf_cs_yhat}){p_end}
{synopt:{cmd:e(chat_var)}}stored censoring nuisance variable when available ({cmd:_grf_cs_chat}){p_end}
{synopt:{cmd:e(nuisance_mode)}}{cmd:estimated_ipcw} or {cmd:precomputed_numer_denom}{p_end}

{marker references}{...}
{title:References}

{pstd}
Cui, Y., M. R. Kosorok, E. Sverdrup, S. Wager, and R. Zhu. 2023.
Estimating Heterogeneous Treatment Effects with Right-Censored Data via
Causal Survival Forests. {it:Journal of the Royal Statistical Society, Series B} 85(2): 179-211.

{pstd}
Athey, S., J. Tibshirani, and S. Wager. 2019.
Generalized Random Forests. {it:Annals of Statistics} 47(2): 1148-1178.
