{smcl}
{* *! version 0.3.0}{...}
{viewerjumpto "Syntax" "grf_causal_survival_forest##syntax"}{...}
{viewerjumpto "Description" "grf_causal_survival_forest##description"}{...}
{viewerjumpto "Nuisance modes" "grf_causal_survival_forest##modes"}{...}
{viewerjumpto "Options" "grf_causal_survival_forest##options"}{...}
{viewerjumpto "Examples" "grf_causal_survival_forest##examples"}{...}
{viewerjumpto "Stored results" "grf_causal_survival_forest##results"}{...}

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

{synoptset 34 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Main}
{synopt:{opt gen:erate(newvar)}}name for CATE predictions{p_end}
{synopt:{opt replace}}overwrite existing generated variables{p_end}

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

{syntab:Nuisance inputs (full-input mode)}
{synopt:{opt whatinput(varname)}}user-supplied propensity nuisance {it:W.hat}{p_end}
{synopt:{opt yhatinput(varname)}}user-supplied outcome nuisance {it:Y.hat}{p_end}
{synopt:{opt shatinput(varname)}}user-supplied survival nuisance {it:S.hat}{p_end}
{synopt:{opt chatinput(varname)}}user-supplied censoring nuisance {it:C.hat}{p_end}
{synopt:{opt numer(varname)}}expert moment override numerator{p_end}
{synopt:{opt denom(varname)}}expert moment override denominator{p_end}

{syntab:Nuisance outputs}
{synopt:{opt whatgenerate(newvar)}}store wrapper nuisance {it:W.hat}{p_end}
{synopt:{opt yhatgenerate(newvar)}}store wrapper nuisance {it:Y.hat}{p_end}
{synopt:{opt shatgenerate(newvar)}}store wrapper nuisance {it:S.hat}{p_end}
{synopt:{opt chatgenerate(newvar)}}store wrapper nuisance {it:C.hat}{p_end}

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
survival outcomes using a causal-survival objective and nuisance moments.

{marker modes}{...}
{title:Nuisance modes}

{phang}
{bf:auto} (default): nuisance surfaces are estimated internally and mapped to
nuisance moments used by the forest objective.

{phang}
{bf:full_input}: provide all of {cmd:whatinput()}, {cmd:yhatinput()},
{cmd:shatinput()}, and {cmd:chatinput()}.

{phang}
{bf:moment_input}: provide {cmd:numer()} and {cmd:denom()} directly.

{pstd}
{cmd:numer()/denom()} is mutually exclusive with nuisance-input options.
Partial nuisance input (for example only {cmd:whatinput()}) is rejected.

{marker options}{...}
{title:Options}

{phang}
{opt numer(varname)} and {opt denom(varname)} specify expert moment overrides.
These must be supplied together.

{phang}
{opt whatinput(varname)} {opt yhatinput(varname)} {opt shatinput(varname)}
and {opt chatinput(varname)} define full-input nuisance mode and must all be
provided together.

{phang}
In the current wrapper implementation, {cmd:S.hat} contributions are used by
the nuisance moment construction for {cmd:target(2)} (survival probability);
for {cmd:target(1)} (RMST), the moment uses event-term weighting.

{phang}
{opt whatgenerate()}, {opt yhatgenerate()}, {opt shatgenerate()}, and
{opt chatgenerate()} write nuisance surfaces used by the wrapper to user
variables.

{marker examples}{...}
{title:Examples}

{pstd}Default nuisance estimation{p_end}
{phang2}{cmd:. grf_causal_survival_forest time status w x1 x2 x3, gen(cs_tau)}{p_end}

{pstd}Full nuisance-input mode{p_end}
{phang2}{cmd:. grf_causal_survival_forest time status w x1 x2 x3, gen(cs_tau2) horizon(4) whatinput(wh) yhatinput(yh) shatinput(sh) chatinput(ch)}{p_end}

{pstd}Moment-input override mode{p_end}
{phang2}{cmd:. grf_causal_survival_forest time status w x1 x2 x3, gen(cs_tau3) numer(num) denom(den)}{p_end}

{marker results}{...}
{title:Stored results}

{pstd}
{cmd:grf_causal_survival_forest} stores the following in {cmd:e()}:

{synoptset 30 tabbed}{...}
{p2col 5 30 34 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}observations used{p_end}
{synopt:{cmd:e(n_events)}}event count{p_end}
{synopt:{cmd:e(n_censored)}}censored count{p_end}
{synopt:{cmd:e(n_trees)}}number of trees{p_end}
{synopt:{cmd:e(horizon)}}horizon used{p_end}
{synopt:{cmd:e(target)}}target type (1=RMST, 2=survival probability){p_end}
{synopt:{cmd:e(ate)}}mean CATE prediction{p_end}
{synopt:{cmd:e(ate_se)}}SE of mean CATE prediction{p_end}

{p2col 5 30 34 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:grf_causal_survival_forest}{p_end}
{synopt:{cmd:e(forest_type)}}{cmd:causal_survival}{p_end}
{synopt:{cmd:e(predict_var)}}CATE output variable{p_end}
{synopt:{cmd:e(variance_var)}}variance output variable (if requested){p_end}
{synopt:{cmd:e(nuisance_mode)}}{cmd:auto}, {cmd:full_input}, or {cmd:moment_input}{p_end}
{synopt:{cmd:e(what_var)}}canonical nuisance variable ({cmd:_grf_cs_what}){p_end}
{synopt:{cmd:e(yhat_var)}}canonical nuisance variable ({cmd:_grf_cs_yhat}){p_end}
{synopt:{cmd:e(shat_var)}}canonical nuisance variable ({cmd:_grf_cs_shat}){p_end}
{synopt:{cmd:e(chat_var)}}canonical nuisance variable ({cmd:_grf_cs_chat}){p_end}
{synopt:{cmd:e(numer_var)}}canonical nuisance variable ({cmd:_grf_cs_numer}){p_end}
{synopt:{cmd:e(denom_var)}}canonical nuisance variable ({cmd:_grf_cs_denom}){p_end}
{synopt:{cmd:e(what_generate)}}user-generated nuisance output (if requested){p_end}
{synopt:{cmd:e(yhat_generate)}}user-generated nuisance output (if requested){p_end}
{synopt:{cmd:e(shat_generate)}}user-generated nuisance output (if requested){p_end}
{synopt:{cmd:e(chat_generate)}}user-generated nuisance output (if requested){p_end}
