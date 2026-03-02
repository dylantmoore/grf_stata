{smcl}
{* *! version 0.1.0}{...}
{title:grf_merge_forests}

{pstd}
{cmd:grf_merge_forests} combines two prediction variables into a single merged
score using convex weights.

{pstd}
Syntax:
{cmd:grf_merge_forests} {it:pred1 pred2} [{ifin}], {opt generate(newvar)}
[{opt weights(numlist)} {opt replace}]

{pstd}
If {cmd:weights()} is omitted, equal 0.5/0.5 weights are used.
