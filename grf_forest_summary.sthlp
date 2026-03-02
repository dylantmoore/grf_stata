{smcl}
{* *! version 0.2.0}{...}
{title:grf_forest_summary}

{pstd}
{cmd:grf_forest_summary} prints metadata for the most recent fitted GRF model
and stores key fields in {cmd:r()} ({cmd:cmd}, {cmd:forest_type}, {cmd:N},
{cmd:n_trees}, {cmd:seed} when available).

{pstd}
Syntax:
{cmd:grf_forest_summary} [{cmd:,} {opt all}]

{pstd}
With {cmd:all}, the command also lists all stored {cmd:e()} scalar and macro
names and returns them as {cmd:r(e_scalars)} and {cmd:r(e_macros)}.
