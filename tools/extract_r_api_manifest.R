#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
tag <- "v2.5.0"
out <- "reviews/r_api_manifest.json"

if (length(args) > 0) {
  i <- 1
  while (i <= length(args)) {
    if (args[[i]] == "--tag") {
      i <- i + 1
      tag <- args[[i]]
    } else if (args[[i]] == "--out") {
      i <- i + 1
      out <- args[[i]]
    } else {
      stop(sprintf("Unknown argument: %s", args[[i]]))
    }
    i <- i + 1
  }
}

base <- sprintf("https://raw.githubusercontent.com/grf-labs/grf/%s/r-package/grf", tag)

fetch_lines <- function(path) {
  url <- sprintf("%s/%s", base, path)
  txt <- tryCatch(
    readLines(url, warn = FALSE),
    error = function(e) stop(sprintf("Failed to fetch %s (%s)", url, e$message))
  )
  list(path = path, url = url, text = txt)
}

files <- list(
  causal_forest = fetch_lines("man/causal_forest.Rd"),
  causal_survival_forest = fetch_lines("man/causal_survival_forest.Rd"),
  survival_forest = fetch_lines("man/survival_forest.Rd"),
  grf_options = fetch_lines("man/grf_options.Rd"),
  plot_rate = fetch_lines("man/plot.rank_average_treatment_effect.Rd"),
  get_scores_r = fetch_lines("R/get_scores.R")
)

contains_pattern <- function(file_obj, pattern) {
  any(grepl(pattern, file_obj$text, fixed = TRUE))
}

option_manifest <- list(
  `orthog.boosting` = list(
    present = contains_pattern(files$causal_forest, "orthog.boosting") ||
      contains_pattern(files$causal_survival_forest, "orthog.boosting"),
    checked_in = c(files$causal_forest$path, files$causal_survival_forest$path)
  ),
  `honesty.prune.method` = list(
    present = contains_pattern(files$causal_forest, "honesty.prune.method") ||
      contains_pattern(files$causal_survival_forest, "honesty.prune.method") ||
      contains_pattern(files$survival_forest, "honesty.prune.method"),
    checked_in = c(files$causal_forest$path, files$causal_survival_forest$path, files$survival_forest$path)
  ),
  `compute.oob.predictions` = list(
    present = contains_pattern(files$causal_forest, "compute.oob.predictions") ||
      contains_pattern(files$causal_survival_forest, "compute.oob.predictions"),
    checked_in = c(files$causal_forest$path, files$causal_survival_forest$path)
  ),
  `causal_survival.W.hat` = list(
    present = contains_pattern(files$causal_survival_forest, "W.hat"),
    checked_in = c(files$causal_survival_forest$path)
  ),
  `causal_survival.Y.hat` = list(
    present = contains_pattern(files$causal_survival_forest, "Y.hat"),
    checked_in = c(files$causal_survival_forest$path)
  ),
  `causal_survival.S.hat` = list(
    present = contains_pattern(files$causal_survival_forest, "S.hat"),
    checked_in = c(files$causal_survival_forest$path)
  ),
  `causal_survival.C.hat` = list(
    present = contains_pattern(files$causal_survival_forest, "C.hat"),
    checked_in = c(files$causal_survival_forest$path)
  ),
  `causal_survival.failure.times` = list(
    present = contains_pattern(files$causal_survival_forest, "failure.times"),
    checked_in = c(files$causal_survival_forest$path)
  )
)

function_manifest <- list(
  `get_scores.causal_survival_forest` = list(
    present = contains_pattern(files$get_scores_r, "get_scores.causal_survival_forest <- function"),
    checked_in = c(files$get_scores_r$path)
  ),
  `plot.rank_average_treatment_effect` = list(
    present = contains_pattern(files$plot_rate, "plot.rank_average_treatment_effect"),
    checked_in = c(files$plot_rate$path)
  ),
  `grf_options` = list(
    present = contains_pattern(files$grf_options, "\\name{grf_options}"),
    checked_in = c(files$grf_options$path)
  )
)

gap_scope <- list(
  gap4_orthog_boosting = list(
    option = "orthog.boosting",
    in_current_upstream_api = option_manifest[["orthog.boosting"]]$present,
    resolution = if (option_manifest[["orthog.boosting"]]$present) "implement" else "non_applicable_in_current_upstream"
  ),
  gap5_honesty_prune_method = list(
    option = "honesty.prune.method",
    in_current_upstream_api = option_manifest[["honesty.prune.method"]]$present,
    resolution = if (option_manifest[["honesty.prune.method"]]$present) "implement" else "non_applicable_in_current_upstream"
  ),
  gap6_causal_survival_nuisance_inputs = list(
    target_function = "causal_survival_forest",
    W_hat_supported = option_manifest[["causal_survival.W.hat"]]$present,
    Y_hat_supported = option_manifest[["causal_survival.Y.hat"]]$present,
    S_hat_supported = option_manifest[["causal_survival.S.hat"]]$present,
    C_hat_supported = option_manifest[["causal_survival.C.hat"]]$present
  )
)

manifest <- list(
  schema_version = 1,
  generated_at_utc = format(as.POSIXct(Sys.time(), tz = "UTC"), "%Y-%m-%dT%H:%M:%SZ"),
  source = list(
    repository = "grf-labs/grf",
    release_tag = tag,
    extraction_mode = "raw_github_docs_and_source",
    files = unname(lapply(files, function(x) list(path = x$path, url = x$url)))
  ),
  options = option_manifest,
  functions = function_manifest,
  gap_scope = gap_scope
)

dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
write_json(manifest, out, auto_unbox = TRUE, pretty = TRUE)

cat(sprintf("Wrote manifest: %s\n", out))
cat(sprintf("  gap4 orthog.boosting present: %s\n", gap_scope$gap4_orthog_boosting$in_current_upstream_api))
cat(sprintf("  gap5 honesty.prune.method present: %s\n", gap_scope$gap5_honesty_prune_method$in_current_upstream_api))
cat(sprintf("  gap6 W.hat/Y.hat/S.hat/C.hat: %s/%s/%s/%s\n",
            gap_scope$gap6_causal_survival_nuisance_inputs$W_hat_supported,
            gap_scope$gap6_causal_survival_nuisance_inputs$Y_hat_supported,
            gap_scope$gap6_causal_survival_nuisance_inputs$S_hat_supported,
            gap_scope$gap6_causal_survival_nuisance_inputs$C_hat_supported))
