#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
manifest_path <- if (length(args) >= 1) args[[1]] else "reviews/r_api_manifest.json"

if (!file.exists(manifest_path)) {
  stop(sprintf("Manifest not found: %s", manifest_path))
}

m <- read_json(manifest_path, simplifyVector = TRUE)

expect_false <- function(flag, label) {
  if (isTRUE(flag)) {
    stop(sprintf("Expected FALSE for %s, got TRUE", label))
  }
}

expect_true <- function(flag, label) {
  if (!isTRUE(flag)) {
    stop(sprintf("Expected TRUE for %s, got FALSE", label))
  }
}

expect_false(m$options[["orthog.boosting"]]$present, "options.orthog.boosting.present")
expect_false(m$options[["honesty.prune.method"]]$present, "options.honesty.prune.method.present")
expect_true(m$options[["causal_survival.W.hat"]]$present, "options.causal_survival.W.hat.present")
expect_false(m$options[["causal_survival.Y.hat"]]$present, "options.causal_survival.Y.hat.present")
expect_false(m$options[["causal_survival.S.hat"]]$present, "options.causal_survival.S.hat.present")
expect_false(m$options[["causal_survival.C.hat"]]$present, "options.causal_survival.C.hat.present")

cat("Parity scope checks passed.\n")
