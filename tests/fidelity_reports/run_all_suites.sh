#!/bin/bash
# Master runner for all 15 R-vs-Stata fidelity test suites
#
# Usage:
#   cd /tmp/grf_stata/tests/fidelity_reports
#   bash run_all_suites.sh           # run all suites
#   bash run_all_suites.sh 01 07     # run only suites 01 and 07
#
# Prerequisites:
#   - R 4.5+ with grf 2.5.0 installed
#   - Stata MP (StataNow 19.5+) at $STATA or on PATH
#   - grf_stata package at /tmp/grf_stata with compiled plugin
#
# Each suite follows a 3-step pipeline:
#   1. R generates data + R reference predictions → CSV files
#   2. Stata loads CSVs, runs matching commands, exports Stata predictions
#   3. R computes Pearson correlations and summarizes results

set -euo pipefail

BASEDIR="$(cd "$(dirname "$0")" && pwd)"
STATA="${STATA:-/Applications/StataNow/StataMP.app/Contents/MacOS/stata-mp}"

# Suite definitions: directory, R script(s), Stata script, comparison script
# Format: suite_id|directory|r_step|stata_step|compare_step
SUITES=(
  "01|01_regression|MULTI_R|run_all_stata.do|compare_all.R"
  "02|02_causal_fit|run_r_tests.R|run_stata_tests.do|INLINE"
  "03|03_ate|run_r.R|run_stata.do|compare.R"
  "04|04_instrumental|run_all_r.R|run_all_stata.do|INLINE"
  "05|05_quantile|MULTI_R|MULTI_DO|analyze_results.R"
  "06|06_probability|MULTI_R|MULTI_DO|compare_all.R"
  "07|07_survival|run_all_r.R|run_all_stata.do|analyze_results.R"
  "08|08_causal_survival|run_r_tests.R|run_stata_tests.do|compute_correlations.R"
  "09|09_multi_arm|run_r_tests.R|run_stata_tests.do|INLINE"
  "10|10_boosted|MULTI_R|run_all_stata.do|compute_correlations.R"
  "11|11_ll_regression|run_all_r.R|run_all_stata.do|compute_correlations.R"
  "12|12_lm_forest|01_gen_data_and_r_tests.R|02_stata_tests.do|03_compare_and_report.R"
  "13|13_multi_regression|run_all_r.R|run_all_stata.do|compute_correlations.R"
  "14|14_blp_calibration|r_reference.R|stata_blp.do|INLINE"
  "15|15_vi_rate_tune|test15_r.R|test15_stata.do|compare15.R"
)

# Parse arguments: if specific suite IDs given, run only those
REQUESTED=("$@")

run_suite() {
  local IFS='|'
  read -r id dir r_step stata_step compare_step <<< "$1"

  # Check if this suite was requested (if any were specified)
  if [ ${#REQUESTED[@]} -gt 0 ]; then
    local found=0
    for req in "${REQUESTED[@]}"; do
      if [ "$req" = "$id" ]; then found=1; break; fi
    done
    if [ $found -eq 0 ]; then return 0; fi
  fi

  local suitedir="$BASEDIR/$dir"
  echo ""
  echo "================================================================"
  echo "  Suite $id: $dir"
  echo "================================================================"

  if [ ! -d "$suitedir" ]; then
    echo "  SKIP: directory $suitedir not found"
    return 0
  fi

  cd "$suitedir"

  # Step 1: Run R reference
  echo "  [1/3] Running R reference..."
  if [ "$r_step" = "MULTI_R" ]; then
    for rscript in test*_*.R; do
      [ -f "$rscript" ] && Rscript "$rscript" > /dev/null 2>&1
    done
  else
    if [ -f "$r_step" ]; then
      Rscript "$r_step" > /dev/null 2>&1
    else
      echo "  WARN: R script $r_step not found, skipping"
    fi
  fi
  echo "  [1/3] R reference complete"

  # Step 2: Run Stata
  echo "  [2/3] Running Stata..."
  if [ "$stata_step" = "MULTI_DO" ]; then
    for dofile in test*_*.do; do
      [ -f "$dofile" ] && "$STATA" -b do "$dofile" 2>/dev/null || true
    done
  else
    if [ -f "$stata_step" ]; then
      "$STATA" -b do "$stata_step" 2>/dev/null || true
    else
      echo "  WARN: Stata script $stata_step not found, skipping"
    fi
  fi
  echo "  [2/3] Stata complete"

  # Step 3: Compare
  echo "  [3/3] Computing correlations..."
  if [ "$compare_step" = "INLINE" ]; then
    echo "  (comparison is inline in the Stata do-file — see results CSV)"
  elif [ -f "$compare_step" ]; then
    Rscript "$compare_step" 2>/dev/null || true
  else
    echo "  WARN: comparison script $compare_step not found"
  fi
  echo "  [3/3] Comparison complete"

  echo "  Suite $id DONE"
}

echo "grf_stata Fidelity Test Runner"
echo "=============================="
echo "Base directory: $BASEDIR"
echo "Stata binary:   $STATA"
echo "Suites:         ${#SUITES[@]}"
if [ ${#REQUESTED[@]} -gt 0 ]; then
  echo "Running only:   ${REQUESTED[*]}"
fi
echo ""

STARTTIME=$(date +%s)

for suite in "${SUITES[@]}"; do
  run_suite "$suite"
done

ENDTIME=$(date +%s)
ELAPSED=$((ENDTIME - STARTTIME))

echo ""
echo "=============================="
echo "All requested suites complete."
echo "Elapsed: ${ELAPSED}s"
echo "=============================="
