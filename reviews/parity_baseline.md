# Parity Baseline (Gaps 4/5/6)

## Upstream source of truth
- Repository: `grf-labs/grf`
- Pinned release tag: `v2.5.0`
- Release line used for this remediation cycle: `2.5.0` (latest stable as of 2026-03-02)
- API extraction artifact: `reviews/r_api_manifest.json`
- Generator: `tools/extract_r_api_manifest.R`

## Vendor status in this repository
- Local vendored core path: `vendor/grf/core`
- This vendored copy does not carry git metadata in-tree, so parity scope is anchored to the pinned upstream API manifest rather than local commit introspection.
- Runtime build remains against the existing vendored core for this cycle; parity gating is enforced by manifest checks.

## Machine-checked scope outcomes
From `reviews/r_api_manifest.json` generated against `v2.5.0`:
- Gap 4 (`orthog.boosting`): **not in current upstream API**.
- Gap 5 (`honesty.prune.method` enum): **not in current upstream API**.
- Gap 6 nuisance surface for `causal_survival_forest`: upstream exposes `W.hat`; does not expose `Y.hat`, `S.hat`, or `C.hat` inputs.

## Reproduction
```bash
Rscript tools/extract_r_api_manifest.R --tag v2.5.0 --out reviews/r_api_manifest.json
Rscript tools/check_parity_scope.R reviews/r_api_manifest.json
```
