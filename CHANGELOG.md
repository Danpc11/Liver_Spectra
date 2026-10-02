# Changelog

## v1.0.0 — 2026-10-02
First public release, accompanying the submitted manuscripts.

Pipeline (scripts 01–16) reproduces every Supplementary Table (S1–S16) and both figure sets
(`results/figures/genome_research`, `results/figures/jhep`) from the raw inputs.

Verified in this release:
* clean end-to-end re-run on one core, < 3 GB RAM, ~35 min;
* automated cross-check of 60 numbers quoted in the manuscripts against the regenerated tables
  (`results/tables/audit_number_checks.csv`): 58 matched, 2 corrected (lineage-intrinsic targets 154
  not 152; distance-stratified TAD test P = 0.26 not ≈ 0.5). See `AUDIT.md`.

Known limitations and planned work are listed in `AUDIT.md` section 3.
