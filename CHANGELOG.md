# Changelog

## v1.0.0 additions (pre-release)
* `15b_mechanical_switch.py`: gene-level, pathway, driver and positional analysis of the F3→F4 switch (Supplementary Table S17; Fig. 3D–H of the JHEP figure set).
* Hepatocytes admitted to the scRNA-seq analyses (≥40 cells per donor, flagged low coverage); bulk-defined hepatocyte class and receptors of hepatocyte-directed drugs (S9c–d, S17i) in Fig. 6C,F,G.
* `18_snrnaseq_gse202379.py`: validation in snRNA-seq of 47 SAF-staged donors (S18a–g; Supplementary Fig. S8).
* GWAS figure moved to Supplementary; circos and two-layer/target panels split into Figs 5 and 6.

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
