# Positional spectra of the human liver transcriptome in MASLD fibrosis

Reproducible pipeline for the manuscript *"The liver reads its genome in neighbourhoods: positional spectra of the
human liver transcriptome reveal an invariant architecture, composition-driven change and cis-coupled gene
neighbourhoods in MASLD fibrosis"*. Every table and figure in the paper is produced by the numbered scripts in
`scripts/`, run in order, from the raw inputs listed below.

```
liver-positional-spectra/
├── README.md               this file (data manifest at the end)
├── requirements.txt        Python 3.12 dependencies (pip install -r requirements.txt)
├── Makefile                `make all` runs scripts 01–10; `make figures` only 09–10
├── data/
│   ├── raw/                INPUTS the user supplies (counts, GEO metadata, positional grid)
│   ├── external/           PUBLIC resources fetched by scripts/00_fetch_external.sh
│   └── curated/            drug_landscape.csv (manual curation of the top candidate targets)
├── scripts/                00–10 (see below) + common.py (paths, constants, marker sets, helpers)
└── results/
    ├── intermediate/       pickles passed between steps (not needed for the paper)
    ├── tables/             Supplementary Tables S1–S11 (CSV) — all GENERATED
    └── figures/            Fig1–Fig7 (PDF + PNG) — all GENERATED
```

## Quick start

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
bash scripts/00_fetch_external.sh          # DoRothEA, Hallmark GMT, liver TADs; then place GSE136103_RAW.tar
make all                                   # ≈ 25 min on one core, < 3 GB RAM
```

Scripts are run from `scripts/` (or via `make`) and locate the repository root automatically. Random seeds are
fixed (`numpy.random.default_rng(0/1)`), so permutation nulls and the shuffled-order control are reproducible to
the printed precision.

## Pipeline

| Step | Script | What it does | Main outputs |
| --- | --- | --- | --- |
| 00 | `00_fetch_external.sh` | downloads DoRothEA (saezlab), Hallmark v7.0 GMT (public mirror), liver TAD partitions (McArthur & Capra); extracts GSE136103 | `data/external/*` |
| 01 | `01_prepare_expression.py` | parses GEO characteristics, defines conditions (Normal/Control/F0–F4), median-of-ratios + log2, cohort correction protecting condition | S1, `expr*.pkl`, `meta.pkl` |
| 02 | `02_spectra.py` | per-sample positional spectra (log, within-cohort z, multitaper), 1/f whitening, permutation null, grid-mask control, universal peaks, condition spectra, stage effects per frequency and band | S2a–e, `W_adj.pkl`, `W_mt.pkl` |
| 03 | `03_specparam.py` | aperiodic offset/exponent and periodic peak count per sample × chromosome; stage models | S3a–d |
| 04 | `04_differential_expression.py` | pyDESeq2 ordinal-stage and F4-vs-Normal; programme enrichment; spatial autocorrelation of the DE statistic; adjacent concordant pairs; contiguous neighbourhoods | S4a–c, S5a–c, `de.pkl` |
| 05 | `05_composition_distance_tads.py` | marker-based composition; composition-adjusted stage effects and ACF; GRCh37 coordinates, orientation, co-expression; liver TAD reconstruction and distance-stratified TAD test; threshold-vs-linear (AIC) and bimodality | S5d–i, S6a–c, S7a–b, `pairs.pkl`, `K_tad.pkl` |
| 06 | `06_single_cell.py` | GSE136103 QC, marker annotation, pseudobulk, within-type cirrhosis-vs-healthy t, within-type ACF and neighbour co-expression | S8a–d |
| 07 | `07_gsea_tf.py` | GSEA pre-ranked (Hallmark) on ordinal and composition-adjusted statistics; DoRothEA A–C TF activity vs stage; shared TFs in coupled pairs; TF × Hallmark overlap | S10a–h |
| 08 | `08_targets_fingerprint.py` | lineage-intrinsic target prioritisation (bulk × composition × single cell, + curated drug landscape); spectral fingerprint with real vs shuffled gene order (leave-one-cohort-out) | S9a–b, S11 |
| 09 | `09_make_figures.py` | Figures 1–5 (Nature style, 183 mm, Okabe–Ito / viridis) | `Fig1`–`Fig5` |
| 10 | `10_make_circos.py` | Figures 6–7 (circos genome map; TF × Hallmark chord) | `Fig6`, `Fig7`, S5j |
| 12 | `12_eqtl_shared_variants.py` | liver eQTL credible sets (eQTL Catalogue QTD000266): shared causal variants in coupled vs uncoupled neighbours | S13a–c |
| 15 | `15_gwas_loci_neighbourhoods.py` | MASLD/cirrhosis GWAS loci (curated + GWAS Catalog with non-liver controls) in coupled neighbourhoods | S16, S16b–d |
| 16 | `16_make_figures_jhep.py` | Figures 1–6 and Table 1 for the Journal of Hepatology version (clinical framing) | `results/figures/jhep/` |
| 14 | `14_eqtl_gtex_signif_pairs.py` | GTEx v8 liver significant eQTL pairs: shared (same-direction) variants in coupled vs uncoupled neighbours, MH stratified by distance | S13d–g |
| 13 | `13_gtex_tissues.py` | spectra in 11 GTEx tissues; tissue specificity of the architecture; replication of biopsy peaks | S14a–f, Extended Data Fig. 2 |
| 11 | `11_mouse_validation.py` | own mouse fatty-liver data (CTL vs UTP, n = 3 + 3): DESeq2, invariant spectrum, spatial coupling, distance decay, syntenic human pairs | S12a–e, Extended Data Fig. 1 |

## Data manifest

### A. Inputs to upload (not generated; `data/raw/`)

| File | Content | Size | Source |
| --- | --- | --- | --- |
| `counts/counts_GSE130970.tsv` | raw gene counts, 26,808 Ensembl genes × 78 samples | 12 MB | GEO GSE130970 (Hoang 2019) |
| `counts/counts_GSE135251.tsv` | raw gene counts × 216 samples | 15 MB | GEO GSE135251 (Govaere 2020) |
| `counts/counts_GSE162694.tsv` | raw gene counts × 143 samples | 10 MB | GEO GSE162694 (Pantano 2021) |
| `metadata/metadata_cruda_GSE*.tsv` (3) | GEO sample characteristics (fibrosis stage, NAS, sex, age) | <100 kB | GEO series matrices |
| `own/mcounts.tsv` | own mouse liver RNA-seq counts (CTL1–3, UTP1–3; ARC-UTP columns present but not used) | 3 MB | this study (to be deposited in GEO) |
| `grid/<GSE>_rejilla_genes.tsv` (5) | positional grid: chr, grid_index, gene_id, gene_name (same gene → same index in every file; union is used) | 2–3 MB | this project |

### B. Public external resources (fetched by `00_fetch_external.sh`; `data/external/`)

| File | Content | Source |
| --- | --- | --- |
| `dorothea_hs.rda` | DoRothEA human regulons (tf, confidence, target, mor) | github.com/saezlab/dorothea |
| `hallmark.gmt` | MSigDB Hallmark v7.0 symbols | public mirror (replace with official MSigDB download) |
| `TAD-stability-heritability-master/data/20binsTADlandscape/Liver_leung2015/` | liver TAD partitions (hg19), 20 bins per domain | github.com/emcarthur/TAD-stability-heritability |
| `GSE136103/*_{matrix.mtx,genes.tsv,barcodes.tsv}.gz` | Ramachandran 2019 human liver scRNA-seq (20 liver samples used) | GEO GSE136103 (`GSE136103_RAW.tar`, 436 MB, manual download) |
| `QTD000266.credible_sets.tsv.gz` | GTEx liver eQTL fine-mapped credible sets | eQTL Catalogue FTP `susie/QTS000015/QTD000266/` |
| `gwas-catalog-download-associations-v1.0-full.tsv` | GWAS Catalog full associations (P < 5e-8 filtered in script 15) | GWAS Catalog downloads |
| `Liver.v8.signif_variant_gene_pairs.txt.gz`, `Liver.v8.egenes.txt.gz` | GTEx v8 liver single-tissue cis-eQTL (significant pairs, eGenes) | GTEx Portal, QTL downloads (GTEx_Analysis_v8_eQTL.tar) |
| `gtex/gene_reads_*_<tissue>.gct.gz` (11) | GTEx gene read counts per tissue (v11; whole blood v10) | GTEx Portal, bulk tissue expression |
| GRCh37 Ensembl 100 gene table | coordinates and strand | bundled in the `pyannotables` package |

### C. Generated tables (`results/tables/`; Supplementary Tables of the manuscript)

| Table | File | Content | Script |
| --- | --- | --- | --- |
| S1 | `Table_S1_samples.csv` | 437 samples: cohort, condition, ordinal stage, sex, age, NAS | 01 |
| S2a | `Table_S2a_consensus_spectrum_universal_peaks.csv` | 10,007 frequencies: period, mean power/null, fraction of samples >3× (overall and per cohort), grid-mask power, universal flag | 02 |
| S2b | `Table_S2b_condition_spectra_{log,z,multitaper}.csv` | consensus spectrum per condition (mean log w, s.e., t vs null, q) | 02 |
| S2c | `Table_S2c_stage_effect_per_frequency_{…}.csv` | ordinal-stage t, p, q per frequency | 02 |
| S2d | `Table_S2d_stage_effect_per_band_{…}.csv` | stage t per period band | 02 |
| S2e | `Table_S2e_band_profile_by_condition_{…}.csv` | mean log w per band and condition | 02 |
| S3a–b | `Table_S3a_specparam_per_sample_chromosome.csv`, `Table_S3b_specparam_peaks_per_sample.csv` | aperiodic offset/exponent, residual s.d., peak count; individual peaks | 03 |
| S3c–d | `Table_S3c_specparam_vs_stage.csv`, `Table_S3d_exponent_per_chromosome_vs_stage.csv` | stage models of the parameters; per-chromosome exponent | 03 |
| S4a–b | `Table_S4a_DESeq2_F4_vs_Normal.csv`, `Table_S4b_DESeq2_ordinal_stage.csv` | DESeq2 results | 04 |
| S4c | `Table_S4c_programme_enrichment_DE.csv` | Fisher enrichment of curated programmes | 04 |
| S5a–c | `Table_S5a_spatial_autocorrelation_DE.csv`, `Table_S5b_adjacent_DE_pairs.csv`, `Table_S5c_contiguous_DE_neighbourhoods.csv` | ACF with permutation null; adjacent pairs vs null; 422 runs ≥3 genes | 04 |
| S5d | `Table_S5d_acf_with_without_composition.csv` | ACF of per-gene stage effect, raw and composition-adjusted | 05 |
| S5e | `Table_S5e_liver_TADs_reconstructed_hg19.csv` | 1,304 reconstructed liver TADs | 05 |
| S5f–h | `Table_S5f_adjacent_pairs_distance_orientation_TAD.csv`, `Table_S5g_concordance_by_distance.csv`, `Table_S5h_concordance_by_orientation.csv` | all adjacent pairs with distance, orientation, family, TAD, statistics, co-expression; summaries | 05 |
| S5i | `Table_S5i_TAD_test_stratified.csv` | same- vs different-TAD difference at equal distance (permutation) | 05 |
| S5j | `Table_S5j_neighbourhood_pathway_links.csv` | neighbourhood → Hallmark leading-edge links (Fig. 6) | 10 |
| S6a–c | `Table_S6a_marker_counts.csv`, `Table_S6b_composition_scores_per_sample.csv`, `Table_S6c_stage_effect_per_gene_with_without_composition.csv` | composition scores and per-gene stage t with/without adjustment | 05 |
| S7a–b | `Table_S7a_threshold_vs_linear.csv`, `Table_S7b_within_stage_bimodality.csv` | ΔAIC of step models; GMM bimodality | 05 |
| S8a–d | `Table_S8a_single_cell_annotation.csv`, `Table_S8b_within_type_neighbour_coexpression.csv`, `Table_S8c_within_type_spatial_autocorrelation.csv`, `Table_S8d_within_type_t_cirrhosis_vs_healthy.csv` | single-cell annotation, co-expression, ACF, per-type t | 06 |
| S9a–b | `Table_S9a_candidate_targets_lineage_intrinsic.csv`, `Table_S9b_integration_bulk_sc_per_gene.csv` | prioritised targets with drug landscape; full integration table | 08 |
| S10a–h | `Table_S10a_GSEA_hallmark_ordinal.csv`, `Table_S10b_GSEA_hallmark_composition_adjusted.csv`, `Table_S10c_dorothea_ABC_regulons_used.csv`, `Table_S10d_DoRothEA_TF_activity_vs_stage.csv`, `Table_S10e_TF_activity_per_sample.csv`, `Table_S10f_shared_TF_coupled_pairs.csv`, `Table_S10g_shared_TFs_by_direction.csv`, `Table_S10h_TF_x_Hallmark_overlap.csv` | GSEA, TF activity, shared regulators, TF × pathway matrix | 07 |
| S11 | `Table_S11_spectral_fingerprint_real_vs_shuffled.csv` | leave-one-cohort-out classification, real vs shuffled gene order | 08 |
| S12a–e | `Table_S12a_mouse_DESeq2_UTP_vs_CTL.csv`, `Table_S12b_mouse_spectrum_universal_peaks.csv`, `Table_S12c_mouse_spatial_autocorrelation.csv`, `Table_S12d_mouse_concordance_by_distance.csv`, `Table_S12e_mouse_syntenic_pairs.csv` | mouse validation: DE, spectrum, ACF, distance decay, syntenic pairs | 11 |

### D. Generated figures (`results/figures/`, PDF vector + PNG 400 dpi)

Fig1 invariant architecture · Fig2 stage effects on the spectrum · Fig3 spatial coupling and TADs ·
Fig4 composition and thresholds · Fig5 single-cell coupling and targets · Fig6 genome circos · Fig7 TF × Hallmark chord ·
Extended Data Fig1 mouse validation.

## Audit

`AUDIT.md` reports a clean re-run of the pipeline, an automated cross-check of 60 numbers quoted in the
manuscripts against the regenerated tables (`results/tables/audit_number_checks.csv`; 58 matched, 2 corrected),
a component-by-component methodological review, and the open items before submission.

## Notes on reproducibility

* Key numbers checked against the manuscript on a clean run: 13,816 genes; 186 universal peaks; stage-associated
  frequencies 6,042 (log) / 0 (z) / 6,345 (multitaper); DE padj<0.05 8,713 (F4 vs Normal) / 7,520 (ordinal);
  adjacent concordant pairs 954/957 vs 692/654; 422 neighbourhoods; composition explains 79.8% of stage-effect
  variance; offset stage t −8.1 → +1.3 after composition; 154/70 TFs before/after composition; fingerprint AUC 0.92
  (real) vs 0.91 (shuffled).
* `pyannotables` ships Ensembl 100 GRCh37; TADs are hg19 — both coordinate systems are GRCh37-consistent.
* The Hallmark GMT is v7.0 from a public mirror; for submission re-run step 07 with the official MSigDB release
  and, optionally, CollecTRI in place of DoRothEA (same column layout: tf, target, mor).
* Single-cell annotation is marker-based (no clustering) and leaves ~50% of cells unassigned by design; see
  Methods.
