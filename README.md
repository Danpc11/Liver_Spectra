# Liver Spectra

### Cellular composition and mechanotransduction signatures across fibrosis stages in MASLD

[![Python](https://img.shields.io/badge/Python-3.12-3776AB?logo=python&logoColor=white)](https://www.python.org/)
[![DOI](https://img.shields.io/badge/DOI-pending-lightgrey)](#citation)

---

## What this is

This repository reproduces every table and figure of a study of how the human liver transcriptome changes from
normal histology to cirrhosis (F4) in metabolic dysfunction-associated steatotic liver disease (MASLD): which part of
that change is associated with cell composition, which part remains after adjusting for composition, and whether a
mechanotransduction signature appears at the transition to cirrhosis.

The study separates **two components** of the fibrosis-associated transcriptome:

| Component | What it is | How it behaves |
| --- | --- | --- |
| **Composition-associated** | changes tracking marker-based scores for hepatocytes, hepatic stellate cells, immune cells and ductular epithelium | largely linear with stage; adjusting for it reduces the variance of gene-level stage statistics by about **80%** |
| **Residual (within-lineage)** | changes that remain after adjustment and are supported within cell populations in single-cell and single-nucleus data | includes a YAP/TAZ mechanotransduction signature with an **additional increase at F4** |

More than half of F3 biopsies have an F4-like composition but lower YAP/TAZ target scores than F4 biopsies, so the F4
signature carries information beyond estimated composition. These are cross-sectional data: they support a model of a
mechanical transition but do not establish causality or irreversibility.

---

## Key findings

1. **A stage-stable positional architecture.** Ordering genes by chromosomal position and computing positional
   spectra in 437 biopsies reveals 186 spatial frequencies above the 1/f background in ≥90% of biopsies, most enriched
   in liver among 11 GTEx tissues and unchanged from normal liver to F4 (Fig. 1).
2. **Composition changes with thresholds.** Composition scores change largely linearly with stage, with additional
   steps at fibrosis onset and, for cholangiocytes, at F3–F4 (Fig. 2).
3. **An F4-associated mechanotransduction signature.** After composition adjustment, 291 genes increase and 547
   decrease at F4 beyond their linear trends, including YAP/TAZ targets (*CCN1*, *CCN2*, *AMOTL2*, *THBS1*), with
   SMAD, NF-κB, STAT3 and HIF1A activity (Fig. 3).
4. **Coordinated neighbouring genes within lineages.** Stage effects are autocorrelated between neighbouring genes
   (lag-1 r = 0.17, z = 20), independent of TADs, and coupled within cell populations, including hepatocyte nuclei
   (z = 17). Stage-coupled pairs share same-direction liver eQTLs more often (odds ratio 1.76) (Fig. 4).
5. **Candidate genes with different trajectories.** Of 154 within-lineage candidates, 86 rise with linear trends
   (*THBS2*, *TREM2*, *IL32*, *LGALS3*) and 17 increase at F4 (*CCN2*/CTGF, *CCN1*, *TIMP1*, *VWF*). Targets of
   hepatocyte-directed drugs decline mainly with hepatocyte composition; THRB also falls within hepatocyte nuclei
   (Figs. 5–6).

---

## Included datasets

| Dataset | Content | Role here | Reference |
| --- | --- | --- | --- |
| [GSE130970](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE130970) | bulk RNA-seq, 78 liver biopsies | discovery cohort | Hoang et al., *Sci Rep* 2019 · [10.1038/s41598-019-48746-5](https://doi.org/10.1038/s41598-019-48746-5) |
| [GSE135251](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE135251) | bulk RNA-seq, 216 liver biopsies | discovery cohort | Govaere et al., *Sci Transl Med* 2020 · [10.1126/scitranslmed.aba4448](https://doi.org/10.1126/scitranslmed.aba4448) |
| [GSE162694](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE162694) | bulk RNA-seq, 143 livers incl. normal histology | discovery cohort | Pantano et al., *Sci Rep* 2021 · [10.1038/s41598-021-96966-5](https://doi.org/10.1038/s41598-021-96966-5) |
| [GSE136103](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE136103) | scRNA-seq, 5 healthy and 5 cirrhotic livers | within-lineage validation (non-parenchymal) | Ramachandran et al., *Nature* 2019 · [10.1038/s41586-019-1631-3](https://doi.org/10.1038/s41586-019-1631-3) |
| [GSE202379](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE202379) | snRNA-seq, 59 samples from 47 SAF-staged donors | within-lineage validation incl. hepatocytes | Gribben et al., *Nature* 2024 · [10.1038/s41586-024-07465-2](https://doi.org/10.1038/s41586-024-07465-2) |
| [GTEx v8 / v10–11](https://gtexportal.org) | liver cis-eQTL; bulk expression of 11 tissues | shared eQTLs; tissue specificity | GTEx Consortium, *Science* 2020 · [10.1126/science.aaz1776](https://doi.org/10.1126/science.aaz1776) |
| [eQTL Catalogue QTD000266](https://www.ebi.ac.uk/eqtl/) | fine-mapped liver eQTL credible sets | shared eQTLs (fine-mapped) | Kerimov et al., *PLoS Genet* 2023 · [10.1371/journal.pgen.1010932](https://doi.org/10.1371/journal.pgen.1010932) |
| [GWAS Catalog](https://www.ebi.ac.uk/gwas/) | genome-wide significant associations | MASLD and cirrhosis loci, comparison traits | Sollis et al., *Nucleic Acids Res* 2023 · [10.1093/nar/gkac1010](https://doi.org/10.1093/nar/gkac1010) |
| [Liver TADs](https://github.com/emcarthur/TAD-stability-heritability) | liver topologically associating domains (hg19) | TAD tests | McArthur & Capra, *Am J Hum Genet* 2021 · [10.1016/j.ajhg.2020.12.008](https://doi.org/10.1016/j.ajhg.2020.12.008) |
| [DoRothEA](https://github.com/saezlab/dorothea) · [MSigDB Hallmark](https://www.gsea-msigdb.org) | TF regulons; Hallmark gene sets | TF activity; GSEA | Garcia-Alonso et al., *Genome Res* 2019 · Liberzon et al., *Cell Syst* 2015 |

The three bulk cohorts are shipped in `data/raw/` as raw counts; all other resources are downloaded (see [Reproduce the analysis](#reproduce-the-analysis)), because of their size and their own licences.

---

## Figures

| | | |
| --- | --- | --- |
| [![Fig. 1](results/figures/jhep/Fig1_design_and_positional_framework.png)](results/figures/jhep/Fig1_design_and_positional_framework.pdf) **Fig. 1** Positional architecture | [![Fig. 2](results/figures/jhep/Fig2_linear_replacement_thresholds.png)](results/figures/jhep/Fig2_linear_replacement_thresholds.pdf) **Fig. 2** Composition and thresholds | [![Fig. 3](results/figures/jhep/Fig3_mechanical_switch.png)](results/figures/jhep/Fig3_mechanical_switch.pdf) **Fig. 3** F4 mechanotransduction signature |
| [![Fig. 4](results/figures/jhep/Fig4_neighbourhood_coupling.png)](results/figures/jhep/Fig4_neighbourhood_coupling.pdf) **Fig. 4** Neighbouring genes within lineages | [![Fig. 5](results/figures/jhep/Fig5_genome_circos.png)](results/figures/jhep/Fig5_genome_circos.pdf) **Fig. 5** Genome-wide map | [![Fig. 6](results/figures/jhep/Fig6_two_layers_and_targets.png)](results/figures/jhep/Fig6_two_layers_and_targets.pdf) **Fig. 6** Residual programmes and candidates |

Vector PDFs and 300-dpi TIFFs are in `results/figures/jhep/`; supplementary figures are in
`results/figures/jhep/supplementary/`. A second figure set formatted for *Genome Research* is in
`results/figures/genome_research/`.

---

## Repository structure

```
scripts/
    common.py                    paths, constants, marker sets, shared statistics
    00_fetch_external.sh         download the public resources
    01–08                        expression, positional spectra, DE, composition and TADs, single cell, GSEA/TF, candidates
    12–15b                       eQTL, GTEx tissues, GWAS loci, the F4 signature
    18_snrnaseq_gse202379.py     single-nucleus validation (47 donors)
    09, 10                       Genome Research figure set
    16, 17                       JHEP main figures (Figs. 1–6, Table 1) and supplementary figures
data/raw/                        bulk counts, GEO metadata, gene-order grid (shipped)
data/external/                   downloaded resources (see data/external/README.md)
data/curated/drug_landscape.csv  curated pharmacology of candidate genes, with reviewed evidence categories
results/tables/                  supplementary tables (CSV)
results/figures/                 jhep/ (main and supplementary) and genome_research/ figure sets
CHANGELOG.md · CITATION.cff · LICENSE
```

---

## Reproduce the analysis

```bash
git clone https://github.com/Danpc11/Liver_Spectra.git && cd Liver_Spectra
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt

make external     # DoRothEA, Hallmark GMT, liver TADs; then add the files listed in data/external/README.md
make all          # full pipeline, raw counts to every table and figure (~60 min, one core, < 3 GB RAM)
```

| Target | Runs |
| --- | --- |
| `make external` | downloads the public resources |
| `make analysis` | the analysis steps (tables only) |
| `make figures` | steps 09, 10, 16 and 17 (needs `results/intermediate/` from a previous run) |
| `make all` | everything, in order |
| `make clean` | deletes generated results |

---

## Limitations

- **Cross-sectional data.** The F4 signature is inferred from threshold models and composition adjustment, not
  observed longitudinally; liver stiffness was not measured. The data do not establish causality, bistability or
  irreversibility.
- **Marker-based composition.** Scores are not cell fractions and can change with activation; the ~80% figure is a
  reduction in the variance of test statistics, not a measured fraction of cell replacement.
- **Cohort structure.** Most normal biopsies come from one cohort, so the step between normal/F0 and F1–F4 combines
  disease onset with cohort composition.
- **Few F4 biopsies in snRNA-seq.** GSE202379 has four biopsy-staged F4 donors; F4 estimates lean on five end-stage
  explants and are reported with and without them.
- **eQTL sharing** is enriched among stage-coupled pairs but not among pairs coupled after composition adjustment.
- **GWAS overlap tests are nominal** (15 tests; none survives Bonferroni correction), and gene mapping assigns
  bystander genes to loci.
- **Candidate genes** are a ranking for follow-up, not a validated discovery set with error control.

---

## Publication

**Cellular composition and mechanotransduction signatures across fibrosis stages in MASLD.** Manuscript submitted to
the *Journal of Hepatology*. Until it is published, please cite the archived software release below.

---

## Citation

The software will be archived on Zenodo at the v1.0.0 release (DOI to be added). A machine-readable citation is in
[`CITATION.cff`](CITATION.cff), which GitHub exposes through *Cite this repository*.

---

## Licence

Source code is distributed under the [MIT licence](LICENSE). The third-party datasets analysed here (GEO, GTEx, eQTL
Catalogue, GWAS Catalog, liver TAD partitions) remain under their own terms; values derived from them in
`results/tables/` should be cited together with the original sources listed in [Included datasets](#included-datasets).

---

## Technical reference

### Pipeline

| Step | Script | What it does | Main outputs |
| --- | --- | --- | --- |
| 00 | `00_fetch_external.sh` | downloads DoRothEA, Hallmark v7.0 GMT (public mirror), liver TAD partitions; extracts GSE136103 | `data/external/*` |
| 01 | `01_prepare_expression.py` | parses GEO characteristics, defines conditions (Normal/Control/F0–F4), median-of-ratios + log2, cohort correction protecting condition | S1, `expr*.pkl`, `meta.pkl` |
| 02 | `02_spectra.py` | per-sample positional spectra (periodogram and multitaper), 1/f whitening, grid-mask control, universal peaks, condition spectra, stage effects per frequency and band | S2a–e, `W_adj.pkl`, `W_mt.pkl` |
| 03 | `03_specparam.py` | aperiodic offset and exponent and periodic peak count per sample × chromosome; stage models | S3a–d |
| 04 | `04_differential_expression.py` | pyDESeq2 ordinal stage and F4 vs Normal; spatial autocorrelation of the DE statistic; adjacent concordant pairs; contiguous neighbourhoods | S4a–c, S5a–c, `de.pkl` |
| 05 | `05_composition_distance_tads.py` | marker-based composition; composition-adjusted stage effects; distance, orientation, co-expression; liver TADs, distance-matched TAD test and TAD containment; spectral parameters with composition; threshold-vs-linear (AIC) | S5d–k, S6a–d, S7a–b, `pairs.pkl`, `K_tad.pkl` |
| 06 | `06_single_cell.py` | GSE136103 QC, marker annotation, pseudobulk, within-population statistics, autocorrelation and neighbour co-expression | S8a–e |
| 07 | `07_gsea_tf.py` | pre-ranked GSEA (Hallmark) on raw and composition-adjusted statistics; DoRothEA TF activity; shared TFs in coupled pairs; TF × Hallmark overlap | S10a–h |
| 08 | `08_targets_fingerprint.py` | within-lineage candidates (bulk × composition × single cell, with curated pharmacology); hepatocyte class and drug targets; spectral classification with real vs shuffled gene order | S9a–d, S11 |
| 09 | `09_make_figures.py` | Genome Research figure set, Figs 1–5 | `results/figures/genome_research/` |
| 10 | `10_make_circos.py` | Genome Research Figs 6–7 (circos genome map; TF × Hallmark chord) | `results/figures/genome_research/`, S5j |
| 12 | `12_eqtl_shared_variants.py` | liver eQTL credible sets (eQTL Catalogue QTD000266): shared variants in coupled vs uncoupled neighbours | S13a–c |
| 13 | `13_gtex_tissues.py` | spectra in 11 GTEx tissues; tissue specificity; replication of biopsy peaks | S14a–f |
| 14 | `14_eqtl_gtex_signif_pairs.py` | GTEx v8 liver significant eQTL pairs: shared same-direction variants, Mantel–Haenszel stratified by distance | S13d–g |
| 15 | `15_gwas_loci_neighbourhoods.py` | MASLD and cirrhosis GWAS loci (curated and GWAS Catalog, with comparison traits) in coordinated neighbourhoods | S16, S16b–d |
| 15b | `15b_mechanical_switch.py` | mechanotransduction scores and thresholds, F4-likeness of composition, per-gene F4 step, GSEA, TF activity, positional scale, cell of origin, candidate timing, hepatocyte programme | S7c, S15, S17a–i |
| 18 | `18_snrnaseq_gse202379.py` | snRNA-seq of 47 donors: annotation, donor pseudobulk, concordance with bulk, hepatocyte class, neighbour coupling and co-expression, F4 signature within populations | S18a–h |
| 16 | `16_figures_jhep.py` | Journal of Hepatology Figs 1–6 and Table 1 | `results/figures/jhep/` |
| 17 | `17_supplementary_figures.py` | supplementary figures in order of citation | `results/figures/jhep/supplementary/` |

### Data manifest

#### A. Inputs (`data/raw/`, shipped)

| File | Content | Size | Source |
| --- | --- | --- | --- |
| `counts/counts_GSE130970.tsv` | raw gene counts, 26,808 Ensembl genes × 78 samples | 12 MB | GEO GSE130970 (Hoang 2019) |
| `counts/counts_GSE135251.tsv` | raw gene counts × 216 samples | 15 MB | GEO GSE135251 (Govaere 2020) |
| `counts/counts_GSE162694.tsv` | raw gene counts × 143 samples | 10 MB | GEO GSE162694 (Pantano 2021) |
| `metadata/metadata_GSE*.tsv` (3) | GEO sample characteristics (fibrosis stage, NAS, sex, age) | <100 kB | GEO series matrices |
| `grid/<GSE>_gene_grid.tsv` (5) | positional grid: chr, grid_index, gene_id, gene_name (union of the five files is used) | 2–3 MB | this project |

#### B. External resources (`data/external/`, downloaded)

| File | Content | Source |
| --- | --- | --- |
| `dorothea_hs.rda` | DoRothEA human regulons | github.com/saezlab/dorothea |
| `hallmark.gmt` | MSigDB Hallmark v7.0 symbols | public mirror (replace with the official MSigDB download) |
| `TAD-stability-heritability-master/…/Liver_leung2015/` | liver TAD partitions (hg19) | github.com/emcarthur/TAD-stability-heritability |
| `GSE136103/*_{matrix.mtx,genes.tsv,barcodes.tsv}.gz` | human liver scRNA-seq | GEO GSE136103 (`GSE136103_RAW.tar`, 436 MB) |
| `GSE202379/GSM*_raw_counts_csv.gz` (59) + series matrix | snRNA-seq of 47 MASLD donors | GEO GSE202379 |
| `QTD000266.credible_sets.tsv.gz` | liver eQTL credible sets | eQTL Catalogue FTP `susie/QTS000015/QTD000266/` |
| `gwas-catalog-download-associations-v1.0-full.tsv` | GWAS Catalog full associations (required by step 15) | GWAS Catalog downloads |
| `Liver.v8.signif_variant_gene_pairs.txt.gz`, `Liver.v8.egenes.txt.gz` | GTEx v8 liver cis-eQTL | GTEx Portal (GTEx_Analysis_v8_eQTL.tar) |
| `gtex/gene_reads_*_<tissue>.gct.gz` (11) | GTEx gene read counts per tissue (v11; whole blood v10) | GTEx Portal |
| GRCh37 Ensembl 100 gene table | coordinates and strand | bundled in the `pyannotables` package |

#### C. Generated tables (`results/tables/`)

| Table | Content | Script |
| --- | --- | --- |
| S1 | samples: cohort, condition, ordinal stage, sex, age, NAS | 01 |
| S2a–e | consensus spectrum and universal peaks; condition spectra; stage effects per frequency and band; band profiles | 02 |
| S3a–d | spectral parameterisation per sample and chromosome; stage models | 03 |
| S4a–c | DESeq2 (F4 vs Normal, ordinal stage); programme enrichment | 04 |
| S5a–c | autocorrelation, adjacent pairs and contiguous neighbourhoods of the DE statistic | 04 |
| S5d–k | composition-adjusted autocorrelation; TADs; adjacent pairs; distance; orientation; TAD test; pathway links (S5j, step 10); TAD containment | 05 |
| S6a–d | marker counts; composition scores; per-gene stage effects with and without composition; spectral parameters with composition | 05 |
| S7a–c | threshold vs linear models; within-stage bimodality; F4-likeness (S7c, step 15b) | 05, 15b |
| S8a–e | scRNA-seq annotation, co-expression, autocorrelation, statistics, mean expression by population | 06 |
| S9a–d | within-lineage candidates with pharmacology; bulk × single-cell integration; hepatocyte class; hepatocyte drug targets | 08 |
| S10a–h | GSEA, TF activity, shared regulators, TF × pathway matrix | 07 |
| S11 | spectral classification with real vs shuffled gene order | 08 |
| S13a–g | shared liver eQTL variants (credible sets; GTEx significant pairs) | 12, 14 |
| S14a–f | GTEx tissue spectra and peak sharing | 13 |
| S15, S15b | mechanotransduction programme scores and threshold tests | 15b |
| S16, S16b–d | GWAS loci (curated; GWAS Catalog by trait; per locus) | 15 |
| S17a–i | F4 signature: per-gene step, GSEA, TF activity, positional organisation, cell of origin, candidate timing, hepatocyte programme | 15b |
| S18a–h | snRNA-seq: composition, within-population effects, concordance, hepatocyte class and targets, coupling, F4 signature, co-expression | 18 |

### Column glossary (Spanish labels kept for provenance)

The pipeline was developed in Spanish and some column values in the generated tables retain Spanish labels. They are
stable identifiers, not free text:

| In the tables | Meaning |
| --- | --- |
| `estadio` | condition / fibrosis stage (`Normal`, `Control`, `F0`–`F4`) |
| `orden` | ordinal stage, 0 (Normal/Control) to 5 (F4) |
| `cohorte`, `muestra`, `sexo`, `edad` | cohort, sample, sex, age |
| `Hepatocito`, `HSC`, `Colangiocito`, `Macrofago`, `Linfocito`, `Endotelio` | hepatocyte, stellate/myofibroblast, cholangiocyte, macrophage, lymphocyte, endothelium |
| `Mesenquima_HSC`, `Fagocito_mononuclear` | single-cell populations: mesenchyme/stellate, mononuclear phagocyte |
| `banda`, `periodo` | period band, period in genes |

The positional grid is the union of five per-cohort grid files; two of them (GSE142530, GSE276114) come from cohorts
not otherwise analysed but contribute gene slots, so all five files are required to reproduce the grid.

### Notes on reproducibility

- Two clean end-to-end re-runs reproduce the released tables; key values include 13,816 genes, 186 universal peaks,
  422 neighbourhoods, a ~80% reduction in the variance of stage statistics after composition adjustment, and an offset
  stage t of −8.1 → +1.3 (R² 0.15 → 0.65) after composition (Table S6d).
- `pyannotables` ships Ensembl 100 GRCh37; TADs are hg19; both are GRCh37-consistent.
- The Hallmark GMT is v7.0 from a public mirror; for submission, re-run step 07 with the official MSigDB release.
- Gene-list strings in a few tables (e.g. S10g, S16) may differ in order between runs; values do not.
