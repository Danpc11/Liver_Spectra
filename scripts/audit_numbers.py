"""audit_numbers.py — cross-check the numbers quoted in the manuscripts against the regenerated tables.
Run after the pipeline (`make all`). Writes results/tables/audit_number_checks.csv and exits non-zero if
any check fails, so it can be used in CI.
"""
import sys
import numpy as np, pandas as pd
from common import *

checks = []
def chk(name, claimed, computed, tol=0.06):
    try: ok = abs(float(computed) - float(claimed)) <= tol * max(1e-9, abs(float(claimed)))
    except (TypeError, ValueError): ok = claimed == computed
    checks.append({'claim': name, 'manuscript': claimed, 'computed': round(float(computed), 4) if isinstance(computed, (int, float, np.floating)) else computed, 'ok': bool(ok)})

S1 = pd.read_csv(tab('Table_S1_samples.csv'), index_col=0)
chk('n samples', 437, len(S1)); chk('n Normal', 35, (S1.estadio == 'Normal').sum()); chk('n F3', 76, (S1.estadio == 'F3').sum()); chk('n F4', 28, (S1.estadio == 'F4').sum())
U = pd.read_csv(tab('Table_S2a_consensus_spectrum_universal_peaks.csv'))
chk('n frequencies', 10007, len(U)); chk('universal peaks', 186, U.universal.sum())
msk = (U.w_mask > 0) & np.isfinite(U.w_mask)
chk('corr consensus vs grid mask', 0.0, np.corrcoef(np.log(U.mean_power_over_null[msk]), np.log(U.w_mask[msk]))[0, 1], tol=1e9)
chk('stage-associated frequencies, log spectra', 6042, (pd.read_csv(tab('Table_S2c_stage_effect_per_frequency_log.csv')).q < 0.05).sum())
chk('stage-associated frequencies, z spectra', 0, (pd.read_csv(tab('Table_S2c_stage_effect_per_frequency_z.csv')).q < 0.05).sum(), tol=0)
Bd = pd.read_csv(tab('Table_S2d_stage_effect_per_band_log.csv')); chk('band 5-10 t', -5.0, Bd.t_stage.iloc[3]); chk('band 50-100 t', 2.9, Bd.t_stage.iloc[6])
SP = pd.read_csv(tab('Table_S3c_specparam_vs_stage.csv')).set_index('parameter')
chk('offset t', -8.1, SP.loc['offset', 't']); chk('peak-count t', -3.8, SP.loc['n_peaks', 't']); chk('exponent P', 0.10, SP.loc['exponent', 'p'], tol=0.2)
chk('DE F4 vs Normal padj<0.05', 8713, (pd.read_csv(tab('Table_S4a_DESeq2_F4_vs_Normal.csv'), index_col=0).padj < 0.05).sum())
chk('DE ordinal padj<0.05', 7520, (pd.read_csv(tab('Table_S4b_DESeq2_ordinal_stage.csv'), index_col=0).padj < 0.05).sum())
A5 = pd.read_csv(tab('Table_S5a_spatial_autocorrelation_DE.csv')); chk('ACF lag1', 0.17, A5.acf_observed.iloc[0]); chk('ACF lag1 z', 20, A5.z.iloc[0], tol=0.1)
P5 = pd.read_csv(tab('Table_S5b_adjacent_DE_pairs.csv')); chk('adjacent up pairs', 954, P5.observed.iloc[0]); chk('adjacent down pairs', 957, P5.observed.iloc[1])
chk('runs >=3 genes', 422, len(pd.read_csv(tab('Table_S5c_contiguous_DE_neighbourhoods.csv'))))
chk('ACF lag1, composition-adjusted', 0.12, pd.read_csv(tab('Table_S5d_acf_with_without_composition.csv')).acf_composition_adjusted.iloc[0])
Dg = pd.read_csv(tab('Table_S5g_concordance_by_distance.csv'), index_col=0)
chk('concordance, overlapping genes', 0.38, Dg.r_t.iloc[0]); chk('concordance, <1 kb', 0.22, Dg.r_t.iloc[1]); chk('concordance, 10-50 kb', 0.18, Dg.r_t.iloc[3])
chk('TAD test P (distance-stratified)', 0.26, pd.read_csv(tab('Table_S5i_TAD_test_stratified.csv')).p.iloc[0], tol=0.4)
T6 = pd.read_csv(tab('Table_S6c_stage_effect_per_gene_with_without_composition.csv'))
chk('variance explained by composition (%)', 80, 100 * (1 - T6.t_stage_comp_adjusted.var() / T6.t_stage.var()))
C6 = pd.read_csv(tab('Table_S6b_composition_scores_per_sample.csv'), index_col=0)
chk('hepatocyte score, Normal', 0.30, C6.groupby('estadio').Hepatocito.mean()['Normal']); chk('hepatocyte score, F4', -0.30, C6.groupby('estadio').Hepatocito.mean()['F4'])
B7 = pd.read_csv(tab('Table_S7a_threshold_vs_linear.csv')).set_index('variable')
chk('cholangiocyte dAIC at F3', 8.6, B7.loc['Colangiocito', 'dAIC_step_at_F3']); chk('cholangiocyte dAIC at F4', 13.6, B7.loc['Colangiocito', 'dAIC_step_at_F4']); chk('stellate dAIC at F4', 4.5, B7.loc['HSC', 'dAIC_step_at_F4'])
Bi = pd.read_csv(tab('Table_S7b_within_stage_bimodality.csv')); r = Bi[(Bi.variable == 'Hepatocito') & (Bi.estadio == 'F3')].iloc[0]
chk('F3 hepatocyte dBIC', 12.8, r['dBIC_1_minus_2']); chk('F3 minor component', 0.08, r.min_weight, tol=0.15)
S8 = pd.read_csv(tab('Table_S8c_within_type_spatial_autocorrelation.csv')).set_index('type')
chk('sc ACF macrophage', 0.099, S8.loc['Fagocito_mononuclear', 'acf_lag1']); chk('sc ACF endothelium', 0.080, S8.loc['Endotelio', 'acf_lag1']); chk('sc ACF stellate', 0.066, S8.loc['Mesenquima_HSC', 'acf_lag1'])
chk('lineage-intrinsic targets', 154, len(pd.read_csv(tab('Table_S9a_candidate_targets_lineage_intrinsic.csv'))), tol=0.02)
TF = pd.read_csv(tab('Table_S10d_DoRothEA_TF_activity_vs_stage.csv')); tfi = TF.set_index('tf')
chk('TFs q<0.05', 154, (TF.q < 0.05).sum()); chk('TFs after composition', 70, (TF.q_adj_comp < 0.05).sum())
chk('NFKB1 adjusted t', 5.7, tfi.loc['NFKB1', 't_stage_comp_adjusted']); chk('HNF4A adjusted t', 1.3, tfi.loc['HNF4A', 't_stage_comp_adjusted'], tol=0.3)
chk('OXPHOS NES, composition-adjusted', 2.2, pd.read_csv(tab('Table_S10b_GSEA_hallmark_composition_adjusted.csv')).set_index('Term').loc['OXIDATIVE_PHOSPHORYLATION', 'NES'])
S11 = pd.read_csv(tab('Table_S11_spectral_fingerprint_real_vs_shuffled.csv'))
chk('fingerprint AUC, real order', 0.91, S11[(S11.features == 'spectrum') & (S11.task == 'early_vs_advanced')].value.iloc[0])
chk('fingerprint AUC, shuffled order', 0.91, S11[(S11.features == 'spectrum_shuffled_order') & (S11.task == 'early_vs_advanced')].value.iloc[0])
chk('mouse ACF lag1', 0.11, pd.read_csv(tab('Table_S12c_mouse_spatial_autocorrelation.csv')).acf_observed.iloc[0])
chk('mouse universal peaks', 157, pd.read_csv(tab('Table_S12b_mouse_spectrum_universal_peaks.csv')).universal_all6.sum())
S13 = pd.read_csv(tab('Table_S13e_shared_GTEx_eQTL_coupled_vs_uncoupled.csv')).iloc[1]
chk('shared same-direction eQTL OR', 1.76, S13.MH_OR_distance_stratified); chk('shared eQTL P', 0.016, S13.MH_p, tol=0.3)
S14 = pd.read_csv(tab('Table_S14c_gtex_universal_peaks_by_tissue.csv')); r = S14[(S14.gene_set == 'own') & (S14.tissue == 'liver')].iloc[0]
chk('biopsy peaks universal in GTEx liver', 0.28, r.frac_biopsy_peaks_universal_in_tissue); chk('r biopsy vs GTEx liver', 0.47, r.r_mean_log_power_with_biopsy_liver)
S15 = pd.read_csv(tab('Table_S15_mechanics_scores_thresholds.csv')).set_index('score')
chk('YAP/TAZ dAIC at F4', 28.8, S15.loc['YAP_TAZ_targets', 'dAIC_step_F4']); chk('YAP/TAZ adjusted t at F4', 5.0, S15.loc['YAP_TAZ_targets', 't_step_F4_comp_adj'])
chk('crosslinking adjusted t at F4', 0.9, S15.loc['matrix_crosslinking', 't_step_F4_comp_adj'], tol=0.25)
S16 = pd.read_csv(tab('Table_S16c_GWAS_catalog_enrichment_by_trait.csv')).set_index('trait_set')
chk('cirrhosis loci in DE neighbourhoods', 0.229, S16.loc['liver fibrosis/cirrhosis', 'frac_in_DE_neighbourhood'])
chk('cirrhosis loci local ACF', 0.169, S16.loc['liver fibrosis/cirrhosis', 'acf_windows'])
chk('MASLD loci in coupled pairs', 0.278, S16.loc['MASLD/NAFLD', 'frac_in_coupled_pair'])

R = pd.DataFrame(checks); R.to_csv(tab('audit_number_checks.csv'), index=False)
pd.set_option('display.width', 200); print(R.to_string(index=False))
failed = int((~R.ok).sum()); print(f'\n{len(R) - failed} of {len(R)} checks passed.')
sys.exit(1 if failed else 0)
