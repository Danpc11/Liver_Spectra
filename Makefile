PY ?= python3
S  := scripts
STEPS := 01_prepare_expression 02_spectra 03_specparam 04_differential_expression \
         05_composition_distance_tads 06_single_cell 07_gsea_tf 08_targets_fingerprint \
         09_make_figures 10_make_circos 11_mouse_validation 12_eqtl_shared_variants \
         13_gtex_tissues 14_eqtl_gtex_signif_pairs 15_gwas_loci_neighbourhoods 15b_mechanical_switch 18_snrnaseq_gse202379 16_figures_jhep 17_supplementary_figures

.PHONY: all external analysis figures audit readme clean help
help:
	@echo "make external  — download the public external resources (see data/external/README.md)"
	@echo "make all       — full pipeline, raw inputs to tables and figures (~35 min)"
	@echo "make analysis  — steps 01-08, 11-15b and 18 (tables only)"
	@echo "make figures   — steps 09, 10 and 16 (needs results/intermediate from a previous run)"
	@echo "make audit     — re-check the numbers quoted in the manuscripts"
	@echo "make clean     — delete generated results"

all:
	@set -e; cd $(S); for s in $(STEPS); do echo ">>> $$s"; $(PY) $$s.py; done

analysis:
	@set -e; cd $(S); for s in 01_prepare_expression 02_spectra 03_specparam 04_differential_expression 05_composition_distance_tads 06_single_cell 07_gsea_tf 08_targets_fingerprint 11_mouse_validation 12_eqtl_shared_variants 13_gtex_tissues 14_eqtl_gtex_signif_pairs 15_gwas_loci_neighbourhoods 15b_mechanical_switch 18_snrnaseq_gse202379; do echo ">>> $$s"; $(PY) $$s.py; done

# The figure scripts read results/intermediate/*.pkl, so run `make all` (or `make analysis`) at least
# once before `make figures`; `make clean` removes those intermediates.
figures:
	@test -f results/intermediate/pairs.pkl || { echo "results/intermediate is empty — run 'make analysis' first"; exit 1; }
	@set -e; cd $(S); for s in 09_make_figures 10_make_circos 16_figures_jhep 17_supplementary_figures; do echo ">>> $$s"; $(PY) $$s.py; done

external:
	bash $(S)/00_fetch_external.sh

readme:
	cd $(S) && $(PY) 99_readme_overview.py

audit:
	cd $(S) && $(PY) audit_numbers.py

clean:
	rm -rf results/intermediate/* results/tables/* results/figures/genome_research/* results/figures/jhep/*
	touch results/intermediate/.gitkeep
