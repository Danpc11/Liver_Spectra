PY=python3
S=scripts
.PHONY: all external figures clean
all: $(S)/01_prepare_expression.py
	cd $(S) && $(PY) 01_prepare_expression.py && $(PY) 02_spectra.py && $(PY) 03_specparam.py && $(PY) 04_differential_expression.py && $(PY) 05_composition_distance_tads.py && $(PY) 06_single_cell.py && $(PY) 07_gsea_tf.py && $(PY) 08_targets_fingerprint.py && $(PY) 09_make_figures.py && $(PY) 10_make_circos.py
external:
	bash $(S)/00_fetch_external.sh
figures:
	cd $(S) && $(PY) 09_make_figures.py && $(PY) 10_make_circos.py
clean:
	rm -rf results/intermediate/* results/tables/* results/figures/*
