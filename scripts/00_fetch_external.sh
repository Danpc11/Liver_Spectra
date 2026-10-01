#!/usr/bin/env bash
# Downloads public external resources. Run from repo root. Requires curl, unzip, tar.
set -euo pipefail
E=data/external; mkdir -p $E
# 1. DoRothEA human regulons (saezlab/dorothea, Bioconductor data package)
curl -L -o $E/dorothea.zip https://codeload.github.com/saezlab/dorothea/zip/refs/heads/master
unzip -o -q $E/dorothea.zip "dorothea-master/data/dorothea_hs.rda" -d $E && mv $E/dorothea-master/data/dorothea_hs.rda $E/ && rm -r $E/dorothea-master
# 2. MSigDB Hallmark v7.0 GMT (public mirror; replace with the official MSigDB download for the final version)
curl -L -o $E/erilu.zip https://codeload.github.com/erilu/bulk-rnaseq-analysis/zip/refs/heads/master
unzip -o -q $E/erilu.zip "bulk-rnaseq-analysis-master/h.all.v7.0.symbols.gmt.txt" -d $E && mv $E/bulk-rnaseq-analysis-master/h.all.v7.0.symbols.gmt.txt $E/hallmark.gmt && rm -r $E/bulk-rnaseq-analysis-master
# 3. Liver TAD partitions (Leung 2015 via McArthur & Capra; hg19)
curl -L -o $E/tad.zip https://codeload.github.com/emcarthur/TAD-stability-heritability/zip/refs/heads/master
unzip -o -q $E/tad.zip "TAD-stability-heritability-master/data/20binsTADlandscape/Liver_leung2015/*" -d $E
# 4. GSE136103 (Ramachandran 2019): download GSE136103_RAW.tar from GEO into data/external/ (436 MB), then extract:
mkdir -p $E/GSE136103
[ -f $E/GSE136103_RAW.tar ] && tar -xf $E/GSE136103_RAW.tar -C $E/GSE136103/ || echo "Place GSE136103_RAW.tar in $E/ and re-run."
echo done
