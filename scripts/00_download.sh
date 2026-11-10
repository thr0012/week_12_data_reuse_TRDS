#!/usr/bin/env bash
# Download public inputs for the DMR escape analysis (run on Easley or any machine with internet).
# Usage: bash 00_download.sh [data_dir]
set -euo pipefail

DATA=${1:-data}
mkdir -p "$DATA"/{genome,annot,methyl,dmrs}

UCSC=https://hgdownload.soe.ucsc.edu/goldenPath/rn6

# 1) Genome (rn6 = Rnor_6.0, the assembly used by Thorson et al. 2020)
wget -nc -P "$DATA/genome" "$UCSC/bigZips/rn6.fa.gz"
wget -nc -P "$DATA/genome" "$UCSC/bigZips/rn6.chrom.sizes"

# 2) Annotation tracks (UCSC table dumps). Missing tracks are skipped, not fatal.
for t in rmsk gap genomicSuperDups simpleRepeat cpgIslandExt; do
  wget -nc -P "$DATA/annot" "$UCSC/database/${t}.txt.gz" || echo "WARN: ${t} not available for rn6, skipping"
done

# 3) Germ cell methylomes, El Omri-Charai et al. 2023 (GSE220440)
#    List the processed files first, then download only GD18 gonocytes and adult sperm.
GEO=https://ftp.ncbi.nlm.nih.gov/geo/series/GSE220nnn/GSE220440
wget -q -O "$DATA/methyl/GSE220440_suppl_listing.html" "$GEO/suppl/" || true
echo "Series-level files:"; grep -o 'href="[^"]*"' "$DATA/methyl/GSE220440_suppl_listing.html" | cut -d'"' -f2 | grep -v '^\?' || true
wget -q -O "$DATA/methyl/GSE220440_family.soft.gz" "$GEO/soft/GSE220440_family.soft.gz" || true
echo
echo "Sample titles (pick GD18 gonocyte and adult sperm GSMs):"
zcat "$DATA/methyl/GSE220440_family.soft.gz" 2>/dev/null | grep -E '^\^SAMPLE|^!Sample_title|^!Sample_supplementary_file' | paste - - - || true
echo
echo "Next: download the GD18 and sperm sample files into $DATA/methyl/gd18 and $DATA/methyl/sperm, e.g."
echo "  wget -P $DATA/methyl/gd18 <ftp link from !Sample_supplementary_file>"
echo
echo "4) DMR tables: download by hand from the article pages into $DATA/dmrs/"
echo "   McBirney et al. 2017 https://doi.org/10.1371/journal.pone.0184306 (F3 atrazine vs control DMR table)"
echo "   Thorson et al. 2020  https://doi.org/10.1371/journal.pone.0239380 (S2 to S6 Tables)"
