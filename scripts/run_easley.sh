#!/usr/bin/env bash
#SBATCH --job-name=dmr_escape
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=03:00:00
#SBATCH --output=dmr_escape_%j.log
# Run from the repository root:  sbatch scripts/run_easley.sh [all]
# Add --partition / --account lines if your allocation needs them.
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"

D=data
COMMON="--genome $D/genome/rn6.fa.gz --gap $D/annot/gap.txt.gz --rmsk $D/annot/rmsk.txt.gz \
  --simple $D/annot/simpleRepeat.txt.gz --cpgi $D/annot/cpgIslandExt.txt.gz \
  --gd18 $D/methyl/gd18/*.cov.gz --sperm $D/methyl/sperm/*.cov.gz --n-perm 1000"
SETS=${1:-F3}
run () { python3 scripts/dmr_escape.py $COMMON "$@"; }

run --dmrs $D/dmrs/McBirney2017_F3_DMRs.tsv --label "McBirney 2017 F3 (atrazine vs control)" --out results/McBirney_F3

if [ "$SETS" = "all" ]; then
  run --dmrs $D/dmrs/McBirney2017_F1_DMRs.tsv --label "McBirney 2017 F1" --out results/McBirney_F1
  run --dmrs $D/dmrs/McBirney2017_F2_DMRs.tsv --label "McBirney 2017 F2" --out results/McBirney_F2
  for t in S2_lean S3_kidney S4_testis S5_puberty S6_multiple; do
    run --dmrs $D/dmrs/Thorson2020_${t}_DMRs.tsv --label "Thorson 2020 ${t}" --out results/Thorson_${t}
  done
  # robustness: repeats must cover >= 50% of DMR bp
  run --dmrs $D/dmrs/McBirney2017_F3_DMRs.tsv --label "McBirney F3, overlap >= 50%" --min-overlap 0.5 \
      --out robustness/McBirney_F3_overlap50
fi
python3 scripts/04_compare.py results
