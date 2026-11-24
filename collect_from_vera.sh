#!/usr/bin/env bash
# Copy the small result files and SLURM logs from the project folder on Vera
# into this repository. Run on Vera, from inside the repository folder:
#   bash collect_from_vera.sh /cephyr/NOBACKUP/groups/n2bin_gu/students/<you>/bio511-genom-project
set -euo pipefail
P=${1:?give the path to your project folder}

mkdir -p results/qc results/assembly results/annotation results/taxonomy logs
cp "$P"/scripts/*.sh scripts/                       # your exact job scripts (overwrite the copies here)
cp "$P"/containers/*.def containers/                # your exact container recipes
cp "$P"/logs/*.out logs/                             # SLURM job logs (show job IDs)
cp "$P"/qc/seqkit_raw.tsv results/qc/
cp "$P"/qc/nanoplot_raw/NanoStats.txt results/qc/
cp "$P"/qc/quast/report.txt results/qc/quast_report.txt
cp "$P"/qc/checkm2/quality_report.tsv results/qc/checkm2_quality_report.tsv
cp "$P"/qc/skani_ani.tsv results/taxonomy/
cp "$P"/refs/reference_genomes.tsv results/taxonomy/
cp "$P"/assembly/flye/assembly_info.txt results/assembly/flye_assembly_info.txt
cp "$P"/assembly/medaka/consensus.fasta results/assembly/ZS1B6-8_polished.fasta
cp "$P"/annotation/prokka/ZS1B6-8.txt results/annotation/prokka_summary.txt
cp "$P"/annotation/prokka/ZS1B6-8.gff results/annotation/
cp "$P"/annotation/barrnap/rRNA.gff results/annotation/barrnap_rRNA.gff
echo "Done. Check with: git status"
