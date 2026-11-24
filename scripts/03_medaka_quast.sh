#!/usr/bin/env bash
#SBATCH -A C3SE408-25-2
#SBATCH -J medaka_quast
#SBATCH -t 03:00:00
#SBATCH -n 8
#SBATCH -o logs/%x_%j.out

MEDAKA=containers/medaka.sif
QUAST=/cephyr/NOBACKUP/groups/n2bin_gu/BIO511-2026/containers/quast.sif

# Polish the Flye assembly with the raw reads
apptainer exec $MEDAKA medaka_consensus \
    -i reads/ZS1B6-8_ONT.fastq.gz \
    -d assembly/flye/assembly.fasta \
    -o assembly/medaka \
    -t 8 \
    --bacteria

# Statistics for draft vs. polished assembly
apptainer exec $QUAST quast.py \
    assembly/flye/assembly.fasta \
    assembly/medaka/consensus.fasta \
    -l "Flye_draft,Medaka_polished" \
    -o qc/quast \
    -t 8
