#!/usr/bin/env bash
#SBATCH -A C3SE408-25-2
#SBATCH -J flye
#SBATCH -t 06:00:00
#SBATCH -n 16
#SBATCH -o logs/%x_%j.out

C=/cephyr/NOBACKUP/groups/n2bin_gu/BIO511-2026/containers/flye.sif

apptainer exec $C flye \
    --nano-hq reads/ZS1B6-8_ONT.fastq.gz \
    --genome-size 3.25m \
    --asm-coverage 50 \
    --out-dir assembly/flye \
    --threads 16
