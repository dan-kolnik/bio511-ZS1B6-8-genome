#!/usr/bin/env bash
#SBATCH -A C3SE408-25-2
#SBATCH -J prokka
#SBATCH -t 02:00:00
#SBATCH -n 8
#SBATCH -o logs/%x_%j.out

C=/cephyr/NOBACKUP/groups/n2bin_gu/BIO511-2026/containers/prokka.sif

apptainer exec $C prokka \
    --outdir annotation/prokka \
    --prefix ZS1B6-8 \
    --locustag ZS1B6 \
    --kingdom Bacteria \
    --genus Wenyingzhuangia \
    --cpus 8 \
    assembly/medaka/consensus.fasta
