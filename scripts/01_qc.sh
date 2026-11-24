#!/usr/bin/env bash
#SBATCH -A C3SE408-25-2
#SBATCH -J qc
#SBATCH -t 01:00:00
#SBATCH -n 4
#SBATCH -o logs/%x_%j.out

C=containers/ont-tools.sif
R=reads/ZS1B6-8_ONT.fastq.gz

apptainer exec $C seqkit stats -a $R > qc/seqkit_raw.tsv
apptainer exec $C NanoPlot --fastq $R -o qc/nanoplot_raw -t 4
