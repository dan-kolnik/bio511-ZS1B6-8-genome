#!/usr/bin/env bash
# =====================================================================
# Script 05: species identity (ANI) and genome quality (CheckM2)
#            for strain ZS1B6-8
# =====================================================================
#SBATCH -A C3SE408-25-2          # project account that pays for the job
#SBATCH -J ani_checkm2           # job name (shown in squeue and in log name)
#SBATCH -t 02:00:00              # maximum run time (hh:mm:ss)
#SBATCH -n 8                     # number of CPU cores
#SBATCH -o logs/%x_%j.out        # log file: %x = job name, %j = job ID

# ---- File and container paths (defined once, reused below) ----------
GENOME=assembly/medaka/consensus.fasta        # our polished assembly
SKANI=containers/skani.sif
CHECKM2=containers/checkm2.sif

# ---- 1. ANI: how similar is our genome to each reference genome? ----
# skani compares whole genomes and reports Average Nucleotide Identity.
#   ANI >= 95%  -> same species;  below 95% -> different species
#   -q    = query (our genome)
#   -r    = reference genomes; GCA_* only, because each GCF_* is a duplicate
#   -s 70 = also report weaker matches (default only reports ANI > 80%)
#   --ci  = add confidence intervals for the ANI values
echo "=== skani ANI ==="
apptainer exec $SKANI skani dist \
    -q $GENOME \
    -r refs/ncbi_dataset/data/GCA_*/*.fna \
    -t 8 -s 70 --ci \
    -o qc/skani_ani.tsv

# ---- 2. CheckM2: genome completeness and contamination --------------
# Completeness  = share of expected marker genes that are present (%)
# Contamination = share of marker genes found more than once (%)
# High-quality genome: completeness > 90% and contamination < 5%
#   --env GLIBC_TUNABLES=... prevents a known TensorFlow "executable stack"
#     error on newer systems; it has no effect otherwise
#   --database_path = where the database is (needed because the container
#     could not save this location itself during the download)
echo "=== CheckM2 ==="
apptainer exec --env GLIBC_TUNABLES=glibc.rtld.execstack=2 $CHECKM2 \
    checkm2 predict \
    --input $GENOME \
    --output-directory qc/checkm2 \
    --database_path databases/CheckM2_database/uniref100.KO.1.dmnd \
    --threads 8 \
    --force                       # overwrite the output folder if it exists

echo "=== All done ==="
