#!/usr/bin/env bash
# =====================================================================
# Script 07:
#   Part A: CAZyme gene clusters (CGCs = candidate PULs) and their
#           predicted substrates for ZS1B6-8, reusing the genes that
#           dbCAN predicted in script 06
#   Part B: CAZyme annotation of three reference genomes, for comparison
# =====================================================================
#SBATCH -A C3SE408-25-2          # project account
#SBATCH -J dbcan_cgc_refs        # job name
#SBATCH -t 03:00:00              # maximum run time
#SBATCH -n 8                     # CPU cores
#SBATCH -o logs/%x_%j.out        # log file: job name + job ID

DBCAN=containers/dbcan.sif
DB=databases/dbcan

# ---- Part A: gene clusters for our genome ---------------------------
# --mode protein   = input is proteins (uniInput.faa from script 06)
# --input_gff      = positions of those same genes (uniInput.gff),
#                    needed to see which genes lie next to each other
# --gff_type prodigal = the GFF was made by Prodigal (inside dbCAN)
echo "=== Part A: CGCs + substrates for ZS1B6-8 ==="
apptainer exec $DBCAN run_dbcan easy_substrate \
    --input_raw_data annotation/dbcan/uniInput.faa \
    --mode protein \
    --input_gff annotation/dbcan/uniInput.gff \
    --gff_type prodigal \
    --output_dir annotation/dbcan_substrate \
    --db_dir $DB \
    --threads 8

# ---- Part B: CAZymes in the three reference genomes -----------------
# A list of "short_name:NCBI_accession" pairs, one per genome
REFS="W_marina_type:GCA_014636515.1 W_fucanilytica:GCA_001697185.1 W_sp_IMCC45467:GCA_054111205.1"

for pair in $REFS; do
    name=${pair%%:*}       # text before the colon  -> short name
    acc=${pair##*:}        # text after the colon   -> accession
    fna=$(ls refs/ncbi_dataset/data/$acc/*.fna)   # that genome's DNA file
    echo "=== Part B: $name ($acc) ==="
    # Same settings as for our genome in script 06, so the counts are
    # directly comparable
    apptainer exec $DBCAN run_dbcan CAZyme_annotation \
        --input_raw_data $fna \
        --mode prok \
        --output_dir annotation/dbcan_refs/$name \
        --db_dir $DB \
        --threads 8
done

echo "=== All done ==="
