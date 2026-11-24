#!/usr/bin/env bash
# =====================================================================
# Script 06: carbohydrate-active enzymes (CAZymes) with dbCAN 5.2.9
#   Part 1: CAZyme annotation (always runs)
#   Part 2: CAZyme gene clusters (CGCs, candidate PULs) + substrate
#           prediction (only runs if the STP.hmm database file exists)
# =====================================================================
#SBATCH -A C3SE408-25-2          # project account
#SBATCH -J dbcan                 # job name
#SBATCH -t 03:00:00              # maximum run time
#SBATCH -n 8                     # CPU cores
#SBATCH -o logs/%x_%j.out        # log file: job name + job ID

GENOME=assembly/medaka/consensus.fasta     # our polished genome (DNA)
DBCAN=containers/dbcan.sif
DB=databases/dbcan

# ---- Part 1: CAZyme annotation -------------------------------------
# dbCAN predicts the genes itself (Prodigal, --mode prok), then assigns
# each protein to CAZy families (GH, GT, PL, CE, AA, CBM) with three
# methods: HMMER (dbCAN.hmm), DIAMOND (CAZy.dmnd) and dbCAN-sub.
# Main result: annotation/dbcan/overview.tsv (one row per CAZyme).
echo "=== Part 1: CAZyme annotation ==="
apptainer exec $DBCAN run_dbcan CAZyme_annotation \
    --input_raw_data $GENOME \
    --mode prok \
    --output_dir annotation/dbcan \
    --db_dir $DB \
    --threads 8

# ---- Part 2: gene clusters + substrates (needs STP.hmm) ------------
# [ -f file ] is true if the file exists.
if [ -f $DB/STP.hmm ]; then
    echo "=== Part 2: CGCs and substrate prediction ==="
    apptainer exec $DBCAN run_dbcan easy_substrate \
        --input_raw_data $GENOME \
        --mode prok \
        --output_dir annotation/dbcan_substrate \
        --db_dir $DB \
        --input_gff gff \
        --gff_type prodigal \
        --threads 8
else
    echo "=== Part 2 skipped: $DB/STP.hmm is missing ==="
fi

echo "=== All done ==="
