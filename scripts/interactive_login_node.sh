#!/usr/bin/env bash
# =====================================================================
# Commands run directly on the Vera LOGIN NODE (not submitted to SLURM).
# These are short tasks or tasks that need internet access, which the
# compute nodes may not have. Kept here as a record; not meant to be run
# as one script. Run from the project folder.
# =====================================================================

# ---- 0. Project folder and link to the course reads ----------------
mkdir -p reads qc assembly annotation logs containers scripts refs databases
ln -s /cephyr/NOBACKUP/groups/n2bin_gu/BIO511-2026/data/ONT_BIO511/ZS1B6-8_ONT.fastq.gz reads/

# ---- 1. Build containers from the recipes in containers/ -----------
cd containers
for t in ont-tools medaka skani checkm2 dbcan barrnap; do
    apptainer build $t.sif $t.def          # long builds were run with nohup ... &
done
cd ..

# ---- 2. Tool versions (for the Methods section) --------------------
F=/cephyr/NOBACKUP/groups/n2bin_gu/BIO511-2026/containers
apptainer exec $F/flye.sif flye --version                  # 2.9.5-b1801
apptainer exec $F/quast.sif quast.py --version             # 5.3.0
apptainer exec $F/prokka.sif prokka --version              # 1.14.6
apptainer exec containers/ont-tools.sif seqkit version     # 2.14.0
apptainer exec containers/ont-tools.sif NanoPlot --version # 1.48.0
apptainer exec containers/medaka.sif medaka --version      # 2.2.2
apptainer exec containers/skani.sif skani --version        # 0.4.0
apptainer exec containers/checkm2.sif checkm2 --version    # 1.1.0
apptainer exec containers/dbcan.sif ls /opt/conda/conda-meta | grep -i dbcan   # 5.2.9

# ---- 3. Reference genomes from NCBI (needs internet) ---------------
apptainer exec containers/skani.sif datasets download genome taxon Wenyingzhuangia \
    --include genome --filename refs/wenyingzhuangia.zip
cd refs && unzip -q wenyingzhuangia.zip && cd ..
apptainer exec containers/skani.sif dataformat tsv genome \
    --inputfile refs/ncbi_dataset/data/assembly_data_report.jsonl \
    --fields accession,organism-name,organism-infraspecific-strain,assminfo-level,assmstats-total-sequence-len,assmstats-gc-percent \
    | grep -E "^Assembly|^GCA_" > refs/reference_genomes.tsv

# ---- 4. Databases (needs internet) ----------------------------------
apptainer exec containers/checkm2.sif checkm2 database --download --path databases/
apptainer exec containers/dbcan.sif run_dbcan database --db_dir databases/dbcan
# STP.hmm failed from the primary server (403); fetched from the AWS mirror:
apptainer exec containers/dbcan.sif run_dbcan database --db_dir databases/dbcan --aws_s3

# ---- 5. rRNA genes with barrnap (Prokka's container lacked nhmmer) ---
mkdir -p annotation/barrnap
apptainer exec containers/barrnap.sif barrnap --kingdom bac --threads 4 \
    --outseq annotation/barrnap/rRNA.fasta \
    assembly/medaka/consensus.fasta > annotation/barrnap/rRNA.gff

# ---- 6. Per-read length and quality table (input for the read figure)
apptainer exec containers/ont-tools.sif seqkit fx2tab -n -i -l -q \
    reads/ZS1B6-8_ONT.fastq.gz | gzip > qc/read_len_qual.tsv.gz

# ---- 7. Checking jobs --------------------------------------------------
# sbatch scripts/0X_name.sh            submit a job
# squeue -u <user>                      waiting / running jobs
# sacct -j <jobid> --format=JobID,JobName,State,Elapsed,MaxRSS
