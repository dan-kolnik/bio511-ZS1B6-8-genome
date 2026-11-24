# Complete genome of the marine flavobacterium *Wenyingzhuangia* sp. ZS1B6-8

Final project in **BIO511 Bioinformatics in Genomics**, Master's Programme in Bioinformatics, University of Gothenburg (autumn 2026).

Author: Dan Kolnik

This repository contains the job scripts, container recipes, Python scripts and main results for the analysis of strain ZS1B6-8, a marine isolate first identified as *Wenyingzhuangia marina* by 16S rRNA Sanger sequencing. The genome was sequenced with Oxford Nanopore. It was then assembled into one circular chromosome, evaluated, annotated, compared with all available *Wenyingzhuangia* genomes, and searched for carbohydrate-active enzymes (CAZymes).

## Main results

| | |
|---|---|
| Reads (barcode06, pass) | 273,383 reads, 545.7 Mb, read N50 4.4 kb |
| Assembly | 1 circular contig, 3,831,794 bp, 139× coverage, GC 31.2% |
| Quality (CheckM2) | 99.99% complete, 0.36% contamination |
| Annotation | 3,122 protein-coding genes, 53 tRNAs, 4 rRNA operons |
| Species identity | highest ANI 82.96% (*W.* sp. IMCC45467), so most likely a novel species |
| CAZymes (dbCAN, ≥ 2 of 3 methods) | 235 CAZymes, 61 CAZyme gene clusters |

## Repository structure

```
bio511-ZS1B6-8-genome/
├── README.md
├── collect_from_vera.sh          # copies results + logs from the Vera project folder into this repo
├── scripts/
│   ├── 01_qc.sh                  # SLURM: SeqKit + NanoPlot
│   ├── 02_flye.sh                # SLURM: Flye assembly
│   ├── 03_medaka_quast.sh        # SLURM: Medaka polishing + QUAST
│   ├── 04_prokka.sh              # SLURM: Prokka annotation
│   ├── 05_ani_checkm2.sh         # SLURM: skani ANI + CheckM2
│   ├── 06_dbcan.sh               # SLURM: dbCAN CAZyme annotation
│   ├── 07_dbcan_cgc_refs.sh      # SLURM: dbCAN gene clusters + reference genomes
│   └── interactive_login_node.sh # record of commands run on the login node
├── containers/                   # Apptainer recipes (.def) for tools not provided by the course
├── python/
│   ├── cazyme_summary.py         # summarises dbCAN results, target enzymes, figures
│   └── make_read_figure.py       # read length / quality figure
├── results/
│   ├── qc/  assembly/  annotation/  taxonomy/  cazymes/  figures/
└── logs/                         # SLURM job logs
```

Raw reads, containers (`.sif`) and databases are not included because of their size. The raw reads are course data (`ONT_BIO511/ZS1B6-8_ONT.fastq.gz`).

## Computing environment

All analyses were run on the **Vera** cluster (C3SE, Chalmers). Every tool ran inside an **Apptainer** container. Flye, QUAST and Prokka used the course-provided containers. The other tools were built from Bioconda with the recipes in `containers/`.

### How each step was run

| Step | Tool (version) | Script | How it was run |
|---|---|---|---|
| Read QC | SeqKit 2.14.0, NanoPlot 1.48.0 | `01_qc.sh` | SLURM job 10547125 |
| Assembly | Flye 2.9.5-b1801 | `02_flye.sh` | SLURM job 10546322 |
| Polishing + stats | Medaka 2.2.2, QUAST 5.3.0 | `03_medaka_quast.sh` | SLURM job 10547198 |
| Annotation | Prokka 1.14.6 | `04_prokka.sh` | SLURM job 10548543 |
| ANI + completeness | skani 0.4.0, CheckM2 1.1.0 | `05_ani_checkm2.sh` | SLURM job 10551258 |
| CAZyme annotation | run_dbcan 5.2.9 | `06_dbcan.sh` | SLURM job 10551346 |
| Gene clusters + references | run_dbcan 5.2.9 | `07_dbcan_cgc_refs.sh` | SLURM job 10596331 |
| Container builds, downloads, version checks | Apptainer, NCBI Datasets | `interactive_login_node.sh` | Login node (downloads need internet) |
| rRNA genes | barrnap 0.9 | `interactive_login_node.sh` | Login node (short task) |
| Per-read length/quality table | SeqKit `fx2tab` | `interactive_login_node.sh` | Login node (short task) |
| CAZyme summary + figures | Python (pandas, seaborn) | `python/cazyme_summary.py` | Laptop, conda env `bio511py` |
| Read figure | Python (pandas, matplotlib) | `python/make_read_figure.py` | See "Use of AI" below |

## How to reproduce

From the project folder on Vera (with `reads/`, `containers/` and `databases/` in place, see `interactive_login_node.sh`):

```bash
sbatch scripts/01_qc.sh
sbatch scripts/02_flye.sh
sbatch scripts/03_medaka_quast.sh      # after 02
sbatch scripts/04_prokka.sh            # after 03
sbatch scripts/05_ani_checkm2.sh       # after 03 + reference genomes + CheckM2 database
sbatch scripts/06_dbcan.sh             # after 03 + dbCAN database
sbatch scripts/07_dbcan_cgc_refs.sh    # after 06
```

Then, with the dbCAN output copied to a computer with pandas, matplotlib and seaborn:

```bash
python python/cazyme_summary.py \
    --overview dbcan/overview.tsv \
               dbcan_refs/W_marina_type/overview.tsv \
               dbcan_refs/W_fucanilytica/overview.tsv \
               dbcan_refs/W_sp_IMCC45467/overview.tsv \
    --names ZS1B6-8 W_marina_type W_fucanilytica W_sp_IMCC45467 \
    --outdir results/cazymes

python python/make_read_figure.py --input qc/read_len_qual.tsv.gz --output results/figures/fig_reads.png
```

## Troubleshooting notes

- **Prokka reported 0 rRNA genes.** The course Prokka container lacks `nhmmer` (HMMER), which barrnap needs. barrnap was run in a separate container and found 4 complete 16S–23S–5S operons.
- **CheckM2 failed with a TensorFlow "cannot enable executable stack" error.** It was solved by running with `--env GLIBC_TUNABLES=glibc.rtld.execstack=2`. The database path is passed explicitly because the container is read-only.
- **The first dbCAN container** installed an old version (4.1.4) on Python 3.15 and lacked `requests`. It was rebuilt with `python=3.11 "dbcan>=5" requests`.
- **dbCAN database download:** `STP.hmm` was refused by the primary server (403) and was downloaded from the AWS mirror (`--aws_s3`).
- **Gene-cluster step in `06_dbcan.sh`:** `--input_gff gff` was read as a file name. It was redone in `07_dbcan_cgc_refs.sh` using the protein and GFF files from the first dbCAN step.
- **Barcode/adapter trimming** was off during basecalling and reads were not trimmed before assembly. See the report's discussion.

## Use of AI

Claude Opus 5.5 (Anthropic) was used throughout the project to:

- explain tools and outputs
- suggest the analysis plan and troubleshoot errors
- draft the job scripts, container recipes and both Python scripts
- produce the read figure (`results/figures/fig_reads.png`), later written up as `python/make_read_figure.py`
- draft the report and the presentation

All commands were run and all outputs were checked by the author.

## Main references

- Flye: Kolmogorov *et al.* (2019) *Nat Biotechnol* 37: 540–546
- Medaka: Oxford Nanopore Technologies, https://github.com/nanoporetech/medaka
- QUAST: Gurevich *et al.* (2013) *Bioinformatics* 29: 1072–1075
- CheckM2: Chklovski *et al.* (2023) *Nat Methods* 20: 1203–1212
- Prokka: Seemann (2014) *Bioinformatics* 30: 2068–2069
- barrnap: https://github.com/tseemann/barrnap
- skani: Shaw & Yu (2023) *Nat Methods* 20: 1661–1665
- dbCAN3: Zheng *et al.* (2023) *Nucleic Acids Res* 51: W115–W121
- SeqKit: Shen *et al.* (2016) *PLoS ONE* 11: e0163962
- NanoPlot: De Coster & Rademakers (2023) *Bioinformatics* 39: btad311
