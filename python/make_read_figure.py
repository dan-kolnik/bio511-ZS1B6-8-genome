#!/usr/bin/env python3
"""
make_read_figure.py
BIO511 genomics project - strain ZS1B6-8

Draws the two-panel read figure (Figure 1 in the report, slide 5):
  a) read-length distribution weighted by bases ("where the bases are"),
     with the read N50 marked
  b) mean read quality against read length (hexagonal bins, log colour scale)

Input: a tab-separated table with one read per line: read ID, length, mean
quality. It was made on Vera with SeqKit:

  apptainer exec containers/ont-tools.sif seqkit fx2tab -n -i -l -q \
      reads/ZS1B6-8_ONT.fastq.gz | gzip > qc/read_len_qual.tsv.gz

Example:
  python make_read_figure.py --input qc/read_len_qual.tsv.gz --output qc/fig_reads.png
"""

import argparse
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")              # draw to a file, no screen needed
import matplotlib.pyplot as plt


def existing_file(file_path):
    """Check that the input file exists."""
    if not os.path.isfile(file_path):
        raise argparse.ArgumentTypeError(f"'{file_path}' is not a valid file")
    return file_path


def read_n50(lengths):
    """N50: length of the read at which the running total of bases, summed
    from the longest read downwards, reaches half of all bases."""
    sorted_lengths = np.sort(lengths)[::-1]
    running_total = np.cumsum(sorted_lengths)
    return sorted_lengths[np.searchsorted(running_total, running_total[-1] / 2)]


def main():
    parser = argparse.ArgumentParser(description="Plot read length and quality for one nanopore sample.")
    parser.add_argument("--input", type=existing_file, required=True,
                        help="seqkit fx2tab output (ID, length, quality); may be gzipped")
    parser.add_argument("--output", default="fig_reads.png", help="output image file")
    args = parser.parse_args()

    reads = pd.read_csv(args.input, sep="\t", header=None, names=["id", "length", "quality"])
    n50 = read_n50(reads["length"].values)
    print(f"{len(reads)} reads, {reads['length'].sum() / 1e6:.1f} Mb, read N50 = {n50} bp")

    plt.rcParams.update({"font.size": 15, "axes.spines.top": False, "axes.spines.right": False})
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.6), gridspec_kw={"wspace": 0.34})

    # ---- a) base-weighted length histogram ---------------------------
    # Log-spaced bins; each read adds its own length (in Mb) to its bin,
    # so bar height = total bases in reads of that length.
    bins = np.logspace(np.log10(60), np.log10(80000), 55)
    ax[0].hist(reads["length"], bins=bins, weights=reads["length"] / 1e6,
               color="#0F7C8C", edgecolor="white", linewidth=0.5)
    ax[0].set_xscale("log")
    ax[0].set_xlabel("Read length (bp)")
    ax[0].set_ylabel("Sequenced bases (Mb)")
    top = ax[0].get_ylim()[1] * 1.12
    ax[0].set_ylim(0, top)
    ax[0].plot([n50, n50], [0, top * 0.87], color="#C9822B", lw=2, ls="--")
    ax[0].text(n50, top * 0.9, f"read N50 = {n50 / 1000:.1f} kb", fontsize=14, ha="center")
    ax[0].set_title("a  Where the bases are", loc="left", fontsize=15, fontweight="bold")

    # ---- b) quality vs length, hexagonal bins ------------------------
    hb = ax[1].hexbin(reads["length"], reads["quality"], xscale="log", gridsize=40,
                      bins="log", cmap="GnBu", mincnt=1, linewidths=0.2)
    ax[1].set_xlabel("Read length (bp)")
    ax[1].set_ylabel("Mean read quality (Phred)")
    colorbar = fig.colorbar(hb, ax=ax[1], pad=0.02)
    colorbar.set_label("Reads (log scale)")
    ax[1].set_title("b  Longer reads are more accurate", loc="left", fontsize=15, fontweight="bold")

    fig.savefig(args.output, dpi=220, bbox_inches="tight")
    print(f"Figure saved to {args.output}")


if __name__ == "__main__":
    main()
