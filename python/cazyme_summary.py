#!/usr/bin/env python3
"""
cazyme_summary.py
BIO511 genomics project - strain ZS1B6-8 (Wenyingzhuangia sp.)

Summarises dbCAN results (overview.tsv) for one or more genomes:
  1. counts CAZymes per class (GH, GT, PL, CE, AA, CBM)
  2. counts CAZymes per family (e.g. GH13, GH2, PL6)
  3. checks the three target enzymes of the project:
       alpha-glucosidase           (EC 3.2.1.20)
       beta-galactosidase          (EC 3.2.1.23)
       N-acetyl-beta-glucosaminidase (EC 3.2.1.52)
  4. saves the tables as CSV files and makes two figures

A protein is counted as a CAZyme only if at least --min-tools of the three
dbCAN methods (HMMER, dbCAN-sub, DIAMOND) agree (default 2, which is the
dbCAN authors' recommendation for reliable predictions).

Example:
  python cazyme_summary.py \
      --overview annotation/dbcan/overview.tsv \
                 annotation/dbcan_refs/W_marina_type/overview.tsv \
      --names ZS1B6-8 W_marina_type \
      --outdir results/cazymes
"""

import argparse
import os
import re

import pandas as pd
import matplotlib
matplotlib.use("Agg")            # draw plots to files, no screen needed
import matplotlib.pyplot as plt
import seaborn as sns


# The six CAZy classes. Anything else dbCAN reports (e.g. SLH, cohesin)
# is ignored, because it is not a CAZyme class.
CAZY_CLASSES = ["GH", "GT", "PL", "CE", "AA", "CBM"]

# Target enzymes: EC number + the CAZy families where this activity is
# usually found (from the CAZy database). Family membership alone does not
# prove the activity, so the EC prediction from dbCAN-sub is reported too.
TARGET_ENZYMES = {
    "alpha-glucosidase":             {"ec": "3.2.1.20", "families": ["GH13", "GH31", "GH97"]},
    "beta-galactosidase":            {"ec": "3.2.1.23", "families": ["GH2", "GH35", "GH42"]},
    "N-acetyl-beta-glucosaminidase": {"ec": "3.2.1.52", "families": ["GH3", "GH20", "GH84"]},
}


def existing_file(file_path):
    """Check that a file exists (same idea as in practical 4)."""
    if not os.path.isfile(file_path):
        raise argparse.ArgumentTypeError(f"'{file_path}' is not a valid file")
    return file_path


def find_column(df, candidates):
    """Return the first column name in df that matches one of the
    candidates (ignoring case and spaces). Different dbCAN versions use
    slightly different column names, so we look for several."""
    simple = {c.lower().replace(" ", ""): c for c in df.columns}
    for cand in candidates:
        key = cand.lower().replace(" ", "")
        if key in simple:
            return simple[key]
    return None


def parse_families(text):
    """Turn a dbCAN annotation string into a list of CAZy families.

    Examples:
      "GH13_8(12-345)+CBM48(400-480)" -> ["GH13", "CBM48"]
      "GH2|CBM32"                     -> ["GH2", "CBM32"]
      "-"                             -> []
    Subfamilies (GH13_8) are merged into their family (GH13).
    """
    if not isinstance(text, str) or text.strip() in ("", "-", "N"):
        return []
    families = []
    for part in re.split(r"[+|,;]", text):
        part = re.sub(r"\(.*?\)", "", part).strip()    # remove (start-end)
        match = re.match(r"^(GH|GT|PL|CE|AA|CBM)(\d+)", part)
        if match:
            family = match.group(1) + match.group(2)
            if family not in families:                # count once per protein
                families.append(family)
    return families


def family_class(family):
    """GH13 -> GH, CBM48 -> CBM."""
    return re.match(r"^[A-Z]+", family).group(0)


def read_overview(path, min_tools):
    """Read one dbCAN overview.tsv and return a table with one row per
    (protein, family) pair, keeping only consensus CAZymes."""
    df = pd.read_csv(path, sep="\t")

    gene_col = find_column(df, ["Gene ID", "GeneID", "Gene_ID"])
    tools_col = find_column(df, ["#ofTools", "#of Tools", "NumberOfTools"])
    ec_col = find_column(df, ["EC#", "EC"])
    # Which column holds the family call: the combined recommendation
    # (dbCAN 5) if present, otherwise the HMMER result.
    fam_col = find_column(df, ["Recommend Results", "Recommend_Results", "dbCAN_hmm", "HMMER"])

    if gene_col is None or tools_col is None or fam_col is None:
        raise SystemExit(f"Unexpected columns in {path}: {list(df.columns)}")

    df = df[pd.to_numeric(df[tools_col], errors="coerce") >= min_tools]

    rows = []
    for _, row in df.iterrows():
        ec = str(row[ec_col]) if ec_col else ""
        for family in parse_families(row[fam_col]):
            rows.append({"gene": row[gene_col], "family": family,
                         "class": family_class(family), "ec": ec})
    return pd.DataFrame(rows, columns=["gene", "family", "class", "ec"])


def target_enzyme_table(hits, genome):
    """For each target enzyme: proteins in the typical families, and
    proteins with the exact EC number predicted by dbCAN-sub."""
    rows = []
    for enzyme, info in TARGET_ENZYMES.items():
        in_family = hits[hits["family"].isin(info["families"])]["gene"].nunique()
        # Match the exact EC number (so 3.2.1.2 does not match 3.2.1.20)
        ec_pattern = r"(?:^|[^\d.])" + re.escape(info["ec"]) + r"(?![\d])"
        with_ec = hits[hits["ec"].str.contains(ec_pattern, regex=True, na=False)]
        rows.append({
            "genome": genome,
            "enzyme": enzyme,
            "EC": info["ec"],
            "typical_families": "/".join(info["families"]),
            "proteins_in_typical_families": in_family,
            "proteins_with_predicted_EC": with_ec["gene"].nunique(),
            "EC_protein_ids": ";".join(sorted(with_ec["gene"].astype(str).unique())),
        })
    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser(description="Summarise dbCAN CAZyme results for one or more genomes.")
    parser.add_argument("--overview", nargs="+", type=existing_file, required=True,
                        help="One or more dbCAN overview.tsv files")
    parser.add_argument("--names", nargs="+", required=True,
                        help="A short name for each genome, in the same order as --overview")
    parser.add_argument("--min-tools", type=int, default=2,
                        help="Minimum number of dbCAN tools that must agree (default: 2)")
    parser.add_argument("--outdir", default="results/cazymes",
                        help="Folder for the output tables and figures")
    args = parser.parse_args()

    if len(args.overview) != len(args.names):
        parser.error("--overview and --names must have the same number of values")
    os.makedirs(args.outdir, exist_ok=True)

    class_counts, family_counts, targets = [], [], []

    for path, name in zip(args.overview, args.names):
        hits = read_overview(path, args.min_tools)
        n_cazymes = hits["gene"].nunique()
        print(f"{name}: {n_cazymes} CAZymes (agreed by >= {args.min_tools} tools)")

        # Proteins per class. A protein with GH + CBM counts in both classes.
        per_class = hits.groupby("class")["gene"].nunique()
        per_class = per_class.reindex(CAZY_CLASSES, fill_value=0)
        class_counts.append(per_class.rename(name))

        family_counts.append(hits.groupby("family")["gene"].nunique().rename(name))
        targets.append(target_enzyme_table(hits, name))

    # ---- Tables ------------------------------------------------------
    class_table = pd.concat(class_counts, axis=1)
    family_table = pd.concat(family_counts, axis=1).fillna(0).astype(int)
    family_table["total"] = family_table.sum(axis=1)
    family_table = family_table.sort_values("total", ascending=False)
    target_table = pd.concat(targets, ignore_index=True)

    class_table.to_csv(os.path.join(args.outdir, "cazyme_classes.csv"))
    family_table.to_csv(os.path.join(args.outdir, "cazyme_families.csv"))
    target_table.to_csv(os.path.join(args.outdir, "target_enzymes.csv"), index=False)

    print("\nCAZymes per class:")
    print(class_table)
    print("\nTarget enzymes:")
    print(target_table.drop(columns="EC_protein_ids").to_string(index=False))

    # ---- Figure 1: CAZymes per class, one bar per genome --------------
    long = class_table.reset_index().melt(id_vars="class", var_name="genome", value_name="proteins")
    plt.figure(figsize=(8, 5))
    sns.barplot(data=long, x="class", y="proteins", hue="genome", palette="Set2")
    plt.title("CAZymes per class")
    plt.xlabel("CAZy class")
    plt.ylabel("Number of proteins")
    plt.savefig(os.path.join(args.outdir, "fig_cazyme_classes.png"), dpi=300, bbox_inches="tight")
    plt.close()

    # ---- Figure 2: heatmap of the 30 most common GH and PL families ---
    # GH (glycoside hydrolases) and PL (polysaccharide lyases) are the
    # enzymes that actually break down polysaccharides.
    degrading = family_table[family_table.index.str.match(r"^(GH|PL)")]
    top = degrading.head(30).drop(columns="total")
    plt.figure(figsize=(2 + 1.2 * top.shape[1], 0.3 * len(top) + 1.5))
    sns.heatmap(top, annot=True, fmt="d", cmap="Blues", cbar_kws={"label": "Proteins"})
    plt.title("Most common GH and PL families")
    plt.xlabel("Genome")
    plt.ylabel("CAZy family")
    plt.savefig(os.path.join(args.outdir, "fig_gh_pl_heatmap.png"), dpi=300, bbox_inches="tight")
    plt.close()

    print(f"\nTables and figures saved in: {args.outdir}")


if __name__ == "__main__":
    main()
