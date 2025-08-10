#!/usr/bin/env python3
import pandas as pd
import argparse

def replace_gene_ids(mapping_file, counts_file, output_file):
    # Load mapping (ignore header, only first two columns: gene_id, gene_name)
    mapping = pd.read_csv(mapping_file, sep="\t", header=0, usecols=[0, 1])
    mapping_dict = dict(zip(mapping.iloc[:, 0], mapping.iloc[:, 1]))

    # Determine delimiter for counts file (TSV or CSV)
    sep = "\t" if counts_file.lower().endswith(".tsv") else ","

    # Load counts
    counts = pd.read_csv(counts_file, sep=sep)

    # Replace IDs in first column
    counts.iloc[:, 0] = counts.iloc[:, 0].map(lambda x: mapping_dict.get(x, x))

    # Save
    counts.to_csv(output_file, sep=sep, index=False)
    print(f"[INFO] Written: {output_file}")

def main():
    parser = argparse.ArgumentParser(
        description="Replace gene IDs with gene names in a counts table."
    )
    parser.add_argument(
        "-m", "--mapping",
        required=True,
        help="TSV file with gene_id and gene_name as first two columns"
    )
    parser.add_argument(
        "-c", "--counts",
        required=True,
        help="Counts file (TSV or CSV) with gene IDs in first column"
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output file with replaced gene names"
    )
    args = parser.parse_args()

    replace_gene_ids(args.mapping, args.counts, args.output)

if __name__ == "__main__":
    main()

