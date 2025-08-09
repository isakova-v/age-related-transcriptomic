#!/usr/bin/env python3
import argparse
import csv
import sys
import os

csv.field_size_limit(sys.maxsize)

def load_id_mapping(valid_ids_file):
    """Load mapping from ID without version → full ID with version."""
    mapping = {}
    with open(valid_ids_file, 'r') as f:
        for line in f:
            full_id = line.strip().split()[0]  # works for TSV/whitespace
            base_id = full_id.split('.')[0]
            mapping[base_id] = full_id
    return mapping

def sanitize_filename(name):
    """Replace problematic filename characters."""
    return "".join(c if c.isalnum() or c in ('_', '-') else "_" for c in name)

def process_csv(input_csv, mapping, out_dir, no_header=False):
    with open(input_csv, 'r', newline='') as f:
        reader = csv.reader(f)
        if not no_header:
            next(reader)  # skip header

        for row in reader:
            if len(row) < 3:
                continue
            tissue = row[0].strip()
            inv_genes = row[1].strip().split(',')
            var_genes = row[2].strip().split(',')

            inv_corrected = [mapping[g] for g in inv_genes if g in mapping]
            var_corrected = [mapping[g] for g in var_genes if g in mapping]

            tissue_safe = sanitize_filename(tissue)

            inv_path = os.path.join(out_dir, f"{tissue_safe}_inv.tsv")
            var_path = os.path.join(out_dir, f"{tissue_safe}_var.tsv")

            with open(inv_path, 'w') as out_inv:
                out_inv.write("\n".join(inv_corrected) + "\n")

            with open(var_path, 'w') as out_var:
                out_var.write("\n".join(var_corrected) + "\n")

def main():
    parser = argparse.ArgumentParser(description="Correct Gencode IDs in CSV and create per-tissue invariant/variable files.")
    parser.add_argument("input_csv", help="CSV file with tissue, invariant genes, variable genes")
    parser.add_argument("valid_ids_tsv", help="TSV file with valid gene IDs (with version)")
    parser.add_argument("out_dir", help="Output directory for tissue-specific files")
    parser.add_argument("--no_header", action="store_true", help="Specify if CSV has no header")
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    mapping = load_id_mapping(args.valid_ids_tsv)

    process_csv(args.input_csv, mapping, args.out_dir, args.no_header)

if __name__ == "__main__":
    main()

