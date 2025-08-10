#!/usr/bin/env python3
import argparse
import csv
import glob
import os

def large_difference(val1, val2, abs_thresh, rel_thresh=None):
    """Check if two values differ by absolute or relative criteria."""
    try:
        v1 = float(val1)
        v2 = float(val2)
    except ValueError:
        return False  # skip if not numbers
    
    abs_diff = abs(v1 - v2)
    rel_diff = abs_diff / max(v1, v2) if max(v1, v2) > 0 else 0

    if abs_diff >= abs_thresh:
        return True
    if rel_thresh is not None and rel_diff >= rel_thresh:
        return True
    return False

def main():
    parser = argparse.ArgumentParser(description="Merge and filter TSVs by transcript length differences.")
    parser.add_argument(
        "inputs", nargs="+",
        help="List of TSV files or wildcards (e.g. '*.tsv')"
    )
    parser.add_argument(
        "-o", "--output", required=True,
        help="Output TSV file"
    )
    parser.add_argument(
        "--utr-thresh", type=int, default=200,
        help="Absolute threshold for UTR length difference (default=200bp)"
    )
    parser.add_argument(
        "--len-thresh", type=int, default=300,
        help="Absolute threshold for CDS/total length difference (default=300bp)"
    )
    parser.add_argument(
        "--rel-thresh", type=float, default=None,
        help="Relative difference threshold (e.g., 0.2 for 20%)"
    )
    args = parser.parse_args()

    # Expand wildcards
    file_list = []
    for pattern in args.inputs:
        file_list.extend(glob.glob(pattern))
    file_list = sorted(set(file_list))

    if not file_list:
        raise FileNotFoundError("No TSV files found matching the input patterns.")

    output_rows = []
    header_written = False

    for filepath in file_list:
        file_base = os.path.splitext(os.path.basename(filepath))[0]

        with open(filepath, newline='') as tsvfile:
            reader = csv.reader(tsvfile, delimiter="\t")
            header = next(reader)

            if not header_written:
                new_header = header[:2] + ["source_file"] + header[2:]
                output_rows.append(new_header)
                header_written = True

            for row in reader:
                if len(row) < 19:
                    continue  # skip malformed rows

                utr5_1, utr3_1 = row[5], row[6]
                cds_1, total_1 = row[9], row[11]
                utr5_2, utr3_2 = row[12], row[13]
                cds_2, total_2 = row[16], row[18]

                keep = (
                    large_difference(utr5_1, utr5_2, args.utr_thresh, args.rel_thresh) or
                    large_difference(utr3_1, utr3_2, args.utr_thresh, args.rel_thresh) or
                    large_difference(cds_1, cds_2, args.len_thresh, args.rel_thresh) or
                    large_difference(total_1, total_2, args.len_thresh, args.rel_thresh)
                )

                if keep:
                    new_row = row[:2] + [file_base.split('.')[0]] + row[2:]
                    output_rows.append(new_row)

    # Write output
    with open(args.output, "w", newline='') as outfile:
        writer = csv.writer(outfile, delimiter="\t")
        writer.writerows(output_rows)


if __name__ == "__main__":
    main()

