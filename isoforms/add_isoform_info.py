#!/usr/bin/env python3
import argparse
import csv

def main():
    parser = argparse.ArgumentParser(description="Merge transcript info from two files.")
    parser.add_argument(
        "-t", "--transcript-info", required=True,
        help="TSV file with transcript info (transcript ID in col1, data in cols 4-10)"
    )
    parser.add_argument(
        "-v", "--variable-transcripts", required=True,
        help="CSV file with variable transcripts (most abundant transcripts in cols 4 and 5)"
    )
    parser.add_argument(
        "-o", "--output", required=True,
        help="Output TSV file"
    )
    args = parser.parse_args()

    # Load transcript info from the first file
    transcript_data = {}
    with open(args.transcript_info, newline='') as tsvfile:
        reader = csv.reader(tsvfile, delimiter="\t")
        header_tsv = next(reader)
        colnames_4_10 = header_tsv[3:10]  # actual names for columns 4-10
        for row in reader:
            transcript_id = row[0]
            if len(row) >= 10:
                transcript_data[transcript_id] = row[3:10]  # store cols 4–10

    # Process variable transcripts file
    output_rows = []
    with open(args.variable_transcripts, newline='') as csvfile:
        reader = csv.reader(csvfile)
        header_csv = next(reader)

        # Create new header: original CSV header + transcript info headers for col4 + for col5
        new_header = (
            header_csv +
            [f"{name}_transcript4" for name in colnames_4_10] +
            [f"{name}_transcript5" for name in colnames_4_10]
        )
        output_rows.append(new_header)

        for row in reader:
            transcript1_id = row[3]
            transcript2_id = row[4]

            transcript1_info = transcript_data.get(transcript1_id, ["NA"] * len(colnames_4_10))
            transcript2_info = transcript_data.get(transcript2_id, ["NA"] * len(colnames_4_10))

            output_rows.append(row + transcript1_info + transcript2_info)

    # Write output TSV
    with open(args.output, "w", newline='') as outfile:
        writer = csv.writer(outfile, delimiter="\t")
        writer.writerows(output_rows)


if __name__ == "__main__":
    main()

