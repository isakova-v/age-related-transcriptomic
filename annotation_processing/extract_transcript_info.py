#!/usr/bin/env python3

# === PROMPT ===
# using gffutils db make a script that for every transcript extracts the following information:
# gene id, gene name, 5' and 3' UTR lengths if present, transcript type, CDS length if present, exon count, total exonic length


import gffutils
import pandas as pd
import argparse

def transcript_features(db_path, output_tsv):
    db = gffutils.FeatureDB(db_path, keep_order=True)
    rows = []

    for transcript in db.features_of_type("transcript"):
        # Gene info
        gene = list(db.parents(transcript, featuretype="gene"))
        if gene:
            gene_id = gene[0].attributes.get("gene_id", [""])[0]
            gene_name = gene[0].attributes.get("gene_name", [""])[0]
        else:
            gene_id = ""
            gene_name = ""

        # Transcript type
        transcript_type = transcript.attributes.get("transcript_type", [""])[0] if "transcript_type" in transcript.attributes else ""

        # Exons
        exons = list(db.children(transcript, featuretype="exon", order_by="start"))
        exon_count = len(exons)
        total_exonic_length = sum(e.end - e.start + 1 for e in exons)

        # 5' UTR length
        found_exon = False
        five_utr = []
        three_utr = []
        for c in db.children(transcript, featuretype=('exon', 'UTR'), order_by="start"):
            if c.featuretype == 'exon': 
                found_exon = True
                continue
            if found_exon:
                three_utr.append(c)
            else:
                five_utr.append(c)
        if transcript.strand == '-':
            five_utr, three_utr = three_utr, five_utr
        five_utr_length = sum(u.end - u.start + 1 for u in five_utr) if five_utr else 0
        three_utr_length = sum(u.end - u.start + 1 for u in three_utr) if three_utr else 0

        # CDS length
        cds = list(db.children(transcript, featuretype="CDS", order_by="start"))
        cds_length = sum(c.end - c.start + 1 for c in cds) if cds else 0

        rows.append((
            transcript.attributes.get("transcript_id", [""])[0],
            gene_id,
            gene_name,
            five_utr_length,
            three_utr_length,
            transcript_type,
            cds_length,
            exon_count,
            total_exonic_length
        ))

    # Save table
    df = pd.DataFrame(rows, columns=[
        "transcript_id", "gene_id", "gene_name",
        "five_prime_UTR_length", "three_prime_UTR_length",
        "transcript_type", "CDS_length",
        "exon_count", "total_exonic_length"
    ])
    df.to_csv(output_tsv, sep="\t", index=False)
    print(f"[INFO] Written: {output_tsv}")

def main():
    parser = argparse.ArgumentParser(
        description="Extract transcript-level features from a gffutils database."
    )
    parser.add_argument(
        "-d", "--db",
        required=True,
        help="Path to gffutils database (.db)"
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output TSV file"
    )
    args = parser.parse_args()

    transcript_features(args.db, args.output)

if __name__ == "__main__":
    main()

