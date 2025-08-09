#!/usr/bin/env python3
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

# --------------------
# Helper functions
# --------------------

def define_age_groups(age):
    """Map age values (months) to age bins."""
    if age == 1:
        return "Adolescent (1m)"
    elif 3 <= age <= 6:
        return "Young (3-6m)"
    elif 9 <= age <= 15:
        return "Middle (9-15m)"
    elif age >= 18:
        return "Old (18m+)"
    else:
        return "Other"

def sort_age_groups(groups, binning):
    """Sort age groups chronologically."""
    if binning:
        order = ["Adolescent (1m)", "Young (3-6m)", "Middle (9-15m)", "Old (18m+)", "Other"]
        return [g for g in order if g in groups]
    else:
        # Numeric sort if not binning
        try:
            return sorted(groups, key=lambda x: float(x))
        except ValueError:
            return sorted(groups)

# --------------------
# Main script
# --------------------

def main():
    parser = argparse.ArgumentParser(description="Compute transcript proportions and variability by age.")
    parser.add_argument("--counts", required=True, help="CSV with transcript counts (columns = samples)")
    parser.add_argument("--gene_list", required=True, help="File with list of gene IDs")
    parser.add_argument("--mapping", required=True, help="TSV with transcript_id (col1), transcript_name (col2), gene_id (col3), gene_name (col4)")
    parser.add_argument("--metadata", required=True, help="CSV with sample metadata; col1=sample id, col5=age (months), col7=sex (m/f)")
    parser.add_argument("--output", required=True, help="Output CSV for variability results")
    parser.add_argument("--plot_dir", required=True, help="Directory to save plots")
    parser.add_argument("--plot_top", type=int, default=10, help="Number of top variable genes to plot")
    parser.add_argument("--use_binning", action="store_true", help="Use age binning instead of raw ages")
    parser.add_argument("--age_filter", nargs="+", help="Specific ages or bins to include")
    parser.add_argument("--sex", choices=["m", "f", "both"], default="both", help="Filter by sex")
    args = parser.parse_args()

    print("📂 Loading data...")

    # Load counts (index = transcript IDs, columns = samples)
    counts_df = pd.read_csv(args.counts, index_col=0)

    # Load gene list
    with open(args.gene_list) as f:
        gene_list = [line.strip() for line in f if line.strip()]

    # Load mapping
    mapping_df = pd.read_csv(args.mapping, sep="\t", header=None)
    mapping_df.columns = ["transcript_id", "transcript_name", "gene_id", "gene_name"]
    gene_name_map = dict(zip(mapping_df["gene_id"], mapping_df["gene_name"]))
    transcript_name_map = dict(zip(mapping_df["transcript_id"], mapping_df["transcript_name"]))

    # Load metadata
    metadata_df = pd.read_csv(args.metadata)
    metadata_df.columns = [f"col{i+1}" for i in range(metadata_df.shape[1])]
    sample_col = "col1"
    age_col = "col5"
    sex_col = "col7"

    # Apply sex filter if needed
    if args.sex != "both":
        before_count = len(metadata_df)
        metadata_df = metadata_df[metadata_df[sex_col] == args.sex]
        print(f"👥 Filtered by sex '{args.sex}': {before_count} -> {len(metadata_df)} samples")

    # Create age_group column
    if args.use_binning:
        metadata_df["age_group"] = metadata_df[age_col].apply(define_age_groups)
    else:
        metadata_df["age_group"] = metadata_df[age_col].astype(str)

    # Apply age filter if specified
    if args.age_filter:
        before_count = len(metadata_df)
        metadata_df = metadata_df[metadata_df["age_group"].isin(args.age_filter)]
        print(f"⏳ Filtered by age: {before_count} -> {len(metadata_df)} samples")

    print("🔄 Processing genes...")

    results = []
    gene_avg_props = {}

    for gene in gene_list:
        # Find transcripts for this gene
        gene_transcripts = mapping_df[mapping_df["gene_id"] == gene]["transcript_id"].tolist()
        if not gene_transcripts:
            continue

        # Subset counts for these transcripts
        sub_counts = counts_df.loc[counts_df.index.intersection(gene_transcripts)]
        if sub_counts.empty:
            continue

        # Merge counts with metadata
        merged = sub_counts.T.merge(metadata_df[[sample_col, "age_group"]],
                                    left_index=True, right_on=sample_col)

        age_groups = sort_age_groups(merged["age_group"].unique(), args.use_binning)
        prop_table = []

        for age_group in age_groups:
            # Get counts for samples in this age group
            group_counts = merged[merged["age_group"] == age_group].set_index(sample_col)
            mean_counts = group_counts[gene_transcripts].mean()
            total = mean_counts.sum()
            if total == 0:
                proportions = [0] * len(mean_counts)
            else:
                proportions = (mean_counts / total * 100).tolist()  # percentage
            prop_table.append(proportions)

        # Create DataFrame of proportions
        prop_df = pd.DataFrame(prop_table, columns=gene_transcripts, index=age_groups)

        # Calculate variability measure (max proportion difference between any two ages)
        variability = max(prop_df.max() - prop_df.min())

        # Store results
        results.append({
            "gene_id": gene,
            "gene_name": gene_name_map.get(gene, ""),
            "max_proportion_difference_percent": variability
        })
        gene_avg_props[gene] = prop_df

    # Sort by variability
    results_df = pd.DataFrame(results)
    results_df = results_df.sort_values(by="max_proportion_difference_percent", ascending=False)
    results_df.to_csv(args.output, index=False)
    print(f"💾 Results written to {args.output}")

    # Plot top N variable genes
    os.makedirs(args.plot_dir, exist_ok=True)
    for gene in results_df.head(args.plot_top)["gene_id"]:
        prop_df = gene_avg_props[gene]
        fig, ax = plt.subplots(figsize=(8, 6))

        # Transcript labels with names
        transcripts = prop_df.columns
        transcript_labels = [f"{transcript_name_map.get(tid, tid)} ({tid})" for tid in transcripts]

        prop_df.plot(kind="bar", stacked=True, ax=ax,
                     color=plt.cm.tab20.colors[:len(transcripts)])

        ax.set_ylabel("Proportion (%)")
        ax.set_xlabel("Age / Age bin")
        ax.set_ylim(0, 100)

        gene_label = f"{gene_name_map.get(gene, '')} ({gene})"
        ax.set_title(f"Isoform proportions - {gene_label}")

        ax.legend(transcript_labels, title="Transcripts", bbox_to_anchor=(1.05, 1), loc="upper left")

        # Add variability value below legend
        var_value = results_df.loc[results_df["gene_id"] == gene, "max_proportion_difference_percent"].values[0]
        ax.text(1.05, 0.5, f"Variability: {var_value:.1f}%", transform=ax.transAxes, fontsize=10)

        plt.tight_layout()
        plot_path = os.path.join(args.plot_dir, f"{gene}.png")
        plt.savefig(plot_path, dpi=300)
        plt.close()

    print(f"📊 Plots saved to {args.plot_dir}")

if __name__ == "__main__":
    main()

