import pandas as pd
import argparse
import matplotlib.pyplot as plt
import os
import numpy

def bin_age(age):
    """
    Convert numeric age into categorical bin.
    Returns a string bin label or None if age does not fit.
    """
    try:
        age = float(age)
    except ValueError:
        return None
    if age == 1:
        return "adolescent"
    elif 3 <= age <= 6:
        return "young"
    elif 9 <= age <= 15:
        return "middle"
    elif age >= 18:
        return "old"
    else:
        return None

def sort_age_groups(groups, use_binning):
    """
    Sort age groups in a chronological order.
    - If binned: adolescent -> young -> middle -> old
    - If not binned: numeric sort
    """
    if use_binning:
        order = ["adolescent", "young", "middle", "old"]
        return sorted(groups, key=lambda g: order.index(g) if g in order else len(order))
    else:
        try:
            return sorted(groups, key=lambda g: float(g))
        except ValueError:
            return sorted(groups)  # fallback to string sort

def main():
    parser = argparse.ArgumentParser(description="Compute transcript proportion variability by age/age bins and plot stacked proportions.")
    parser.add_argument("--counts", required=True, help="Transcript normalized counts CSV (rows=transcripts, cols=samples)")
    parser.add_argument("--gene_list", required=True, help="File with list of gene IDs (one per line)")
    parser.add_argument("--mapping", required=True, help="Gene-transcript mapping TSV (col1=transcript_id, col3=gene_id)")
    parser.add_argument("--metadata", required=True, help="Metadata CSV (col1=sample_id, col5=age, col7=sex)")
    parser.add_argument("--sex", choices=["m", "f", "both"], default="both", help="Filter by sex: m, f, or both")
    parser.add_argument("--tissue", default=None, help="Tissue to use")
    parser.add_argument("--output", required=True, help="Output CSV for results")
    parser.add_argument("--plot_top", type=int, default=10, help="Number of top variable genes to plot")
    parser.add_argument("--plot_dir", default="plots", help="Directory to save plots")
    parser.add_argument("--use_binning", action="store_true", help="Whether to bin ages into categories")
    parser.add_argument("--bins_to_use", default=None, help="Comma-separated list of bins/ages to keep (e.g. young,old or 1,3,6)")
    args = parser.parse_args()

    print("📂 Reading input files...")
    gene_list = pd.read_csv(args.gene_list, header=None)[0].tolist()
    mapping_df = pd.read_csv(args.mapping, sep="\t", header=None, usecols=[0, 2], names=["transcript_id", "gene_id"])
    if not gene_list[0].startswith("ENS"):
        gene_name_mapping_df = pd.read_csv(args.mapping, sep="\t", header=None, usecols=[2, 3], names=["gene_id", "gene_name"])
        genes_name_map = {}
        for i, g in gene_name_mapping_df.iterrows():
            genes_name_map[g['gene_name']] = g['gene_id']
        gene_list = list(filter(lambda x: x is not None, map(lambda x: genes_name_map[x] if x in genes_name_map else None, gene_list)))
    elif "." not in gene_list[0]:
        genes_id_map = {}
        for g in mapping_df['gene_id']:
            genes_id_map[g.split('.')[0]] = g
        gene_list = list(filter(lambda x: x is not None, map(lambda x: genes_id_map[x] if x in genes_id_map else None, gene_list)))
        
    if not gene_list:
        print("Gene list is empty, will not plot anything")
        exit(0)
    print("Gene list contains %d genes, first of which are %s" % (len(gene_list), str(gene_list[:5])))   
   
    metadata_df = pd.read_csv(args.metadata)
    counts_df = pd.read_csv(args.counts, index_col=0)
    
    # Load mapping
    mapping_full = pd.read_csv(args.mapping, sep="\t", header=None)
    mapping_full.columns = ["transcript_id", "transcript_name", "gene_id", "gene_name"]
    gene_name_map = dict(zip(mapping_full["gene_id"], mapping_full["gene_name"]))
    transcript_name_map = dict(zip(mapping_full["transcript_id"], mapping_full["transcript_name"]))

    # Identify columns for sample_id, age, sex
    sample_col = metadata_df.columns[0]
    age_col = metadata_df.columns[4]
    sex_col = metadata_df.columns[6]
    tissue_col = metadata_df.columns[2]

    # Ensure sample IDs are strings
    metadata_df[sample_col] = metadata_df[sample_col].astype(str)
    counts_df.columns = counts_df.columns.astype(str)

    print(f"ℹ Using binning: {args.use_binning}")
    if args.use_binning:
        metadata_df["age_group"] = metadata_df[age_col].apply(bin_age)
    else:
        metadata_df["age_group"] = metadata_df[age_col].astype(str)

    # Drop rows with no age group assigned
    metadata_df = metadata_df[metadata_df["age_group"].notna()]

    # Filter by sex
    if args.sex != "both":
        print(f"ℹ Filtering samples for sex: {args.sex}")
        metadata_df = metadata_df[metadata_df[sex_col].str.lower() == args.sex.lower()]
        
    if args.tissue is not None:
        print(f"ℹ Filtering samples for tissue: {args.tissue}")
        metadata_df = metadata_df[metadata_df[tissue_col].str.split('_').str[0] == args.tissue]
        
    # Filter to specific bins/ages if provided
    if args.bins_to_use:
        bins_keep = set(args.bins_to_use.split(","))
        print(f"ℹ Keeping only age groups: {', '.join(bins_keep)}")
        metadata_df = metadata_df[metadata_df["age_group"].isin(bins_keep)]

    # Filter mapping to relevant genes
    mapping_df = mapping_df[mapping_df["gene_id"].isin(gene_list)]

    results = []
    gene_avg_props = {}

    print("🔄 Processing genes...")
    for gene in gene_list:
        # Get transcripts for this gene that are present in the counts table
        transcripts = mapping_df.loc[mapping_df["gene_id"] == gene, "transcript_id"].tolist()
        transcripts = [t for t in transcripts if t in counts_df.index]
        if len(transcripts) < 2:
            continue  # Skip genes with <2 transcripts

        # Extract counts for these transcripts
        gene_counts = counts_df.loc[transcripts]

        # Merge counts with metadata to attach age group info
        merged = gene_counts.T.merge(metadata_df[[sample_col, "age_group"]], left_index=True, right_on=sample_col)

        # Skip if only one age group is present
        if merged["age_group"].nunique() < 2:
            continue

        # Compute proportions within each sample (per gene)
        merged[transcripts] = merged[transcripts].div(merged[transcripts].sum(axis=1), axis=0)

        # Convert proportions to percentages and average per age group
        avg_props = merged.groupby("age_group")[transcripts].mean() * 100

        # Sort rows by chronological order
        age_order = sort_age_groups(avg_props.index.tolist(), args.use_binning)
        avg_props = avg_props.loc[age_order]

        # Store for plotting later
        gene_avg_props[gene] = avg_props

        # Compute variability score = max transcript shift between age groups
        max_diff = 0
        max_transcript = None
        for t in transcripts:
            diff = avg_props[t].max() - avg_props[t].min()
            if diff > max_diff:
                max_diff = diff
                max_transcript = t
                
        most_expressed_transcripts = []
        for t in transcripts:
            expr = numpy.mean(avg_props[t])
            most_expressed_transcripts.append((t, expr))
        most_expressed_transcripts = list(sorted(most_expressed_transcripts, key=lambda x: x[1], reverse=True))
            
        results.append({
            "gene_id": gene,
            "gene_name": gene_name_map.get(gene, ""),
            "max_proportion_difference_percent": max_diff,
            "top_transcript": most_expressed_transcripts[0][0],
            "second_transcript": most_expressed_transcripts[1][0]
        })

    # Save results
    results_df = pd.DataFrame(results).sort_values("max_proportion_difference_percent", ascending=False)
    results_df.to_csv(args.output, index=False)
    print(f"✅ Analysis complete. Results saved to {args.output}")

    # Plotting
    print(f"📊 Generating plots...")
    os.makedirs(args.plot_dir, exist_ok=True)
    #top_genes = results_df[results_df['max_proportion_difference_percent'] >= 40]["gene_id"].tolist()
    top_genes = [] #results_df["gene_id"].tolist()
    
    for gene in top_genes:
        avg_props = gene_avg_props[gene]
        avg_props.plot(kind="bar", stacked=True, figsize=(8, 5))
        gene_name = gene_name_map.get(gene, '')
        var_value = str(int(results_df.loc[results_df["gene_id"] == gene, "max_proportion_difference_percent"].values[0]))
        gene_label = f"{gene_name} ({gene})"
        plt.title(f"{gene_label} - {var_value}")
        plt.ylabel("Proportion (%)")
        plt.xlabel("Age group")
        plt.ylim(0, 100)
        plt.legend(title="Transcript ID", bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()

        plt.savefig(os.path.join(args.plot_dir, f"{var_value}_{gene_name}_proportions.png"))
        plt.close()

    print(f"📂 Plots saved to {args.plot_dir}")

if __name__ == "__main__":
    main()

