# Mouse Aging Transcriptomics Dashboard

A Streamlit app for rapid, exploratory analysis of age-related transcriptomic patterns in mouse tissues. It helps you:

- Load counts and metadata, harmonize IDs, and choose normalization.
- Visualize gene trends vs. age, PCA/UMAP embeddings, and detect sample outliers.
- Discover single-peak age clusters and annotate them with GO terms (per tissue).
- Run per-gene age regressions.
- Explore GO term presence and enrichment across tissues.
- Cluster age trajectories (KMeans) and perform GO enrichment on clusters.

> Note: an “Isoforms (splicing)” tab is scaffolded but currently commented out.

---

## 0) Reproducibility micro-benchmark

![Reproducibility](https://github.com/isakova-v/age-related-transcriptomic/tree/main/images/image.jpg](https://github.com/isakova-v/age-related-transcriptomic/blob/main/images/salmon_vs_HTSeq.png)

> Motivation:
there are multiple approaches and pipelines for analying RNA-seq data. It may be interesting to compare some different popular approaches and 
check how distinguish are they.

> Realisation:
we used STAR+HTSeq from the original article (alignment and unique count estimation) and salmon (which produces pseudocounts)

> Results:
methods' convergence is far from ideal but isn't catastrophic. At the level of the raw counts the pearson correlation between methods was in range 0.7-0.9.

For more information please visit the [page](https://github.com/isakova-v/age-related-transcriptomic/tree/main/htseq_salmon_fc_bench) of reproducibility readme.

## 1) Quick start

### Prerequisites
- Python ≥ 3.9 (3.10+ recommended)
- Command line access to run Streamlit
- Internet access for fetching GO libraries via **gseapy** (Enrichr)

### Installation
```bash
# create & activate a virtualenv (example)
python -m venv .venv
source .venv/bin/activate  # on Windows: .venv\Scripts\activate

# install dependencies
pip install -U pip
pip install streamlit pandas numpy plotly scikit-learn umap-learn gseapy scipy pyyaml
```

## Minimal project layout
```bash
project/
├─ app.py                 # this application
├─ config.yaml            # your configuration (see below)
├─ data/
│  ├─ metadata.csv        # sample metadata
│  ├─ counts_matrix.csv   # or .tsv — genes × samples
│  └─ id_mapping.csv      # optional, for gene ID harmonization
└─ go/                    # optional local GO assets if you keep any
```

### Run
```bash
streamlit run app.py
```
The app reads settings from config.yaml. You can keep your data anywhere; just point to them in the config.

## 2) Configuration (config.yaml)
Below is a comprehensive example. Omit or tweak sections you don’t need.
```yaml
# Page/UI
ui:
  page_title: "Mouse Aging Transcriptomics"

ui_defaults:
  filters:
    tissue_default: "__ALL__"      # "__ALL__" or list: ["liver","brain"]
    sex_default: "__ALL__"         # "__ALL__" or list: ["m","f"]
    age_months_default: [0, 36]    # optional slider default

# Data paths
paths:
  metadata: "data/metadata.csv"    # must include 'sample_name'
  datasets:                        # name -> counts path (genes × samples)
    Main: "data/counts_matrix.csv"
    Alt : "data/counts_matrix_alt.tsv"
  id_mapping: "data/id_mapping.csv"  # optional; see formats below

# Normalization & gene ID cleaning
normalization:
  mode_default: "log1p(CPM)"       # one of: raw | CPM | log1p(CPM)
  clean_gene_ids:
    strip_version: true            # "ENSMUSG... .1" → "ENSMUSG..."
    preserve_case: true            # keep case; if false → upper-case

# Embedding (PCA/UMAP)
embedding:
  top_variable_genes: 2000         # default slider upper bound caps apply
  umap:
    n_neighbors: 15
    min_dist: 0.10
  outliers:
    z_threshold: 3.0               # for outlier detection in embeddings

# Heatmap (single-peak clusters)
heatmap_single_peak:
  bin_width: 3.0                   # months
  dominance_ratio: 1.5             # max/second-max z-score ratio
  min_genes_in_cluster: 15
  top_clusters_to_show: 6

# Per-gene regression
regression:
  min_samples_per_gene: 8

# GO Enrichment (via gseapy/Enrichr)
organism:
  name: "mouse"                    # used by default if not overridden below
go_enrichment:
  default_library: "GO_Biological_Process_2023"
  organism: "mouse"                # fallback if organism.name missing
  min_genes_term_presence: 1
  fdr_line: 0.05
  libraries:                       # label -> [enrichr_library_name, organism]
    GO: ["GO_Biological_Process_2023", "mouse"]
    GO_CC: ["GO_Cellular_Component_2023", "mouse"]
    GO_MF: ["GO_Molecular_Function_2023", "mouse"]

# Patterns (trajectory clustering)
patterns_clustering:
  by_tissue_default: true
  bin_width: 3.0
  kmeans_k: 12
  min_bins_per_gene: 4
  zscore_trajectories: true
  random_state: 0
```
## 3) Input data formats

### 3.1 Metadata (`paths.metadata`)
CSV/TSV with at least a **sample identifier**. The app tries hard to standardize column names. It will create:
- `sample_name`: taken from existing columns like `sample`, `run`, `biosample`, etc.
- `tissue`: from `tissue` or `source_name`/`organ` if present (underscores stripped to the first token).
- `sex`: normalized to `{"m","f"}` when possible; otherwise left as given or NaN.
- `age_months`: parsed numeric months derived from `age_raw` or other “age-like” fields.

**Age parsing** supports:
- Raw months: `12`, `12.5`, `12 m`, `12mo`, `12 months`
- Weeks/days/years: `8 w`, `30 d`, `1 y`, `1.5 years` → converted to months
- Ranges: `3-5 months` → average
- Postnatal day: `P21` → 21 days → months

**Minimum required column**: something that can become `sample_name`.

---

### 3.2 Counts (`paths.datasets`)
A **genes × samples** table (first column are gene IDs/symbols). CSV or TSV is auto-detected.

- The first column is renamed to `gene`.
- Gene names are “soft-cleaned” according to config: strip dot-versions (`X.Y` → `X`) and optionally change case.
- All other column names (samples) will have dot-suffixes stripped if configured (to match `metadata.sample_name`).
- Duplicated columns (after cleaning) are summed; genes become the index.

**Important:** The app keeps only columns (samples) present in both metadata and counts after cleaning.

---

### 3.3 Optional ID mapping (`paths.id_mapping`)
Use when your counts’ gene IDs are not official symbols. The app builds a mapper to **official symbols (MGI)** from any of these source columns if present:

- Source candidates: `gene`, `symbol`, `alias`, `synonym`, `ensembl`, `ensembl_gene`, `uniprot`, `uniprot_id`, `uni`
- Target column (required): one of `symbol`, `gene_symbol`, `mgi_symbol`

**Example (`id_mapping.csv`):**
```csv
ensembl_gene,symbol,alias
ENSMUSG00000064370,Atp5f1a,Atp5a
ENSMUSG00000064371,Cox4i1,COXIV
Q9D...,Ndufs2,
```

## 4) Normalization

Choose in the sidebar:
- **raw**: as-is
- **CPM**: counts per million (library size)
- **log1p(CPM)**: quick variance-stabilized view

> For publication-grade pipelines consider TPM or DESeq2/edgeR outside this app.

---

## 5) UI overview (tabs)

### A) **Overview**
- Sample counts per tissue (histogram).
- Age distribution (months).
- Metadata preview.

### B) **Gene trends**
- Enter a **gene ID/symbol** (must exist in the counts index after mapping).
- See mean expression vs. age per tissue (lines with hover lists of samples/values).
- Boxplots by age points for a selected tissue.

### C) **PCA & UMAP**
- Select top-variable genes and PCA/UMAP parameters.
- Interactive 2D embedding colored by tissue, symbol by sex (if present).
- **Outlier detection**:
  - By **tissue** clusters (distance-to-centroid z-score > threshold).
  - By **age bins** (customizable bin width).
  - Lists sample names identified as outliers per group.

### D) **Heatmap**
- Bins ages and computes **z-scores** per gene across bins.
- Finds genes with a **single dominant peak** (max/second-max ≥ threshold).
- Clusters genes by peak bin; shows stacked heatmaps by peak age.
- For each cluster, reports **top-6 (GO term, tissue)** pairs using hypergeometric enrichment on tissue-specific backgrounds.

### E) **Age regression**
- For one tissue at a time, fits per-gene **linear regression**:
  `expression ~ age_months`
- Displays table (slope, p-value, R², n) and a volcano-like scatter.

### F) **GO groups**
- Loads a GO library (Enrichr via gseapy).
- Shows which **terms have genes present** in your dataset.
- Drill-down into a selected term: size in library vs. matches in your data.
- **Enrichment across tissues**: builds invariant gene sets per tissue (by CV/mean/non-zero thresholds) and tests a selected term across tissues; shows FDR bar plot and downloadable table.

### G) **Patterns (ML)**
- Cluster **gene trajectories over age** (KMeans) inside one tissue or across all.
- Bins ages, averages per bin, interpolates gaps, optional z-score per gene.
- Shows cluster means, then **details** for a selected cluster:
  - Individual trajectories, heatmap
  - Download gene list
  - **GO enrichment** for the cluster (selectable library)
  - (If clustering by tissue) plot **term trajectory** (term∩cluster) over age in the selected tissue.

> The **Isoforms** tab is provided in commented code for future work on splicing proportions by age.

---

## 6) Caching & performance

- Expensive steps are cached:
  - `get_go_lib_cached`: Enrichr library fetching via **gseapy**.
  - Computations for enrichment and clustering steps use `@st.cache_*`.
- If you change **data files** or **config**, use Streamlit’s “Rerun” and/or “Clear cache” from the app menu.

**Tips:**
- Start with a smaller dataset to validate `config.yaml` and formats.
- Use fewer top-variable genes and moderate UMAP neighbors for speed.
- Narrow filters (tissue/age) to accelerate downstream steps.

---

## 7) Common issues & troubleshooting

- **“Metadata file not found” / “File not found: …”**  
  Check `paths` in `config.yaml`.

- **“No datasets specified in config.”**  
  Add at least one entry under `paths.datasets`.

- **“No overlap between metadata sample_name and counts_std columns.”**  
  Ensure sample names in metadata match counts column names **after** cleaning (dot-suffix stripping). If your counts columns carry extra suffixes (e.g., `.gencode.vM19`), set `clean_gene_ids.strip_version: true` and confirm the metadata uses the base names.

- **GO library errors**  
  Requires internet to fetch Enrichr libraries via **gseapy**. Check `go_enrichment.default_library` and `organism`. You can switch libraries in tabs where a selector is offered.

- **Age parsing issues**  
  If ages don’t parse into months, verify the raw strings and consider pre-cleaning or providing a consistent `age` column.

- **Few samples / empty results**  
  Some analyses require minimum samples per gene or bins with multiple points. Loosen filters, widen bin width, or aggregate tissues.

---

## 8) Reproducibility notes

This app focuses on **exploration**. For formal analyses:
- Re-run pipelines using fixed seeds (UMAP, KMeans) and pinned packages.
- Export the **filters, bin widths, and library choices** you used.
- Reproduce differential analyses with established tools (e.g., DESeq2/edgeR, limma-voom) or statistical models suited to your design (mixed models for repeated measures).

---

## 9) Developer notes

- **ID mapping** is applied early to standardize gene indices to official symbols. Duplicated symbols are summed.
- **Hypergeometric tests** are implemented in-house (`scipy.stats.hypergeom`) with **BH-FDR** correction.
- **Outlier detection** uses centroid distance within groups → z-score thresholding.
- **Config-driven** UI defaults keep the app adaptable to new datasets.

---

## 10) Example minimal inputs

**metadata.csv**
```csv
sample_name,tissue,sex,age_raw
S1,liver,m,3 mo
S2,liver,f,6 months
S3,brain,m,P21
S4,brain,f,1 year
```

### counts_matrix.csv

A table of gene expression counts: **genes × samples**.
	•	First column: gene (gene IDs or symbols)
	•	Remaining columns: counts for each sample (column names must match sample_name in metadata after cleaning)

 Example:
 ```csv
gene,S1,S2,S3,S4
Atp5f1a,100,200,50,75
Ndufs2,20,15,10,9
Cox4i1,300,280,150,120
```

### id_mapping.csv (optional)

Maps IDs in the counts file to official gene symbols.
	•	At least one column with source IDs (e.g., ensembl_gene, uniprot)
	•	One column with official gene symbols (symbol, gene_symbol or mgi_symbol)

Example:
```csv
ensembl_gene,symbol
ENSMUSG00000064370,Atp5f1a
ENSMUSG00000012345,Ndufs2
ENSMUSG00000099999,Cox4i1
```

## 11) License & citation

This dashboard is built using open-source Python packages:
	•	Streamlit — interactive web app framework
	•	pandas, NumPy, SciPy — data processing and statistical analysis
	•	scikit-learn, umap-learn — dimensionality reduction and clustering
	•	plotly — interactive plotting
	•	gseapy — GO enrichment analysis (via Enrichr)

If you use results or methods from this app in a publication:
	•	Cite Enrichr and gseapy
	•	Cite the Gene Ontology resource and any specific libraries used for enrichment

 ## 12) Contact

Please open an issue in the repository or submit a pull request with your changes.

Happy exploring! 🐭
