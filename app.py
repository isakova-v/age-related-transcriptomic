import re
from pathlib import Path

import numpy as np
import pandas as pd
import streamlit as st
import plotly.express as px
from scipy.stats import linregress
from sklearn.decomposition import PCA

# ------------------------------
# Helpers
# ------------------------------

def _std_colnames(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df.columns = (
        df.columns.str.strip()
        .str.lower()
        .str.replace(" +", "_", regex=True)
        .str.replace("[^0-9a-zA-Z_]+", "_", regex=True)
    )
    return df

def parse_age_to_months(x):
    if pd.isna(x):
        return None
    if isinstance(x, (int, float)):
        return float(x)
    m = re.search(r"([0-9]+\.?[0-9]*)\s*month", str(x).lower())
    if m:
        return float(m.group(1))
    m = re.search(r"([0-9]+\.?[0-9]*)\s*week", str(x).lower())
    if m:
        return float(m.group(1)) / 4.345
    m = re.search(r"([0-9]+\.?[0-9]*)\s*year", str(x).lower())
    if m:
        return float(m.group(1)) * 12
    return None

def load_metadata(file) -> pd.DataFrame:
    df = pd.read_csv(file)
    df = _std_colnames(df)
    rename_map = {}
    for c in df.columns:
        if c.startswith("characteristics_age") or c.startswith("characteristics__age") or c.endswith("_age"):
            rename_map[c] = "age_raw"
        elif ("sex" in c) or c.startswith("characteristics__sex") or c.startswith("characteristics_sex"):
            rename_map[c] = "sex"
        elif c == "source_name":
            rename_map[c] = "tissue"
    if "sample_name" not in df.columns:
        for alt in ["sample", "sampleid", "sra_sample", "biosample", "run", "raw_file"]:
            if alt in df.columns:
                rename_map[alt] = "sample_name"
                break
    df = df.rename(columns=rename_map)

    if "tissue" not in df.columns and "organ" in df.columns:
        df = df.rename(columns={"organ": "tissue"})

    if "tissue" in df.columns:
        df["tissue"] = df["tissue"].astype(str).str.split("_").str[0]

    if "sex" not in df.columns:
        cand = [c for c in df.columns if "sex" in c]
        if cand:
            df = df.rename(columns={cand[0]: "sex"})
        else:
            df["sex"] = pd.Series([np.nan] * len(df))

    if "age_raw" not in df.columns:
        for c in df.columns:
            if c.startswith("characteristics_age") or c.endswith("_age"):
                df = df.rename(columns={c: "age_raw"})
                break
    df["age_months"] = df.get("age_raw", pd.Series([np.nan]*len(df))).map(parse_age_to_months)
    if "tissue" not in df.columns:
        if "source_name" in df.columns:
            df["tissue"] = df["source_name"].astype(str)
        else:
            df["tissue"] = "unknown"
    if "sex" in df.columns:
        df["sex"] = df["sex"].replace({"female":"f","male":"m","F":"f","M":"m"}).astype(str)
    return df

def load_counts(file) -> pd.DataFrame:
    """
    Считываем counts без изменения имён колонок (чтобы совпали с sample_name в метаданных),
    и убираем суффиксы после первой точки: A10_...S10.gencode.vM19 -> A10_...S10
    """
    if isinstance(file, (str, Path)):
        with open(file, "r", encoding="utf-8", errors="ignore") as f:
            probe = f.read(1024)
    else:
        probe = file.read(1024); file.seek(0)
    delim = "\t" if probe.count("\t") > probe.count(",") else ","
    df = pd.read_csv(file, sep=delim)

    if df.columns[0] != "gene":
        df = df.rename(columns={df.columns[0]: "gene"})
    new_cols = []
    for i, c in enumerate(df.columns):
        if i == 0:
            new_cols.append("gene")
        else:
            base = str(c).split(".", 1)[0]
            new_cols.append(base)
    df.columns = new_cols
    dup_mask = pd.Series(df.columns).duplicated()
    if dup_mask.any():
        df = df.groupby(axis=1, level=0).sum()
        if "gene" in df.columns:
            df = df.set_index("gene")
    else:
        df = df.set_index("gene")
    return df


def cpm_normalize(counts: pd.DataFrame) -> pd.DataFrame:
    counts = counts.clip(lower=0)
    lib_size = counts.sum(axis=0)
    lib_size = lib_size.replace(0, np.nan)
    return counts.div(lib_size, axis=1) * 1e6

# ------------------------------
# UI
# ------------------------------
st.set_page_config(page_title="Mouse Aging Transcriptomics", layout="wide")

st.title("🐭 Mouse Aging Transcriptomics Dashboard")

# Автоматическая загрузка файлов при старте
metadata = load_metadata("../metadata.csv")
counts = load_counts("../counts/old_counts.csv")

norm = st.sidebar.selectbox("Normalization", ["raw", "CPM", "log1p(CPM)"])

# Harmonize samples
common_samples = [s for s in metadata["sample_name"].astype(str) if s in counts.columns]
metadata = metadata.set_index("sample_name").loc[common_samples].reset_index()
counts = counts[common_samples]
if len(common_samples) == 0:
    st.error("Нет пересечения sample_name из метаданных с колонками counts. Проверь, что имена совпадают (без суффиксов после точки).")
    st.stop()

if norm == "CPM":
    E = cpm_normalize(counts)
elif norm == "log1p(CPM)":
    E = np.log1p(cpm_normalize(counts))
else:
    E = counts

# Filters
st.sidebar.header("Filters")
tissues = sorted(metadata["tissue"].astype(str).unique())
tissue_sel = st.sidebar.multiselect("Tissue", tissues, default=tissues)
sex_col_exists = ("sex" in metadata.columns)
sexes = sorted(metadata["sex"].dropna().astype(str).unique()) if sex_col_exists else []
sex_sel = st.sidebar.multiselect("Sex", sexes, default=sexes) if sex_col_exists else []

min_age, max_age = (
    float(np.nanmin(metadata["age_months"])) if metadata["age_months"].notna().any() else 0.0,
    float(np.nanmax(metadata["age_months"])) if metadata["age_months"].notna().any() else 36.0,
)
age_rng = st.sidebar.slider("Age (months)", min_value=0.0, max_value=max(1.0, max_age), value=(min_age, max_age), step=0.5)

mask = (
    metadata["tissue"].astype(str).isin(tissue_sel)
    & (metadata["sex"].astype(str).isin(sex_sel) if sex_col_exists and len(sex_sel) > 0 else True)
 )

if metadata["age_months"].notna().any():
    mask &= metadata["age_months"].between(age_rng[0], age_rng[1])

meta_f = metadata.loc[mask].copy()
E_f = E[meta_f["sample_name"].astype(str)]

# ------------------------------
# Tabs
# ------------------------------

tab_overview, tab_gene, tab_pca, tab_heatmap, tab_reg = st.tabs([
    "Overview",
    "Gene trends",
    "PCA",
    "Heatmap",
    "Age regression",
])

# --- Overview ---
with tab_overview:
    c1, c2 = st.columns([1,1])
    with c1:
        st.subheader("Samples distribution")
        color_col = "sex" if "sex" in meta_f.columns else None
        fig = px.histogram(meta_f, x="tissue", color=color_col, barmode="group")
        st.plotly_chart(fig, use_container_width=True)
    with c2:
        st.subheader("Age distribution (months)")
        if meta_f["age_months"].notna().any():
            fig = px.histogram(meta_f, x="age_months", nbins=20, color="tissue")
            st.plotly_chart(fig, use_container_width=True)
        else:
            st.info("No numeric age parsed in metadata.")
    st.subheader("Metadata preview")
    st.dataframe(meta_f.head(50), use_container_width=True)

# --- Gene trends ---
with tab_gene:
    st.subheader("Expression vs age")
    gene = st.text_input("Gene ID/symbol (must match counts index)")
    group_by_tissue = st.checkbox("Color by tissue", value=True)
    if gene:
        if gene in E_f.index:
            df_long = (
                E_f.loc[[gene]].T
                .merge(meta_f.set_index("sample_name"), left_index=True, right_index=True)
                .reset_index(names="sample")
                .rename(columns={gene: "expr"})
            )
            if df_long["age_months"].notna().any():
                color = "tissue" if group_by_tissue else "sex" if "sex" in df_long.columns else None
                fig = px.scatter(df_long, x="age_months", y="expr", color=color, hover_data=["sample", "tissue", "sex"])
                st.plotly_chart(fig, use_container_width=True)
            else:
                st.info("Age (months) not available to plot trends.")
        else:
            st.warning(f"Gene '{gene}' not found in counts index.")
    st.caption(f"Normalization: {norm}")

# --- PCA ---
with tab_pca:
    st.subheader("PCA of samples")
    n_comp = st.slider("Components", 2, 10, 2)
    # simple variance filtering for speed
    var = E_f.var(axis=1)
    topk = st.slider("Top variable genes", 100, int(max(200, min(5000, len(var)))), 1000)
    top_genes = var.sort_values(ascending=False).head(topk).index
    X = E_f.loc[top_genes].T.fillna(0.0)
    pca = PCA(n_components=n_comp, random_state=0)
    Xp = pca.fit_transform(X)
    dfp = meta_f.copy()
    dfp["PC1"], dfp["PC2"] = Xp[:, 0], Xp[:, 1]
    fig = px.scatter(dfp, x="PC1", y="PC2", color="tissue", symbol=("sex" if "sex" in dfp.columns else None), hover_data=["sample_name", "age_months"])
    st.plotly_chart(fig, use_container_width=True)
    st.caption(f"Explained variance PC1={pca.explained_variance_ratio_[0]:.2%}, PC2={pca.explained_variance_ratio_[1]:.2%}")

# --- Heatmap ---
with tab_heatmap:
    st.subheader("Heatmap: genes across samples")
    # choose genes by variance or list
    mode = st.radio("Select genes", ["Top variable", "By list"], horizontal=True)
    if mode == "Top variable":
        k = st.slider("Number of genes", 20, 200, 50)
        genes = E_f.var(axis=1).sort_values(ascending=False).head(k).index
    else:
        txt = st.text_area("Paste gene IDs (one per line)")
        genes = [g for g in map(str.strip, txt.splitlines()) if g]
    if len(genes) > 0:
        dfH = E_f.loc[E_f.index.intersection(genes)]
        # z-score per gene for nicer heatmap
        Z = dfH.subtract(dfH.mean(axis=1), axis=0).div(dfH.std(axis=1).replace(0, np.nan), axis=0)
        fig = px.imshow(Z, aspect="auto", color_continuous_scale="RdBu_r", origin="lower",
                        labels=dict(x="Samples", y="Genes", color="z-score"))
        fig.update_yaxes(tickmode='array', tickvals=list(range(len(Z.index))), ticktext=Z.index.tolist())
        st.plotly_chart(fig, use_container_width=True)
        st.caption("Values are z-scored per gene.")

# --- Age regression ---
with tab_reg:
    st.subheader("Per-gene linear regression: expression ~ age (by tissue)")
    if meta_f["age_months"].notna().any():
        tissue_for_reg = st.selectbox("Tissue", sorted(meta_f["tissue"].unique()))
        min_n = st.slider("Min samples per gene", 5, 50, 8)
        # subset
        idx = meta_f[meta_f["tissue"] == tissue_for_reg].index
        samples = meta_f.loc[idx, "sample_name"].tolist()
        Ef = E[samples]
        ages = meta_f.loc[idx, "age_months"].values
        # run regression (vectorized-ish)
        results = []
        for gene, y in Ef.iterrows():
            yv = y.values.astype(float)
            ok = ~np.isnan(yv) & ~np.isnan(ages)
            if ok.sum() >= min_n:
                lr = linregress(ages[ok], yv[ok])
                results.append((gene, lr.slope, lr.pvalue, lr.rvalue**2, ok.sum()))
        if results:
            out = pd.DataFrame(results, columns=["gene", "slope", "pvalue", "r2", "n"]).sort_values("pvalue")
            st.dataframe(out.head(200), use_container_width=True)
            fig = px.scatter(out, x="slope", y=-np.log10(out["pvalue"].replace(0, np.nextafter(0,1))), hover_name="gene")
            fig.update_layout(xaxis_title="Slope per month", yaxis_title="-log10(p-value)")
            st.plotly_chart(fig, use_container_width=True)
        else:
            st.info("Not enough samples to run regression with current filters.")
    else:
        st.info("Age (months) not available in metadata.")

st.markdown("---")
st.caption("Tip: counts should be genes×samples, metadata must include 'sample_name'. This app computes CPM for quick exploration; for publication-grade results consider TPM/DESeq2/edgeR pipelines.")
