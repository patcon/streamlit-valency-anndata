import streamlit as st
import tempfile
import os

import valency_anndata as val
from anndata import AnnData


# ----------------------------
# Utilities
# ----------------------------

def find_bad_statement_columns(df, adata: AnnData):
    """
    Detect columns in uns["statements"] that fail h5ad serialization.
    """
    bad = []

    for col in df.columns:
        ad = adata.copy()
        ad.uns = {}

        test_df = df[[col]].copy()
        test_df.index = test_df.index.astype(str)
        test_df.columns = test_df.columns.astype(str)

        ad.uns["statements"] = test_df.to_dict(orient="split")

        try:
            with tempfile.NamedTemporaryFile(suffix=".h5ad", delete=True) as f:
                ad.write_h5ad(f.name)
        except Exception:
            bad.append(col)

    return bad


def sanitize_for_export(adata: AnnData):
    """
    Make AnnData safe for downstream JS / h5ad consumers.
    """
    adata = adata.copy()

    # --- fix uns["statements"] ---
    if "statements" in adata.uns:
        df = adata.uns["statements"]
        bad_cols = find_bad_statement_columns(df, adata)

        for col in bad_cols:
            adata.uns["statements"][col] = (
                adata.uns["statements"][col]
                .map(lambda x: "" if x is None else str(x))
            )

    # --- force cluster labels to strings ---
    categorical_cols = [
        c for c in adata.obs.columns
        if c.startswith("kmeans_")
    ]

    for col in categorical_cols:
        adata.obs[col] = (
            adata.obs[col]
            .astype(object)
            .infer_objects(copy=False)
            .fillna(-2)
            .astype(str)
        )

    return adata


# ----------------------------
# Streamlit UI
# ----------------------------

st.set_page_config(page_title="Polis → h5ad exporter", layout="centered")

st.title("🗳️ Polis → AnnData (.h5ad)")
st.caption("Run Polis-style processing and download a clean h5ad file")

polis_url = st.text_input(
    "Polis report URL",
    value="https://pol.is/report/r29kkytnipymd3exbynkd",
)

st.subheader("Projections")

run_pacmap = st.checkbox("Include PaCMAP + kmeans", value=True)
run_localmap = st.checkbox("Include LocalMAP + kmeans", value=True)

st.subheader("Metrics")

run_qc = st.checkbox(
    "Calculate QC metrics (calculate_qc_metrics)",
    value=False,
)

st.divider()

if st.button("Run pipeline and export", type="primary"):
    if not polis_url.strip():
        st.error("Please provide a Polis report URL.")
        st.stop()

    with st.spinner("Loading Polis data…"):
        adata = val.datasets.polis.load(
            polis_url,
            translate_to="en",
        )

    with st.spinner("Running Polis PCA + kmeans…"):
        val.tools.recipe_polis(adata)

    if run_qc:
        with st.spinner("Calculating QC metrics…"):
            val.preprocessing.calculate_qc_metrics(
                adata,
                inplace=True,
            )

    if run_pacmap or run_localmap:
        with st.spinner("Running additional projections…"):
            if run_pacmap:
                val.tools.pacmap(
                    adata,
                    layer="X_masked_imputed_mean",
                )
                val.tools.kmeans(
                    adata,
                    mask_obs="cluster_mask",
                    init="random",
                    use_rep="X_pacmap",
                    key_added="kmeans_pacmap",
                )

            if run_localmap:
                val.tools.localmap(
                    adata,
                    layer="X_masked_imputed_mean",
                )
                val.tools.kmeans(
                    adata,
                    mask_obs="cluster_mask",
                    init="random",
                    use_rep="X_localmap",
                    key_added="kmeans_localmap",
                )

    with st.spinner("Sanitizing AnnData for export…"):
        adata_export = sanitize_for_export(adata)

    with tempfile.NamedTemporaryFile(
        suffix=".h5ad",
        delete=False,
    ) as f:
        export_path = f.name
        adata_export.write_h5ad(export_path)

    with open(export_path, "rb") as f:
        st.success("Done!")
        st.download_button(
            label="⬇️ Download h5ad file",
            data=f,
            file_name="polis_export.h5ad",
            mime="application/octet-stream",
        )

    st.markdown(
        """
After downloading, you can **[visit the app](https://y6p9xj.csb.app/)** to visualize the export.
"""
    )

    st.markdown(
        """
For an alternative, **advanced** workflow, the export can also be generated directly via a Python notebook  
(**[Advanced](https://colab.research.google.com/drive/1h-bp6FsiCFHH0qH8s34jOSZY5-BuGPkH)**).
"""
    )

    os.remove(export_path)

st.markdown(
    """
<hr style="margin-top: 3rem; margin-bottom: 1rem;" />
<sub><a href="https://github.com/patcon/streamlit-valency-anndata">Code</a></sub>
""",
    unsafe_allow_html=True,
)
