import streamlit as st
import tempfile
import os

import valency_anndata as val


# ----------------------------
# Streamlit UI
# ----------------------------

st.set_page_config(page_title="Polis → h5ad exporter", layout="centered")

st.title("🗳️ Polis → AnnData (.h5ad)")
st.caption("Run Polis-style processing and download a clean h5ad file — powered by [valency-anndata](https://patcon.github.io/valency-anndata/)")

default_url = st.query_params.get("report", "https://pol.is/report/r29kkytnipymd3exbynkd")
default_lang = st.query_params.get("lang", "")

polis_url = st.text_input(
    "Polis report URL",
    value=default_url,
)

translate_to = st.text_input(
    "Translate to language (2-letter code, e.g. en)",
    value=default_lang,
    max_chars=2,
)

st.subheader("Projections")

run_pacmap = st.checkbox("Include PaCMAP + kmeans", value=True)
run_localmap = st.checkbox("Include LocalMAP + kmeans", value=True)

st.subheader("Metrics")

run_qc = st.checkbox(
    "Calculate QC metrics (calculate_qc_metrics)",
    value=True,
)

st.divider()

if st.button("Run pipeline and export", type="primary"):
    if not polis_url.strip():
        st.error("Please provide a Polis report URL.")
        st.stop()

    with st.spinner("Loading Polis data…"):
        load_kwargs = dict()
        if translate_to.strip():
            load_kwargs["translate_to"] = translate_to.strip()
        adata = val.datasets.polis.load(
            polis_url,
            **load_kwargs,
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

    with tempfile.NamedTemporaryFile(
        suffix=".h5ad",
        delete=False,
    ) as f:
        export_path = f.name

    with st.spinner("Exporting h5ad…"):
        val.write(export_path, adata)

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
<br/>
<sub>Tip: You can pre-fill the form via query string, e.g. <code>?report=URL&amp;lang=en</code></sub>
""",
    unsafe_allow_html=True,
)
