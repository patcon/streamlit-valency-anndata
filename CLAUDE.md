# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Development

```bash
uv venv
uv pip install -r requirements.txt
uv run streamlit run app.py
```

Visit http://localhost:8501. There are no tests or linting configured.

## Architecture

This is a single-file Streamlit app (`app.py`) that wraps the [`valency-anndata`](https://patcon.github.io/valency-anndata/) library into a web UI.

**Data flow:**
1. User provides a Polis report URL (pre-populated via `?report=` query param)
2. Optional: translation language (`?lang=` query param, supports locale codes like `zh-tw`)
3. `val.datasets.polis.load()` fetches the Polis report data
4. `val.tools.recipe_polis()` runs PCA + kmeans
5. Optional: QC metrics, PaCMAP projection, LocalMAP projection (each with their own kmeans pass)
6. `val.write()` exports to `.h5ad` (HDF5-based AnnData format)
7. User downloads the file for visualization at https://y6p9xj.csb.app/

**Key `valency_anndata` API used in app.py:**
- `val.datasets.polis.load(url, translate_to=...)` — fetches and parses Polis data
- `val.tools.recipe_polis(adata)` — standard Polis PCA + kmeans pipeline
- `val.preprocessing.calculate_qc_metrics(adata, inplace=True)`
- `val.tools.pacmap(adata, layer=...)` / `val.tools.localmap(adata, layer=...)`
- `val.tools.kmeans(adata, mask_obs=..., init=..., use_rep=..., key_added=...)`
- `val.write(path, adata)`

The `valency-anndata` dependency is installed directly from GitHub (`main` branch), so changes upstream are picked up on reinstall.
