# lncRNA micropeptides — code (GitHub)

This repository holds the **analysis code** for the paper: manuscript figure scripts, supplement utilities, NetMHC merge helpers, and the smaller inputs that fit in git. Paths resolve through **`repo_paths.py`** from the repo root so clones are self-contained.

**Large files** are often copied from a private working directory or from Zenodo; a copy log may live at **`data/MANIFEST_COPIED.txt`**. (Wide NetMHC spreadsheets, merged `*_with_iedb.tsv`, IEDB dumps, and bulky SmProt/TCGA sources are usually gitignored.) File-by-file expectations: **`data/README.md`** and **`docs/ZENODO.md`**.

## What lives where

- **`manuscript/`** — scripts whose defaults write under **`figures/`** (catalog 1B, 2–4A, merged NetMHC 5–6, TTN-AS1 Fig 6 entrypoint). Supplement-only wide-XLS Fig 5A stays here for legacy/supplement use.
- **`supplement/`** — extra NetMHC panels, sensitivity grids, IEDB helpers, exploratory plots, regenerate orchestrators.
- **`pipeline/`** — SmProt / TCGA prep and filters; not the main figure entrypoints.
- **`scripts/`** — shared **`netmhc_sb_core.py`** and **`merge_netmhcpan_xls_with_iedb.py`**.
- **`notebooks/process_scratch.ipynb`** — TCGA matrices → filters → `data/primary_exp_*_lncRNA.csv` (edit paths inside if your raw downloads live elsewhere).

Top-level drivers (run from **repo root** so `import repo_paths` works):

| Script | What it does |
|--------|----------------|
| `generate_catalog_figures.py` | Fig 1B + 2–4A |
| `generate_netmhc_figure_bundle.py` | Canonical merged NetMHC Fig 5–6 |
| `rebuild_netmhc_merged_tsvs.py` | Rebuild merged TSVs from wide `*.xls` + IEDB |
| `generate_tr_lncrna_identification.py` | z-scores + R limma for Tr-lncRNA tables |
| `generate_netmhc_supplement.py` | Optional wide 5A–5E, sensitivity, Fig 6 sweeps |
| `export_publication_figures.py` | Canonical 1B–6: PNG + **PDF + TIFF** under `figures/publication/` (see folder README) |
| `supplement/regenerate_manuscript_netmhc_figures.py` | Full clean tree under `figures/manuscript_netmhc/` |

More command-level notes: `docs/figure_catalog.md`, `docs/netmhc_figure_commands.md`, `docs/iedb_tools_api.md`, `docs/figure6_ttn_as1_parameters.md`, `docs/CODE_REVIEW.md`, `data/netmhc/README_netmhc.md`.

## Setup

```bash
python -m venv .venv
.venv\Scripts\activate   # Windows
pip install -r requirements.txt
Rscript install_r_dependencies.R   # once: limma + jsonlite for tr_limma_de.R
```

`tr_limma_de.R` is at repo root for the limma step.

## Zenodo

Dataset deposit (large tables, figures bundle, and related files): **[https://doi.org/10.5281/zenodo.20167452](https://doi.org/10.5281/zenodo.20167452)**.

What ships on GitHub vs the archive, and how to refresh uploads: **`docs/ZENODO.md`**.

## License

Root **`LICENSE`** is MIT unless you replace it. Dataset rights belong in Zenodo metadata or a separate data notice if you need that — see **`docs/ZENODO.md`**.
