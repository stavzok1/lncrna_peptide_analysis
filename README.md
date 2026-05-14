# lncRNA micropeptides — code (GitHub)

This repository holds the **analysis code** for the paper: manuscript figure scripts, supplement utilities, NetMHC merge helpers, and the smaller inputs that fit in git. Paths resolve through **`repo_paths.py`** from the repo root so clones are self-contained.

**Large files** are often copied from a private working directory or from Zenodo; a copy log may live at **`data/MANIFEST_COPIED.txt`**. (Wide NetMHC spreadsheets, merged `*_with_iedb.tsv`, IEDB dumps, and bulky SmProt/TCGA sources are usually gitignored.) File-by-file expectations: **`data/README.md`** and **`docs/ZENODO.md`**.

## What lives where

- **`manuscript/`** — scripts whose defaults write canonical panels under **`figures/`** (catalog 1B, 2–4A, merged NetMHC **5A–5E** for the main coding cohort, TTN-AS1 Fig 6 entrypoint). Supplement-only wide-XLS Fig 5A stays here for legacy/supplement use.
- **`figures/supplementary/`** — non-canonical mirrors grouped by topic (PCA supplement; NetMHC random-fragment cohort for Fig 5; **Fig 5–6 NetMHC sensitivity + combination grid + TTN SB sweeps** under `netmhc_fig5_fig6_supplement/`). See **`figures/supplementary/README.md`**.
- **`supplement/`** — extra NetMHC panels, sensitivity grids, IEDB helpers, exploratory plots, regenerate orchestrators.
- **`pipeline/`** — SmProt / TCGA prep and filters; not the main figure entrypoints.
- **`scripts/`** — shared **`netmhc_sb_core.py`** and **`merge_netmhcpan_xls_with_iedb.py`**.
- **`notebooks/process_scratch.ipynb`** — TCGA matrices → filters → `data/primary_exp_*_lncRNA.csv` (edit paths inside if your raw downloads live elsewhere).

Top-level drivers (run from **repo root** so `import repo_paths` works):

| Script | What it does |
|--------|----------------|
| `generate_catalog_figures.py` | Fig 1B + 2–4A (both peptide modes when FASTA exists) |
| `generate_canonical_manuscript_figures.py` | **Single driver:** main-text Fig 1B (sklearn 2D default), OpenTSNE 1B supplement, tcga-matrix 2–3D, 4A, NetMHC 5–6 core-only; optional Fig 6 unique with ``--write-fig6-unique-supplement`` |
| `generate_netmhc_figure_bundle.py` | Canonical merged NetMHC Fig 5–6; **Fig 6 instances only** by default (add ``--also-write-unique`` for `*_unique_*` under `figures/`). Use ``--canonical-main-text-only`` for proportional-whole cohort + no random-fragment mirrors |
| `regenerate_all_figures.py` | **Full pipeline:** catalog → NetMHC bundle → Fig 5–6 supplement → canonical extras → publication export (``--skip-*``); Fig 6 unique: ``--include-fig6-unique`` |
| `rebuild_netmhc_merged_tsvs.py` | Rebuild merged TSVs from wide `*.xls` + IEDB |
| `generate_tr_lncrna_identification.py` | z-scores + R limma for Tr-lncRNA tables |
| `generate_netmhc_fig5_fig6_supplement.py` | **Fig 5–6 NetMHC supplement only** (five subfolders under `figures/supplementary/netmhc_fig5_fig6_supplement/`) |
| `generate_netmhc_supplement.py` | Optional wide 5A–5E (IC50-from-BA); unless `--skip-sensitivity`, runs **`generate_netmhc_fig5_fig6_supplement.py`** |
| `export_publication_figures.py` | Canonical 1B–6: PNG + **PDF + TIFF** under `figures/publication/`; PCA supplement after 1B unless `--skip-supplement-pca`; Fig 6 unique export only with ``--include-fig6-unique-split`` |
| `supplement/regenerate_manuscript_netmhc_figures.py` | Full clean tree under `figures/manuscript_netmhc/` |

More command-level notes: **`docs/figure_generation_overview.md`** (orchestrator map), `docs/figure_catalog.md`, `docs/netmhc_figure_commands.md`, `docs/iedb_tools_api.md`, `docs/figure6_ttn_as1_parameters.md`, `docs/CODE_REVIEW.md`, `data/netmhc/README_netmhc.md`.

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
