# lncRNA micropeptides — code (GitHub)

This repository holds the **analysis code** for the paper: manuscript figure scripts, supplement utilities, NetMHC merge helpers, and the smaller inputs that fit in git. Paths resolve through **`repo_paths.py`** from the repo root so clones are self-contained.

**Large files** are often copied from a private working directory or from Zenodo; a copy log may live at **`data/MANIFEST_COPIED.txt`**. (Wide NetMHC spreadsheets, merged `*_with_iedb.tsv`, IEDB dumps, and bulky SmProt/TCGA sources are usually gitignored.) File-by-file expectations: **`data/README.md`** and **`docs/ZENODO.md`**.

## What lives where

| Folder | What it is | Figure files? |
|--------|------------|---------------|
| **`manuscript/`** | Main-text figure **scripts** (Fig 1B–4A, merged NetMHC 5–6). Defaults write canonical PNGs to **`figures/`**. | No — code only |
| **`supplement/`** | Supplement / sensitivity / legacy wide-XLS **scripts** and small helpers (IEDB fetch, tables). | No — code only; outputs go under **`figures/supplementary/`**, **`figures/`**, or **`data/netmhc/figures/`** |
| **`figures/supplementary/`** | **Supplementary figure outputs** (PNG/CSV), grouped by topic. Index: **`figures/supplementary/README.md`**. | Yes — this is where supplement *figures* live |
| **`pipeline/`** | SmProt / TCGA prep (filters, significant-lnc lists). | No |
| **`scripts/`** | Shared NetMHC **library** code (see below). | No |
| **`tr_lncrna_output/`** | Limma / DE tables and some legacy plots from the R/Python DE step. | Some PNGs (DE diagnostics) |

**`scripts/` (shared NetMHC logic, not figure scripts):**

- **`merge_netmhcpan_xls_with_iedb.py`** — builds merged long tables `data/netmhc/*_with_iedb.tsv` from wide NetMHCpan `*.xls` + IEDB peptide_table CSVs. Used by **`rebuild_netmhc_merged_tsvs.py`** and indirectly before any **merged** Fig 5–6 run. Manuscript/supplement plot scripts **read** the TSVs; they do not call the merger unless you run it yourself.
- **`netmhc_sb_core.py`** — one definition of SB gates (immuno, processing, EL, IC50 defaults). Imported by **`manuscript/plot_fig5*.py`**, **`manuscript/plot_figure6_*.py`**, and **`supplement/*sensitivity*`**, **`*combination_grid*`**, etc., so main text and supplement use the same thresholds.

**`notebooks/process_scratch.ipynb`** — TCGA matrices → `data/primary_exp_*_lncRNA.csv` (paths are local to your machine).

---

## Which driver script to run

Think in three layers: **(A) main text**, **(B) supplement figures**, **(C) optional extras**. You rarely need every top-level script.

### A — Main-text figures (usually enough for the paper)

| Script | Role | Needed? |
|--------|------|---------|
| **`generate_canonical_manuscript_figures.py`** | One-shot **main text**: Fig 1B (sklearn), tcga-matrix 2–3D, 4A, merged NetMHC 5–6 (instances), OpenTSNE 1B under `figures/supplementary/embedding/`. | **Yes** — primary entry for reviewers cloning the repo |
| **`generate_netmhc_figure_bundle.py`** | Merged NetMHC **5A–5E + Fig 6** only (also called from the canonical driver). Use alone if you only refresh NetMHC. | **Yes** (directly or via canonical driver) |
| **`generate_catalog_figures.py`** | Broader catalog: **both** peptide modes (`tcga_matrix` + `all_smprot_filtered`) for Fig 2–3. Overlaps canonical driver for tcga-matrix. | Only if you need the all-filtered universe on GitHub |

**Prerequisite for NetMHC:** merged `*_with_iedb.tsv` in `data/netmhc/` → `python rebuild_netmhc_merged_tsvs.py` (or copy from Zenodo).

### B — Supplement figures (sensitivity / grids)

| Script | Role | Needed? |
|--------|------|---------|
| **`generate_netmhc_fig5_fig6_supplement.py`** | Only the five folders under **`figures/supplementary/netmhc_fig5_fig6_supplement/`** (1D+LOO, Cartesian grids, TTN sweeps). | **Yes** if those supplement panels are in the paper |
| **`generate_netmhc_supplement.py`** | **Deprecated wrapper.** Legacy wide-XLS Fig 5A–5E only with `--include-wide-xls-fig5`; otherwise just calls `generate_netmhc_fig5_fig6_supplement.py`. Prefer that script or `regenerate_all_figures.py` instead. | **Legacy only** — emits `DeprecationWarning` on run |

### C — Optional / local-only

| Script | Role | Needed? |
|--------|------|---------|
| **`regenerate_all_figures.py`** | Runs A + B + OpenTSNE/canonical pass + **`export_publication_figures.py`** in order. Convenience for a full refresh. | Optional convenience |
| **`export_publication_figures.py`** | Regenerates panels and writes **PDF + TIFF** under `figures/publication/` (gitignored; for journal submission locally). | Optional — not on GitHub |
| **`supplement/regenerate_manuscript_netmhc_figures.py`** | Rebuilds an archival matrix under **`figures/manuscript_netmhc/`** (all SB modes × instances/unique). For method comparison, not the flat `figures/` layout. | Optional — developer/archive tool |
| **`rebuild_netmhc_merged_tsvs.py`** | Data prep only (merge XLS + IEDB). | When merged TSVs are missing |
| **`generate_tr_lncrna_identification.py`** | z-scores + R limma → `tr_lncrna_output/`, canonical gene list. | When (re)building DE inputs |

**Typical minimal refresh:** `generate_canonical_manuscript_figures.py --strict` then, if needed, `generate_netmhc_fig5_fig6_supplement.py --strict`.

**Typical full refresh:** `regenerate_all_figures.py --strict` (add `--skip-export` if you do not need PDF/TIFF).

More detail: **`docs/figure_generation_overview.md`**, **`docs/figure_catalog.md`**, **`docs/netmhc_figure_commands.md`**, **`data/netmhc/README_netmhc.md`**.

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
