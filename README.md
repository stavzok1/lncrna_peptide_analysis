# lncRNA micropeptides — analysis code (GitHub bundle)

Trimmed project: **canonical figure code**, **supplement / robustness code**, **data-prep pipeline**, **shared NetMHC merge + SB logic**, **docs**, **small inputs**, and **committed catalog / manuscript PNGs** under **`figures/`** (regenerate with orchestrators when numbers change). Large NetMHC / IEDB tables are gitignored (see `data/README.md`).

**Inputs / data you must provide:** see **`data/README.md`** for a checklist by orchestrator (`generate_catalog_figures.py`, `generate_netmhc_figure_bundle.py`, supplements, merge rebuild). Large files are often copied from your full `UNDEFINED` tree or Zenodo; a copy log may live at **`data/MANIFEST_COPIED.txt`**.

## Layout

| Path | Contents |
|------|----------|
| **`repo_paths.py`** | `REPO_ROOT`, `DATA`, `FIGURES`, `NETMHC_DATA`, … — all scripts build paths from here (no hard-coded user directories). |
| **`manuscript/`** | Default figure scripts whose outputs go to **`figures/`** (or are mirrored there): catalog Fig **1B**, **2–4A**, merged Fig **5–6**, `manuscript_figure6_ttn_as1.py` entrypoint. Legacy **wide-XLS Fig 5A** (IC50-from-BA) stays here for supplement use. |
| **`supplement/`** | Legacy wide-XLS Fig 5B–5E, cohort sensitivity, SB combination grid, Fig 6 SB sweeps, regenerate-matrix orchestrator, IEDB/NetMHC4.1 helpers, exploratory plots. |
| **`pipeline/`** | SmProt / TCGA exports, filters, NetMHC prep, summaries — not main figure panels. |
| **`notebooks/`** | **`process_scratch.ipynb`** — mirror of **`../process_scratch.ipynb`** (parent UNDEFINED tree): TCGA matrix → filters → stage / M_stage → `data/primary_exp_*_lncRNA.csv`. |
| **`scripts/`** | **`netmhc_sb_core.py`**, **`merge_netmhcpan_xls_with_iedb.py`** only (shared libraries). |
| **`generate_catalog_figures.py`** | Orchestrates **`manuscript/`** for Fig **1B** + **2–4A**. |
| **`generate_netmhc_figure_bundle.py`** | Orchestrates **canonical** NetMHC Fig **5–6** (merged **5A–5E** + Fig **6** only; no wide-XLS 5A). |
| **`generate_netmhc_supplement.py`** | Optional **wide 5A–5E**, sensitivity, combo grid, Fig 6 sensitivity. |
| **`figures/`** | **Shipped PNGs + CSV/txt sidecars** for Fig **1B**, **2–4A** (and `tcga_matrix/` / `all_smprot_filtered/`). **`figures/manuscript_netmhc/`** — merged NetMHC **Fig 5–6** bundle from **`supplement/regenerate_manuscript_netmhc_figures.py`** (copy or regenerate here, then commit). |

## Setup

```bash
python -m venv .venv
.venv\Scripts\activate   # Windows
pip install -r requirements.txt
Rscript install_r_dependencies.R   # R: limma + jsonlite for tr_limma_de.R (once per machine)
```

Run commands from the **repository root** (`import repo_paths` relies on that).

## Regenerate figures

| Command | Role |
|---------|------|
| `python generate_catalog_figures.py` | Fig **1B** + **2–4A** → `figures/` (1B + shared panels), `figures/tcga_matrix/`, `figures/all_smprot_filtered/` |
| `python generate_netmhc_figure_bundle.py` | Canonical NetMHC Fig **5–6** (merged cohort panels + TTN; needs merged `*_with_iedb.tsv`; see `data/README.md`) |
| `python rebuild_netmhc_merged_tsvs.py` | Rebuild merged TSVs from wide `*.xls` + IEDB (see `--help`; default compares row counts to existing files) |
| `python generate_tr_lncrna_identification.py` | Tr-lncRNA **z-scores** + **limma** from `data/primary_exp_*_lncRNA.csv` (needs **R** + **limma**); see `data/README.md` |
| `python generate_netmhc_supplement.py` | Supplement / sensitivity / legacy wide cohort (optional flags) |
| `python supplement/regenerate_manuscript_netmhc_figures.py --help` | Full clean tree under `figures/manuscript_netmhc/` |

Docs: `docs/figure_catalog.md`, `docs/netmhc_figure_commands.md`, `docs/iedb_tools_api.md`, `docs/figure6_ttn_as1_parameters.md`, `docs/ZENODO.md`, **`docs/CODE_REVIEW.md`** (bundle code audit / checklist), `data/netmhc/README_netmhc.md`.

## R

`tr_limma_de.R` lives at repo root for the limma DE step.

## License

This repository includes a root **`LICENSE`** file (**MIT** by default). Replace the copyright line with your name or institution if you prefer. If you need **non-code** terms (e.g. data-only CC BY), use a second license file or Zenodo metadata for the data deposit—see **`docs/ZENODO.md`**.

## Zenodo (DOI archive)

Short guide: **`docs/ZENODO.md`** — GitHub release integration, what to do about **gitignored** large NetMHC/IEDB files, and how to refresh figure folders from the terminal before a new upload.
