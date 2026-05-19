# lncRNA micropeptides — code (GitHub)

This repository holds the **analysis code** for the paper: manuscript figure scripts, supplement utilities, NetMHC merge helpers, and the smaller inputs that fit in git. Paths resolve through **`repo_paths.py`** from the repo root so clones are self-contained.

**Large files** are often copied from a private working directory or from the [Zenodo dataset](https://doi.org/10.5281/zenodo.20167452); a copy log may live locally at **`data/MANIFEST_COPIED.txt`**. (Wide NetMHC spreadsheets, merged `*_with_iedb.tsv`, IEDB dumps, and bulky SmProt/TCGA sources are gitignored — see **`data/README.md`**.)

## What lives where

| Folder | Role |
|--------|------|
| **`manuscript/`** | Main-text figure **scripts** (defaults write to **`figures/`**) |
| **`supplement/`** | Supplement / sensitivity **scripts** (outputs under **`figures/supplementary/`**, etc.) |
| **`figures/`** | **Main-text** figure PNGs at repo root |
| **`figures/supplementary/`** | **Supplementary** figure PNGs/CSVs — see **`figures/supplementary/README.md`** |
| **`pipeline/`** | SmProt / TCGA prep |
| **`scripts/`** | Shared NetMHC merge + SB logic (`merge_netmhcpan_xls_with_iedb.py`, `netmhc_sb_core.py`) |
| **`tr_lncrna_output/`** | Limma / DE tables |

**Why are `generate_*.py` orchestrators at the repo root?** They are entry points meant to be run as `python generate_…py` from the repository root (same cwd as `repo_paths.py`). **Plotting** code lives in **`manuscript/`** (main-text panels) and **`supplement/`** (supplement / sensitivity panels); orchestrators only call into those folders. Shared NetMHC merge/SB logic stays in **`scripts/`**; SmProt/TCGA prep in **`pipeline/`**.

## Figure regeneration

Three top-level orchestrators (run from **repo root**):

| Script | Outputs |
|--------|---------|
| **`generate_canonical_manuscript_figures.py`** | **Main text only** — `figures/` (Fig 1–4A, merged NetMHC 5–6 instances) |
| **`generate_supplementary_figures.py`** | **All supplement panels** — `figures/supplementary/` (OpenTSNE, PCA, Fig 2–3 modes, NetMHC sensitivity, etc.) |
| **`regenerate_all_figures.py`** | **Both** of the above, then optional PDF/TIFF export (`--skip-export` to omit). Add **`--include-fig6-unique`** for Fig 6 unique supplement panels + matching publication export. |

Granular drivers (NetMHC subsets, publication export, optional Fig 6 unique) are listed in **`docs/figure_generation_overview.md`**.

**Documentation map:**

| Topic | Document |
|-------|----------|
| Per-figure definitions, paths, thresholds | **`docs/figure_catalog.md`** |
| Orchestrator commands and order | **`docs/figure_generation_overview.md`** |
| Supplement folder layout | **`figures/supplementary/README.md`** |
| NetMHC merge / CLI | **`docs/netmhc_figure_commands.md`**, **`data/netmhc/README_netmhc.md`** |

NetMHC merged tables (`data/netmhc/*_with_iedb.tsv`) must exist before Fig 5–6 scripts run — see **`rebuild_netmhc_merged_tsvs.py`** and **`data/README.md`**.

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

Large tables and optional archives are on Zenodo (link above); this repo ships code, small bundled `data/`, and figure PNGs per **`.gitignore`**.

## License

Root **`LICENSE`** is MIT unless you replace it. Dataset rights are described in the Zenodo record metadata.
