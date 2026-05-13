# lncRNA micropeptides — analysis code (GitHub bundle)

Trimmed project: **canonical figure code**, **supplement / robustness code**, **data-prep pipeline**, **shared NetMHC merge + SB logic**, **docs**, and **small inputs**. Large NetMHC / IEDB tables are gitignored (see `data/README.md`).

## Layout

| Path | Contents |
|------|----------|
| **`repo_paths.py`** | `REPO_ROOT`, `DATA`, `FIGURES`, `NETMHC_DATA`, … — all scripts build paths from here (no hard-coded user directories). |
| **`manuscript/`** | Default figure scripts whose outputs go to **`figures/`** (or are mirrored there): catalog Fig 2–4A plotters, wide Fig 5A, merged Fig 5A–5E, Fig 6 TTN coverage, `manuscript_figure6_ttn_as1.py` entrypoint. |
| **`supplement/`** | Legacy wide-XLS Fig 5B–5E, cohort sensitivity, SB combination grid, Fig 6 SB sweeps, regenerate-matrix orchestrator, IEDB/NetMHC4.1 helpers, exploratory plots. |
| **`pipeline/`** | SmProt / TCGA exports, filters, NetMHC prep, summaries — not main figure panels. |
| **`scripts/`** | **`netmhc_sb_core.py`**, **`merge_netmhcpan_xls_with_iedb.py`** only (shared libraries). |
| **`generate_catalog_figures.py`** | Orchestrates **`manuscript/`** for Fig 2–4A. |
| **`generate_netmhc_figure_bundle.py`** | Orchestrates **canonical** NetMHC Fig 5–6 (merged + wide 5A + Fig 6). |
| **`generate_netmhc_supplement.py`** | Optional wide 5B–5E, sensitivity, combo grid, Fig 6 sensitivity. |

## Setup

```bash
python -m venv .venv
.venv\Scripts\activate   # Windows
pip install -r requirements.txt
```

Run commands from the **repository root** (`import repo_paths` relies on that).

## Regenerate figures

| Command | Role |
|---------|------|
| `python generate_catalog_figures.py` | Fig 2–4A → `figures/tcga_matrix/`, `figures/all_smprot_filtered/`, shared `figures/` |
| `python generate_netmhc_figure_bundle.py` | Canonical NetMHC Fig 5–6 (needs merged `*_with_iedb.tsv`; see `data/README.md`) |
| `python generate_netmhc_supplement.py` | Supplement / sensitivity / legacy wide cohort (optional flags) |
| `python supplement/regenerate_manuscript_netmhc_figures.py --help` | Full clean tree under `figures/manuscript_netmhc/` |

Docs: `docs/figure_catalog.md`, `docs/netmhc_figure_commands.md`, `data/netmhc/README_netmhc.md`.

## R

`tr_limma_de.R` lives at repo root for the limma DE step.

## License

Add a `LICENSE` when you publish the GitHub repo.
