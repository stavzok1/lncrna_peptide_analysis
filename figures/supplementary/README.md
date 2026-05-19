# Supplementary figure outputs

**All supplementary PNGs/CSVs live here** (or under `data/netmhc/figures/` for NetMHC sidecars). The Python code that *creates* them is in **`supplement/`** and **`manuscript/`** at the repo root — not in this folder.

## Layout

| Path | Contents |
|------|----------|
| `tcga_matrix/` | Fig 2 peptide-fraction bars, Fig 3A–3B, tables (TCGA-matrix peptide set) |
| `all_smprot_filtered/` | Same panels for the all-filtered SmProt universe (when FASTA exists) |
| `embedding/` | OpenTSNE Fig 1B alternate (`figS1b_opentsne4_*`) |
| `pca/` | PCA supplement on the same matrix as Fig 1B (`figS_pca_*`) |
| `netmhc/coding_fragments_random_sample/` | Random-fragment coding cohort mirrors for Fig 5C / 5E |
| `netmhc_fig5_fig6_supplement/` | NetMHC **sensitivity** bundle (1D+LOO, Cartesian grids, TTN sweeps) — five subfolders |
| `figure6_ttn_as1/` | Optional Fig 6 **unique** split panels (not main-text instances) |

Main-text canonical PNGs stay at **`figures/`** root (e.g. `fig5_merged_5a_*.png`, `fig6_ttn_as1_split_instances_*.png`).

## How to regenerate

See **`docs/figure_generation_overview.md`**. Short version:

- Catalog side files: `python generate_catalog_figures.py`
- NetMHC sensitivity tree only: `python generate_netmhc_fig5_fig6_supplement.py` (after merged `*_with_iedb.tsv` exist; do not use deprecated `generate_netmhc_supplement.py` for this alone)
- Full chain: `python regenerate_all_figures.py`
