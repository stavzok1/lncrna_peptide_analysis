# Supplementary figure outputs

All **supplementary** PNGs and CSVs live here. **Main-text** panels are at **`figures/`** root. Plotting code is in **`manuscript/`** and **`supplement/`** at the repo root.

**Regenerate everything in this tree:**  
`python generate_supplementary_figures.py --strict`

See **`docs/figure_generation_overview.md`** for step-by-step script mapping and **`docs/figure_catalog.md`** for panel definitions.

---

## Layout by figure

### Figure 1 — embedding (TCGA stage lncRNA matrix)

| Folder | Panels |
|--------|--------|
| **`embedding/`** | OpenTSNE supplement: dims **1–2** and **3–4** (`figS1b_opentsne4_*_dims12_*`, `*_dims34_*`) |
| **`pca/`** | Sample PCA PC1–2 and PC3–4 (`figS_pca_stage_lncrna_samples_pc*_*.png`) |

Main-text Fig 1B (sklearn combined A|B) is **`figures/fig1b_*_panels_AB.png`**, not in this tree.

### Figures 2 and 3 — peptide-set modes

| Folder | Panels |
|--------|--------|
| **`tcga_matrix/`** | Fig **2** peptide-fraction bars (`peptide_fraction/`), Fig **3A** 1-mer frequency, Fig **3B** volcano (+ CSV sidecars) |
| **`all_smprot_filtered/`** | Same for the all-filtered SmProt universe; includes **Fig 3C** all-filtered heatmap when the FASTA exists |

Canonical **Fig 3C** (TCGA-matrix) and **Fig 3D** are at **`figures/`** root.

### Figure 5 — NetMHC

| Folder | Panels |
|--------|--------|
| **`netmhc/coding_fragments_random_sample/`** | Random coding-fragment cohort **5C** / **5D** / **5E** mirrors |
| **`netmhc_fig5_fig6_supplement/`** | Sensitivity and SB threshold grids (five subfolders — see below) |

Main-text Fig **5A–5E** are at **`figures/fig5_merged_*`** and **`figures/fig5de_merged_whole_*`**.

**`netmhc_fig5_fig6_supplement/`** subfolders:

| Subfolder | Content |
|-----------|---------|
| `fig5_merged_cohort_1d_sensitivity_loo/` | Fig 5 merged cohort: 1D SB sweeps + leave-one-out |
| `fig5_merged_cohort_cartesian_sb_grid/` | Fig 5 merged cohort: Cartesian SB grid |
| `fig6_ttn_wide_netmhc_sb_sweeps/` | Fig 6 TTN: NetMHC wide XLS sweeps only |
| `fig6_ttn_merged_iedb_1d_sensitivity_loo/` | Fig 6 TTN: merged IEDB+NetMHC 1D sweeps + LOO |
| `fig6_ttn_merged_iedb_cartesian_sb_grid/` | Fig 6 TTN: merged IEDB+NetMHC Cartesian grid |

Regenerate this subtree only: `python generate_netmhc_fig5_fig6_supplement.py --strict`

### Figure 6 — optional unique panels

| Folder | Panels |
|--------|--------|
| **`figure6_ttn_as1/`** | Fig 6 **unique** split panels (optional; `--include-fig6-unique` on `generate_supplementary_figures.py`) |

Main-text Fig 6 **instances** are **`figures/fig6_ttn_as1_split_instances_*.png`**.

---

## Granular commands

| Goal | Command |
|------|---------|
| Full supplement tree | `python generate_supplementary_figures.py --strict` |
| NetMHC sensitivity only | `python generate_netmhc_fig5_fig6_supplement.py --strict` |
| PCA only | `python supplement/plot_supplement_pca_stage_samples.py` |
| Fig 2–3 both modes only | `python generate_catalog_figures.py --supplement-only --strict` |
| Main + supplement + export | `python regenerate_all_figures.py --strict` |
