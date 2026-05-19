# Supplementary figure outputs

All **supplementary** PNGs and CSVs live here. **Main-text** panels are at **`figures/`** root (e.g. **`fig1_tsne_stage_lncrna_samples_dims12_panels_AB.png`** = Figure 1 combined panels 1A|1B). Plotting code is in **`manuscript/`** and **`supplement/`** at the repo root.

**Regenerate everything in this tree:** `python generate_supplementary_figures.py --strict`

See **`docs/figure_generation_overview.md`** and **`docs/figure_catalog.md`**.

---

## Layout by figure

### Figure 1 — embedding supplement

| Folder | Panels |
|--------|--------|
| **`embedding/`** | OpenTSNE supplement (`figS1_opentsne4_*_dims12_*`, `*_dims34_*`) |
| **`pca/`** | Sample PCA PC1–2 and PC3–4 (`figS_pca_*`) |

Main-text **Figure 1** (sklearn combined 1A|1B) is **`figures/fig1_tsne_stage_lncrna_samples_dims12_panels_AB.png`**.

### Figures 2 and 3 — peptide-set modes

| Folder | Panels |
|--------|--------|
| **`tcga_matrix/`** | Fig **2** peptide-fraction bars, Fig **3A** / **3B** |
| **`all_smprot_filtered/`** | Same for all-filtered SmProt; includes **Fig 3C** all-filtered heatmap when FASTA exists |

Canonical **Fig 3C** (TCGA-matrix) and **Fig 3D** are at **`figures/`** root.

### Figure 5 — NetMHC

| Folder | Panels |
|--------|--------|
| **`netmhc/coding_fragments_random_sample/`** | Random-fragment **5C** / **5D** / **5E** mirrors |
| **`netmhc_fig5_fig6_supplement/`** | Sensitivity / SB grids (five subfolders) |

Main-text Fig **5A–5E** are at **`figures/fig5_merged_*`** and **`figures/fig5de_merged_whole_*`**.

### Figure 6 — optional unique panels

| Folder | Panels |
|--------|--------|
| **`figure6_ttn_as1/`** | Fig 6 **unique** split panels (`--include-fig6-unique`) |

Main-text Fig 6 **instances** are **`figures/fig6_ttn_as1_split_instances_*.png`**.

---

## Granular commands

| Goal | Command |
|------|---------|
| Full supplement tree | `python generate_supplementary_figures.py --strict` |
| NetMHC sensitivity only | `python generate_netmhc_fig5_fig6_supplement.py --strict` |
| Main + supplement + export | `python regenerate_all_figures.py --strict` |
| Above + Fig 6 unique | `python regenerate_all_figures.py --strict --include-fig6-unique` |
