# Figure generation overview

Entry map for orchestrators. Per-figure paths and thresholds: **`docs/figure_catalog.md`**. Supplement folder index: **`figures/supplementary/README.md`**. NetMHC data: **`data/netmhc/README_netmhc.md`**.

---

## Three main orchestrators

| Command | Scope |
|---------|--------|
| `python generate_canonical_manuscript_figures.py` | **Main text** → `figures/` (Fig 1B–4A, merged NetMHC 5–6 instances). **No** OpenTSNE / PCA / supplement trees. |
| `python generate_supplementary_figures.py` | **All supplement** → `figures/supplementary/` (OpenTSNE, PCA, Fig 2–3 both modes, NetMHC random-fragment mirrors, NetMHC sensitivity tree; optional Fig 6 unique with `--include-fig6-unique`). |
| `python regenerate_all_figures.py` | **Both** above + `export_publication_figures.py` (PDF/TIFF under `figures/publication/`, gitignored). Use `--skip-export` to omit step 3. |

Add `--strict` to any of the above to exit non-zero on the first failure.

**Prerequisites:** inputs under `data/` as in **`data/README.md`** and **`docs/figure_catalog.md`** (expression matrices, SmProt FASTAs, merged `*_with_iedb.tsv`, etc.). Large files may be gitignored — copy from Zenodo first.

---

## What `generate_supplementary_figures.py` runs

| Step | Output | Underlying script(s) |
|------|--------|----------------------|
| OpenTSNE Fig 1B | `figures/supplementary/embedding/` | `manuscript/plot_figure1b_tsne_stage_lncrna.py` (`opentsne4`) |
| PCA | `figures/supplementary/pca/` | `supplement/plot_supplement_pca_stage_samples.py` |
| Fig 2–3 (both peptide modes) | `tcga_matrix/`, `all_smprot_filtered/` | `generate_catalog_figures.py --supplement-only` |
| NetMHC random-fragment mirrors | `netmhc/coding_fragments_random_sample/` | `generate_netmhc_figure_bundle.py --supplement-mirrors-only` |
| NetMHC sensitivity / grids | `netmhc_fig5_fig6_supplement/` | `generate_netmhc_fig5_fig6_supplement.py` |
| Fig 6 unique (optional) | `figure6_ttn_as1/` | `manuscript/plot_figure6_ttn_as1_allele_coverage.py` |

Skip flags: `--skip-netmhc-fig5-fig6-supplement`, `--skip-netmhc-random-fragments`.

---

## Granular orchestrators

| Command | Role |
|---------|------|
| `python generate_netmhc_figure_bundle.py` | Main-text NetMHC 5–6 + random-fragment mirrors (full bundle). `--canonical-main-text-only` = main text only; `--supplement-mirrors-only` = random-fragment mirrors only. |
| `python generate_netmhc_fig5_fig6_supplement.py` | NetMHC sensitivity tree only (`netmhc_fig5_fig6_supplement/`). |
| `python generate_catalog_figures.py` | Legacy **combined** catalog (main + supplement overlap). Prefer canonical + supplementary drivers. `--supplement-only` = supplement Fig 2–3 trees only. |
| `python generate_netmhc_supplement.py` | **Deprecated** — legacy wide-XLS Fig 5A–5E with `--include-wide-xls-fig5`; otherwise forwards to `generate_netmhc_fig5_fig6_supplement.py`. |
| `python export_publication_figures.py` | Re-runs plotting + PDF/TIFF export (local; gitignored). |
| `python rebuild_netmhc_merged_tsvs.py` | Merge wide `*.xls` + IEDB → `*_with_iedb.tsv`. |
| `python supplement/regenerate_manuscript_netmhc_figures.py` | Archival tree under `figures/manuscript_netmhc/`. |

NetMHC Fig 5–6 supplement subfolders (1D+LOO, Cartesian grids, TTN sweeps) are documented in **`docs/figure_catalog.md`** § Figure 5–6 supplement.

---

## Suggested workflows

**Main text only:**  
`python generate_canonical_manuscript_figures.py --strict`

**Supplement only:**  
`python generate_supplementary_figures.py --strict`

**Everything (PNG + optional publication export):**  
`python regenerate_all_figures.py --strict`  
Add `--include-fig6-unique` for Fig 6 unique supplement panels and publication mirrors.

**NetMHC sensitivity tree only:**  
`python generate_netmhc_fig5_fig6_supplement.py --strict`  
(requires merged TSVs; run after main-text NetMHC or `rebuild_netmhc_merged_tsvs.py`)
