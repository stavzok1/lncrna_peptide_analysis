# Figure generation overview

Entry map for orchestrators. Per-figure paths and thresholds: **`docs/figure_catalog.md`**. Supplement folder index: **`figures/supplementary/README.md`**. NetMHC data: **`data/netmhc/README_netmhc.md`**.

---

## Three main orchestrators

| Command | Scope |
|---------|--------|
| `python generate_canonical_manuscript_figures.py` | **Main text** → `figures/` (Fig 1–4A, merged NetMHC 5–6 instances). |
| `python generate_supplementary_figures.py` | **All supplement** → `figures/supplementary/` (OpenTSNE Fig 1, PCA, Fig 2–3 modes, NetMHC sensitivity, etc.). |
| `python regenerate_all_figures.py` | **Both** above + `export_publication_figures.py` (PDF/TIFF under `figures/publication/`, gitignored). Use `--skip-export` to omit step 3; **`--include-fig6-unique`** for Fig 6 unique under `figures/supplementary/figure6_ttn_as1/`. |

Add `--strict` to exit non-zero on the first failure.

---

## What `generate_supplementary_figures.py` runs

| Step | Output | Script(s) |
|------|--------|-----------|
| OpenTSNE Fig 1 supplement | `figures/supplementary/embedding/` | `manuscript/plot_figure1_tsne_stage_lncrna.py` (`opentsne4`, prefix `figS1_opentsne4_…`) |
| PCA | `figures/supplementary/pca/` | `supplement/plot_supplement_pca_stage_samples.py` |
| Fig 2–3 (both peptide modes) | `tcga_matrix/`, `all_smprot_filtered/` | Inlined: `plot_tr_de_peptide_fractions_by_transition.py`, `plot_aa_frequency_tcga_vs_proteome.py`, `plot_dipeptide_volcano_lnc_vs_proteome.py`, `plot_figure3cd` (`--only-all-smprot-filtered-3c`) |
| NetMHC random-fragment mirrors | `netmhc/coding_fragments_random_sample/` | `generate_netmhc_figure_bundle.py --supplement-mirrors-only` |
| NetMHC sensitivity / grids | `netmhc_fig5_fig6_supplement/` | `generate_netmhc_fig5_fig6_supplement.py` |
| Fig 6 unique (optional) | `figure6_ttn_as1/` | `plot_figure6_ttn_as1_allele_coverage.py` (`--include-fig6-unique`) |

---

## Granular orchestrators

| Command | Role |
|---------|------|
| `python generate_netmhc_figure_bundle.py` | Main-text NetMHC 5–6 + random-fragment mirrors (full bundle). `--canonical-main-text-only` or `--supplement-mirrors-only` for subsets. |
| `python generate_netmhc_fig5_fig6_supplement.py` | NetMHC sensitivity tree only. |
| `python generate_netmhc_supplement.py` | **Deprecated** — legacy wide-XLS Fig 5A–5E with `--include-wide-xls-fig5`. |
| `python export_publication_figures.py` | PDF + TIFF export (local; gitignored). |
| `python rebuild_netmhc_merged_tsvs.py` | Merge wide `*.xls` + IEDB → `*_with_iedb.tsv`. |
| `python supplement/regenerate_manuscript_netmhc_figures.py` | Archival tree under `figures/manuscript_netmhc/`. |

---

## Suggested workflows

**Main text only:** `python generate_canonical_manuscript_figures.py --strict`

**Supplement only:** `python generate_supplementary_figures.py --strict`

**Everything:** `python regenerate_all_figures.py --strict`

**Everything + Fig 6 unique supplement:** `python regenerate_all_figures.py --strict --include-fig6-unique`

**NetMHC sensitivity only:** `python generate_netmhc_fig5_fig6_supplement.py --strict`
