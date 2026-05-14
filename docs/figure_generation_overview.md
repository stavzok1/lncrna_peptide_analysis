# Figure generation overview

This page is the **entry map** for orchestrators: what to run, in what order, and where outputs land. For a **single command** that runs the usual full chain, see **`regenerate_all_figures.py`**. For per-figure definitions, thresholds, and file names, see **`docs/figure_catalog.md`**. NetMHC merge and WSL notes: **`data/netmhc/README_netmhc.md`**. Supplementary tree index: **`figures/supplementary/README.md`**.

---

## 1. Canonical main-text figures (single driver)

**Command:** `python generate_canonical_manuscript_figures.py`  
**Add:** `--strict` to fail on any subprocess error.

**Builds:**

- **Fig 1B** under `figures/` (default sklearn **2D t-SNE only** — two PNGs; no PC3/PC4 panels).
- **Fig 2** (tcga-matrix mode), **Fig 3A–3C** (TCGA-matrix 3C only via `--only-tcga-matrix-3c`), **Fig 3D**, **Fig 4A** — catalog scripts; canonical `fig2b` / `fig3a` / `fig3b` copies at `figures/` root; mode-specific tables and PNGs under `figures/supplementary/tcga_matrix/`.
- **Fig 5–6 NetMHC** via `generate_netmhc_figure_bundle.py --canonical-main-text-only` (merged proportional-whole cohort, Fig 6 split **instances** under `figures/`). **Fig 6** is drawn with ``plot_figure6_ttn_as1_allele_coverage.py`` whose **default gating is ``iedb_sb``** (IEDB CSV merged on ``stable_key`` + cohort-style SB bundle, including **IEDB IC50** when that column exists on the merged rows—not “local IC50 only” unless you pass ``--gating netmhc``).
- **Supplement mirrors:** OpenTSNE Fig 1B panels under `figures/supplementary/embedding/` (unless skipped); Fig 6 **unique** split panels under `figures/supplementary/figure6_ttn_as1/` **only** if you pass ``--write-fig6-unique-supplement``.

**Typical prerequisite:** expression matrix and SmProt inputs already in `data/` as described in the figure catalog.

---

## 2. Full SmProt catalog (both peptide modes)

**Command:** `python generate_catalog_figures.py`  
**Options:** `--only tcga_matrix` | `all_smprot_filtered` | `both` (default `both`); `--strict`.

**Builds:** Fig 1B, then for each mode **Fig 2** + **Fig 3A–3B** under `figures/supplementary/<mode>/`, then shared **Fig 3C–3D** + **Fig 4A** (see `manuscript/plot_figure3cd_dipeptide_log2fc_heatmaps.py` for where the all-filtered 3C panel goes).

Use this when you need the **all SmProt-filtered** universe alongside TCGA-matrix, not only the main-text subset from §1.

---

## 3. Publication exports (PDF + TIFF)

**Command:** `python export_publication_figures.py`  
**Options:** `--strict`, `--only tcga_matrix` | `all_smprot_filtered` | `both`, `--skip-netmhc`, `--skip-supplement-pca`, `--publication-tiff-kind line|color`.

Regenerates the same panels as the catalog / supplement scripts and writes **`figures/publication/...`** mirroring paths under `figures/` (including `figures/supplementary/...`). **Fig 6:** by default only **instances** under flat `figures/`; add **`--include-fig6-unique-split`** to also write **unique** under `figures/supplementary/figure6_ttn_as1/` (publication mirrors follow each PNG path).

---

## 4. General supplementary figures (non–Fig 5–6 NetMHC bundle)

Anything under **`figures/supplementary/`** that is not the dedicated NetMHC Fig 5–6 bundle (§5) is indexed in **`figures/supplementary/README.md`** — PCA, random-fragment NetMHC cohort mirrors, embedding alternates, Fig 6 unique split, z-stratum log2FC histograms, tcga_matrix / all_smprot_filtered catalog side files, etc.

Those panels are produced by the **canonical** driver (§1), **catalog** driver (§2), **`export_publication_figures.py`** (§3), or by running the listed scripts directly.

---

## 5. NetMHC — canonical merged Fig 5–6

**Command:** `python generate_netmhc_figure_bundle.py`  
**Main-text subset:** `python generate_netmhc_figure_bundle.py --canonical-main-text-only` (invoked from §1).

**Builds:** merged `*_with_iedb.tsv` workflow (unless skipped), **Fig 5A–5C** + **5D–5E** mirrors under `figures/` (and `data/netmhc/figures/` per script defaults), random-fragment cohort mirrors under `figures/supplementary/netmhc/coding_fragments_random_sample/` when the full bundle runs. **Fig 6** split panels default to **instances only** under `figures/`; pass **`--also-write-unique`** to also emit `*_unique_*` companions next to `*_instances_*`.

**Prerequisite:** NetMHC wide XLS + IEDB merge inputs in `data/netmhc/` (see **`data/netmhc/README_netmhc.md`**).

---

## 6. NetMHC — Fig 5–6 supplement (1D sensitivity, Cartesian grids, TTN sweeps)

**Dedicated folder:** **`figures/supplementary/netmhc_fig5_fig6_supplement/`**

Five **non-overlapping** roles (same default gate values as the underlying scripts; nothing is “maxed” beyond those defaults):

| Subfolder | Script | Role |
|-----------|--------|------|
| `fig5_merged_cohort_1d_sensitivity_loo/` | `supplement/netmhc_sb_sensitivity_robustness.py` | **Fig 5** merged cohort: **one dimension at a time** + leave-one-filter-out + fold-change vs baseline (**not** the full Cartesian product). |
| `fig5_merged_cohort_cartesian_sb_grid/` | `supplement/plot_fig5_netmhc_sb_combination_grid.py` | **Fig 5** merged cohort: **Cartesian product** of immuno × proc × EL × IC50 + sig/coding sharing panels. |
| `fig6_ttn_wide_netmhc_sb_sweeps/` | `supplement/plot_figure6_ttn_as1_sb_sensitivity.py` | **Fig 6** TTN: **NetMHC wide XLS only** (BA_rank / IC50-from-BA / optional EL); **no** IEDB merge. |
| `fig6_ttn_merged_iedb_1d_sensitivity_loo/` | `supplement/netmhc_ttn_merged_iedb_sb_sensitivity_robustness.py` | **Fig 6** TTN XLS + IEDB merge: **1D sweeps + LOO** (parallel in role to Fig 5’s ``netmhc_sb_sensitivity_robustness.py``). |
| `fig6_ttn_merged_iedb_cartesian_sb_grid/` | `supplement/plot_fig6_ttn_merged_iedb_sb_combination_grid.py` | **Fig 6** TTN: merged IEDB+NetMHC **Cartesian** SB grid (CSVs + **(a)** four fold heatmaps + **(b–e)** 1D sweeps at baseline + combined suite PNG). |

**Command (this folder only):**  
`python generate_netmhc_fig5_fig6_supplement.py`  
Add `--strict` if any step should abort the run.

**Command (wide legacy 5A–5E + this bundle):**  
`python generate_netmhc_supplement.py`  
Optional: `--include-wide-xls-fig5` (IC50-from-BA wide XLS cohort). Use `--skip-sensitivity` to skip §6 entirely.

**Recommendation:** run **`generate_netmhc_figure_bundle.py`** (§5) first so merged tables exist for the Fig 5 supplement scripts.

---

## 7. Optional clean NetMHC “matrix” under `figures/manuscript_netmhc/`

**Command:** `python supplement/regenerate_manuscript_netmhc_figures.py --help`  
Rebuilds a structured tree of SB modes, instances/unique, sensitivity copies, etc. This is for **archival / comparison**, not the same as the flat `figures/` canonical mirrors. See **`docs/figure_catalog.md`**.

---

## Full regeneration (“everything”, from repo root)

**Prerequisites:** inputs under `data/` as described in **`data/README.md`**, **`data/netmhc/README_netmhc.md`**, and **`docs/figure_catalog.md`** (expression matrices, SmProt FASTAs, merged `*_with_iedb.tsv`, NetMHC XLS, TTN IEDB companion CSV, etc.). Large files may be gitignored—copy from Zenodo or your working tree first.

**One command (same order as below):**  
`python regenerate_all_figures.py --strict`  
Optional skips: `--skip-catalog`, `--skip-netmhc-bundle`, `--skip-netmhc-fig5-fig6-supplement`, `--skip-canonical`, `--skip-export`.  
Optional Fig 6 **unique** (canonical + publication): ``--include-fig6-unique``.

**Core pipeline (recommended order):**

1. **SmProt catalog — both peptide modes + shared Fig 3C–3D + 4A:**  
   `python generate_catalog_figures.py --strict`
2. **NetMHC canonical Fig 5–6 + random-fragment mirrors:**  
   `python generate_netmhc_figure_bundle.py --strict`
3. **NetMHC Fig 5–6 supplement tree** (1D + LOO cohort, Cartesian Fig 5, TTN NetMHC sweeps, TTN merged IEDB 1D+LOO, TTN merged IEDB Cartesian):  
   `python generate_netmhc_fig5_fig6_supplement.py --strict`
4. **Main-text extras** (OpenTSNE Fig 1B supplement, optional Fig 6 **unique** via ``--write-fig6-unique-supplement``, and a **proportional-whole-only** NetMHC pass via `--canonical-main-text-only`):  
   `python generate_canonical_manuscript_figures.py --strict`  
   *Safe after steps 1–3; overlaps parts of 1–2 but refreshes embedding / optional unique Fig 6.*
5. **Publication PDF + TIFF:**  
   `python export_publication_figures.py --strict`  
   *(Add ``--include-fig6-unique-split`` if you also want Fig 6 unique publication exports.)*

**Optional add-ons:**

- **Legacy wide-XLS Fig 5A–5E (IC50-from-BA):**  
  `python generate_netmhc_supplement.py --strict --include-wide-xls-fig5`  
  (still runs step 3’s bundle unless you pass `--skip-sensitivity`.)
- **PCA supplement** (also run by `export_publication_figures.py` unless `--skip-supplement-pca`):  
  `python supplement/plot_supplement_pca_stage_samples.py --out-dir figures/supplementary/pca`
- **Z-stratum log2FC exploratory histograms:**  
  `python supplement/plot_z_stratum_logfc_histograms.py`
- **Clean archival tree under `figures/manuscript_netmhc/`:**  
  `python supplement/regenerate_manuscript_netmhc_figures.py --help` (see flags such as `--clean`).

---

## Suggested minimal refresh (main text only)

1. `python generate_canonical_manuscript_figures.py --strict`  
2. `python export_publication_figures.py --strict --skip-netmhc` *or* full export if NetMHC panels unchanged.
