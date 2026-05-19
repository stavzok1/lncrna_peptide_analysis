# Figure catalog

Structured notes for selected manuscript-style figures: what they show, how they are generated, inputs and outputs, inferential procedures (if any), and limitations. **Which script to run first:** see **`docs/figure_generation_overview.md`**.

**Parameter index:** for a table of *where* each figure’s thresholds and CLI are documented (including **`analysis_params.md`** deep notes for Fig 2 bar logic), see **`docs/analysis_params.md`** → *Where each figure’s parameters are documented*.

**Script layout in this repository:** canonical figure code is under **`manuscript/`** (outputs default to **`figures/`** at repo root; catalog **mode-specific** tables and PNGs default to **`figures/supplementary/<mode>/`**). Supplement / sensitivity / legacy wide cohort code is under **`supplement/`**. Shared NetMHC merge + SB logic is under **`scripts/`** (`merge_netmhcpan_xls_with_iedb.py`, `netmhc_sb_core.py`). SmProt / TCGA prep utilities are under **`pipeline/`**. Orchestrators at repo root call into `manuscript/`. Table entries below list **basename** for readability; run as e.g. ``python manuscript/plot_tr_de_peptide_fractions_by_transition.py``.

**Output root:** catalog scripts for **Fig 2** and **Fig 3A–3B** write under **`figures/supplementary/<peptide_mode>/`**, where **`peptide_mode`** is **`tcga_matrix`** or **`all_smprot_filtered`**, except **Figure 1B** (sample embedding) which writes directly under **`figures/`**. Canonical copies **`fig2b_stage_E_L_combined.png`**, **`fig3a.png`**, and **`fig3b.png`** (TCGA-matrix mode) are also written at **`figures/`** root.

**Figure 1B** — sample embedding of TCGA primary samples on the lncRNA stage matrix (`plot_figure1b_tsne_stage_lncrna.py`; default **2D sklearn t-SNE only**, two PNGs under **`figures/`**). **Figure 2** — peptide-fraction bars (TCGA limma–z). **Figure 3** — composition vs `known_proteins.fasta` (**3A** 1-mer, **3B** volcano, **3C** log2FC dipeptide heatmaps: **TCGA-matrix** file at repo **`figures/`** root; **all-filtered** 3C under **`figures/supplementary/all_smprot_filtered/`** when that FASTA exists, **3D** Tr log2FC heatmap at repo root). **Figure 4A** — TIS vs Ribo-seq p-values for **501** analyzed MPs (`data/significant_lnc_peptides.tsv`; `plot_figure4a_tis_vs_ribo_tr_mps.py`, no flags). For **TCGA-matrix-only** 3C (omit all-filtered 3C), use ``plot_figure3cd_dipeptide_log2fc_heatmaps.py --only-tcga-matrix-3c`` (as in ``generate_canonical_manuscript_figures.py``).

---

## Figure 1B — t-SNE panels (TCGA primary, lncRNA stage matrix)

**Paths:** `figures/fig1b_tsne_stage_lncrna_samples_dims12_*.png` (**two** panels: cancer type, AJCC stage) with default `--embedding sklearn2_pca34`. With `--embedding opentsne4`, also `*_dims34_*.png` (four panels total). Alternate **OpenTSNE** supplement defaults to `figures/supplementary/embedding/figS1b_opentsne4_tsne_stage_lncrna_samples_*.png` when using `generate_canonical_manuscript_figures.py` (override basename with `plot_figure1b_tsne_stage_lncrna.py --filename-prefix`).

**Script:** `plot_figure1b_tsne_stage_lncrna.py` — input **`data/primary_exp_stage_lncRNA.csv`**. **Default** (`--embedding sklearn2_pca34`): **sklearn** `TSNE` **n_components=2** (Barnes–Hut) on `X_in` (after optional `--n-pca` gene-level truncation). **Four t-SNE panels:** `--embedding opentsne4` (requires **`pip install opentsne`**). Same matrix columns as the limma / peptide-fraction pipeline. Use **`--filename-prefix`** and **`--out-dir`** to write non-default embeddings without overwriting canonical basenames.

---

## Supplement — PCA PC1–PC2 and PC3–PC4 (same matrix as Fig 1B)

**Paths:** `figures/supplementary/pca/figS_pca_stage_lncrna_samples_pc*_*.png` (four panels; default output of `supplement/plot_supplement_pca_stage_samples.py`).

**Script:** `supplement/plot_supplement_pca_stage_samples.py` — **sample PCA** on standardized expression (optional `--n-pca` gene truncation before PCA, default off). Independent of Fig 1B t-SNE; use for a straightforward linear low-dimensional view.

---

## Figure 2 — Peptide-producing fraction among DE ∩ high-|z| genes, by cancer

**Paths:** `figures/supplementary/<peptide_gene_set>/peptide_fraction/*.png` and `figures/supplementary/<peptide_gene_set>/tr_de_peptide_fraction_by_cancer.csv`

| `peptide_gene_set` | Peptide TSV |
|--------------------|-------------|
| `tcga_matrix` (default) | `data/smprot_filtered_tcga_expr_genes.tsv` |
| `all_smprot_filtered` | `data/smprot_filtered.tsv` |

**Script:** `plot_tr_de_peptide_fractions_by_transition.py` — optional **`--figures-dir`** to override the `figures/supplementary/` parent.

---

## Figure 3 — lncRNA peptide composition vs **`data/known_proteins.fasta`**

### Figure 3A — Pooled 1-mer frequency

**Paths:** `figures/supplementary/<peptide_set>/aa_frequency_*_vs_known_proteins.{png,csv,txt}`

**Script:** `plot_aa_frequency_tcga_vs_proteome.py` — default **`tcga_matrix`**; **`--peptide-set all_smprot_filtered`** needs **`data/smprot_all_filtered_peptides.faa`**.

### Figure 3B — Dipeptide volcano

**Paths:** `figures/supplementary/<peptide_set>/dipeptide_volcano_lnc_vs_proteome*.{png,csv,txt}`

**Script:** `plot_dipeptide_volcano_lnc_vs_proteome.py`

### Figure 3C — Dipeptide log2FC heatmaps (TCGA-matrix and all-filtered, **separate files**)

**Paths:**

- `figures/fig3c_dipeptide_log2fc_tcga_matrix_vs_proteome.png` — always (TCGA-matrix FASTA).
- `figures/supplementary/all_smprot_filtered/fig3c_dipeptide_log2fc_all_smprot_filtered_vs_proteome.png` — when
  `data/smprot_all_filtered_peptides.faa` exists (unless a custom `--out-dir` is set, in which case both 3C panels use that directory).

Each file has its own symmetric color scale (99th percentile of |log2 ratio|, floor 0.5).
**Script:** `plot_figure3cd_dipeptide_log2fc_heatmaps.py` — same +0.5 / 400-cell smoothing and
`log2(peptide/proteome)` as `plot_dipeptide_mp_figure2.py`.

### Figure 3D — Dipeptide log2FC: Tr-lncRNA MPs (TCGA-matrix) vs proteome

**Path:** `figures/fig3d_dipeptide_log2fc_tr_lncrna_tcga_vs_proteome.png`

**Script:** `plot_figure3cd_dipeptide_log2fc_heatmaps.py` (second output). Restricts
`data/smprot_tcga_filtered_peptides.faa` to **canonical** Tr gene symbols
(`canonical_significant_lncRNAs.txt` or `limma_z_intersection_genes.txt`).

**Legacy multi-panel figure:** `plot_dipeptide_mp_figure2.py` → `tr_lncrna_output/figures/fig2_dipeptide_mp_composition.png`

### Figure 4A — TIS vs Ribo-seq p-values (canonical Tr-lncRNA MPs)

**Path:** `figures/fig4a_tr_lncrna_mp_tis_vs_riboseq_pvalues.png`

**Script:** `plot_figure4a_tis_vs_ribo_tr_mps.py` — default: all rows in **`data/significant_lnc_peptides.tsv`**
(~501 exportable MPs; NetMHC cohort), both **TISPvalue** and **RiboPvalue** ≤ 0.05. Optional `--cohort tr_lncrna`
filters by **GeneSymbol or GeneID** (Ensembl) vs **`significant_lncs.csv`** / canonical Tr genes (`smprot_gene_match.py`).
Log-scaled
axes **10⁻¹²–10⁻¹** on both dimensions (ticks do not extend below 10⁻¹²). Shading: **green**
(TIS ≤ 10⁻⁴, Ribo ≥ 10⁻⁴), **blue** (Ribo ≤ 10⁻⁴, TIS ≥ 10⁻⁴), **violet** (both &lt; 10⁻⁴);
polygons and p = 10⁻⁴ reference lines extend to **10⁻¹**. Orange points; framed labels for
**PTPRG-AS1**, **LINC00326**, and **LINC00958** only inside those shaded regions; most
TIS-extreme **LINC00958** MP as red dot with arrow callout; extra plain-text labels for the
top **N** combined-significance MPs (`--top-extreme-labels`, default 28).

---

## Index

| ID | Short name | Script |
|----|--------------|--------|
| 1B | t-SNE panels (stage lncRNA matrix) | `plot_figure1b_tsne_stage_lncrna.py` |
| S | PCA PC1–2 & PC3–4 (supplement) | `supplement/plot_supplement_pca_stage_samples.py` |
| 2 | Peptide fraction by cancer × transition | `plot_tr_de_peptide_fractions_by_transition.py` |
| 3A | 1-mer AA vs proteome | `plot_aa_frequency_tcga_vs_proteome.py` |
| 3B | Dipeptide volcano vs proteome | `plot_dipeptide_volcano_lnc_vs_proteome.py` |
| 3C–3D | Dipeptide log2FC heatmaps (3C split files / 3D Tr) | `plot_figure3cd_dipeptide_log2fc_heatmaps.py` |
| 4A | TIS vs Ribo-seq p scatter (canonical Tr MPs) | `plot_figure4a_tis_vs_ribo_tr_mps.py` |
| — | Fig 1B + both modes + 3C–3D + 4A (orchestrator) | `generate_catalog_figures.py` |
| — | Main-text 1B (sklearn) + tcga-matrix 2–4A + NetMHC core + supplements (OpenTSNE 1B; optional Fig 6 unique) | `generate_canonical_manuscript_figures.py` |
| — | Fig 5–6 NetMHC supplement tree (5 subfolders under `figures/supplementary/netmhc_fig5_fig6_supplement/`) | `generate_netmhc_fig5_fig6_supplement.py` |
| — | Fig 6 TTN merged IEDB 1D + LOO | `supplement/netmhc_ttn_merged_iedb_sb_sensitivity_robustness.py` |
| — | Fig 6 TTN merged IEDB Cartesian SB grid | `supplement/plot_fig6_ttn_merged_iedb_sb_combination_grid.py` |

---

## Figure 5 — NetMHCpan cohort (manuscript: merged IEDB+NetMHC, **SB default**)

**Manuscript stance:** Figure 5 in this catalog refers to the **merged** ``*_with_iedb.tsv`` pipeline
with **``sb_mode=full``** (default) and thresholds shared with Figure 6 IEDB gating:
**``FIG5_IEDB_*``** in ``scripts/netmhc_sb_core.py`` (immunogenicity **> 0.1**, processing **> 1.5**,
EL **< 1%** by default, IC50 **< 150 nM** using the IEDB BA-IC50 column when present, else local
``BA_score``). Merge inputs with ``scripts/merge_netmhcpan_xls_with_iedb.py``.

**Sig lnc cohort (all panels):** ``data/netmhc/netmhcpan_sig_lnc_with_iedb.tsv``.

**Coding cohorts (merged TSVs; two **5C** outputs in `data/netmhc/figures/`, distinct **repo mirrors**):**

| Cohort | Merged TSV | Default **output stem** (`--out-dir` default `data/netmhc/figures`) | Repo mirror for PNG copies |
|--------|------------|----------------------------------------------------------------------|----------------------------|
| **Proportional whole** (length-matched parents) | ``data/netmhc/netmhcpan_coding_proportional_whole_with_iedb.tsv`` | ``fig5_merged`` (same run as canonical **5A–5B**) | ``figures/`` (with **5A–5C**) |
| **Random coding-fragment control** | ``data/netmhc/netmhcpan_coding_control_with_iedb.tsv`` | ``figS5c_random_fragments`` (``plot_fig5abc_netmhc_sb_triple.py --panels c``) | ``figures/supplementary/netmhc/coding_fragments_random_sample/`` (**5C** only) |

**5B** is **significant lncRNA only** (one canonical stem per run); the second coding cohort does not re-emit 5B.

**Scripts:** ``manuscript/plot_fig5abc_netmhc_sb_triple.py`` (5A–5C) and ``manuscript/plot_fig5de_merged_iedb_sb_per_allele.py`` (5D–5E; run twice with ``--coding-tsv`` + ``--output-stem`` + optional ``--repo-mirror-dir`` for proportional-whole vs random-fragment merged coding). Default plot metrics use **SB row instances** (``--fig5a-y-metric`` / ``--sharing-y-metric`` / ``--count-metric`` can switch to **unique**).

**Repo mirrors:** PNGs are written under ``data/netmhc/figures/`` and copied into **``figures/``** or **``figures/supplementary/...``** unless ``--no-repo-mirror`` is passed. Use ``--repo-mirror-dir DIR`` to override the default mirror root (canonical merged triple defaults to ``figures/``).

| Panel | Uses coding cohort? | Mirrored PNG location (merged **full** SB; default / supplement stems) |
|-------|---------------------|------------------------------------------------------------------------|
| **5A** | No (sig lnc SB only) | ``figures/fig5_merged_5a_epitopes_vs_allele_frequency.png`` |
| **5B** | No | ``figures/fig5_merged_5b_epitope_sharing.png`` |
| **5C** | **Yes** — proportional whole (default ``--coding-tsv``) | ``figures/fig5_merged_5c_epitope_sharing.png`` |
| **5C** | **Yes** — random-fragment merged coding (second run, ``--panels c``) | ``figures/supplementary/netmhc/coding_fragments_random_sample/figS5c_random_fragments_5c_epitope_sharing.png`` |
| **5D** | No (sig lnc SB only; same plot in both merged runs) | ``figures/fig5de_merged_whole_5d_sig_per_allele.png`` (canonical); duplicate under ``figures/supplementary/netmhc/coding_fragments_random_sample/figS5de_random_fragments_5d_sig_per_allele.png`` |
| **5E** | **Yes** — proportional-whole merged coding | ``figures/fig5de_merged_whole_5e_coding_per_allele.png`` |
| **5E** | **Yes** — random-fragment merged coding | ``figures/supplementary/netmhc/coding_fragments_random_sample/figS5de_random_fragments_5e_coding_per_allele.png`` |

**Wide XLS (legacy / supplement):** ``manuscript/plot_netmhc_epitopes_vs_hla_frequency.py`` (**5A**, IC50-from-BA on ``netmhcpan_sig_lnc.xls``); ``supplement/plot_figure5b_epitope_sharing_across_alleles.py`` / ``plot_figure5de_epitopes_per_allele.py`` for **5B–5E** — **IC50-from-BA only** on ``*.xls``; **not** the manuscript default for merged panels. Run via **`python generate_netmhc_supplement.py --include-wide-xls-fig5`**.

CSV companions for 5A–5C live alongside the PNGs under ``data/netmhc/figures/`` (not mirrored to ``figures/`` by default).

**Orchestrator:** ``generate_netmhc_figure_bundle.py`` — merged **5A–5C** (IEDB+NetMHC SB default) + second **5C** (random-fragment coding cohort; mirror under ``figures/supplementary/netmhc/``) + two merged **5D–5E** runs (whole → ``figures/``; fragments → supplementary folder) + **Fig 6** TTN split panels. **Supplement sensitivity:** ``generate_netmhc_fig5_fig6_supplement.py``. **Legacy wide XLS 5A–5E:** ``generate_netmhc_supplement.py --include-wide-xls-fig5`` (deprecated wrapper).

**Optional clean tree:** ``python supplement/regenerate_manuscript_netmhc_figures.py`` (see script ``--help`` for ``--clean``, purges, skips) → ``figures/manuscript_netmhc/``. ``--purge-repo-figures-netmhc`` removes top-level ``figures/fig5*`` / ``fig6*`` **and** files under ``figures/supplementary/netmhc/``.

**`figures/manuscript_netmhc/` layout** (typical top-level folders after a full ``--clean`` regenerate):

- ``fig5_wide_ic50_lt_150nm/{instances,unique}/`` — wide XLS cohort (**5A–5E**, IC50-from-BA only).
- ``fig5_merged/sb_{full,no_ic50,ic50_only}/{instances,unique}/`` — merged **5A–5C** per SB mode; includes two **5C** stems (proportional-whole coding + coding-control ``--panels c`` run).
- ``fig5_merged_de/sb_*/{instances,unique}/`` — merged **5D–5E** per-allele bars (**5D** sig-lnc duplicated; **5E** differs by coding cohort stem: ``*_proportional_whole`` vs ``*_random_fragments``).
- ``fig6/`` — TTN panels (e.g. ``netmhc_default/``, ``iedb_sb_*``).
- ``sensitivity/cohort/`` — ``netmhc_sb_sensitivity_robustness.py`` outputs; ``sensitivity/fig6_netmhc_sb_sweeps/`` — Fig 6 NetMHC-only sweeps.

Wide cohort scripts accept ``--no-repo-mirror`` so outputs stay under this tree without extra copies in repo-root ``figures/``. **Copy-paste commands** for orchestrators and one-off runs: **`docs/netmhc_figure_commands.md`**.

**Figure 5 — sensitivity & robustness (IEDB+NetMHC SB):**

- **`supplement/netmhc_sb_sensitivity_robustness.py`** — one-dimensional sweeps, leave-one-filter-out,
  main CSV + **`sb_threshold_sensitivity_robustness_fold_change_vs_baseline.csv`**.
  Default curves use **SB row instances**; ``--plot-metric unique`` restores unique-peptide y-axes.
  For a **repo-local supplement bundle** (recommended), run **`generate_netmhc_fig5_fig6_supplement.py`**
  → outputs under **`figures/supplementary/netmhc_fig5_fig6_supplement/fig5_merged_cohort_1d_sensitivity_loo/`**.
  Ad-hoc default when you run the script alone remains **`data/netmhc/figures/`** (``--out-dir``).

**Figure 5 — supplement (Cartesian SB grid on merged cohort TSVs):**

- **`supplement/plot_fig5_netmhc_sb_combination_grid.py`** — **full Cartesian product** of SB thresholds
  (immuno × processing × EL × IC50), fold-change table, 3×3 sharing grids, heatmap slice.
  **Orchestrated default** (``generate_netmhc_fig5_fig6_supplement.py``):
  **`figures/supplementary/netmhc_fig5_fig6_supplement/fig5_merged_cohort_cartesian_sb_grid/`**. When run standalone,
  default folder is **`data/netmhc/figures/fig5_netmhc_sb_combinations/`** (override with ``--out-dir``).

> **Naming note:** an earlier draft put combination-grid outputs under `fig6_sb_combinations/`.
> Per this catalog, **Figure 6 is reserved for TTN-AS1**; use the **`fig5_netmhc_sb_combinations`**
> path going forward.

---

## Figure 6 — TTN-AS1 (smPEP 108065): single-peptide allele / epitope coverage

**Scope:** one parent peptide (**TTN-AS1**, 79 aa; default sequence embedded unless
`--parent-fasta` is passed) with NetMHCpan wide output **`data/netmhc/netmhcpan_ttn_as1_108065.xls`**.

**Scripts:** canonical CLI entry **`manuscript/manuscript_figure6_ttn_as1.py`** (forwards to
`manuscript/plot_figure6_ttn_as1_allele_coverage.py`).

- **Default coverage metric:** **instances** — heatmap, histogram, and top coverage track use
  **total SB peptide×allele hits** overlapping each residue (`--coverage-output instances`, default).
  **Unique** epitope / allele counts are optional: **`--coverage-output unique`**, **`both`**, or
  **`--also-write-unique`** with default `instances` (writes a second `*_unique*` file set).
  Output paths insert **`_instances`** or **`_unique`** before the extension (unless the stem
  already ends with that tag), so split panels become e.g. `fig6_ttn_as1_split_instances_a.png`.
- With **`instances`**, the **third row overlay** (two curves on the same axes) plots **SB hit
  instances** vs **distinct alleles** per site (aligned with the instance-mode heatmap and top
  track). The **middle** row remains **unique SB 9-mers** per site (`epitope_cov`). In **`unique`**
  mode, all three rows use unique epitope / allele metrics consistently.
- Default combined PNG: **`figures/fig6_ttn_as1_allele_coverage.png`** (override with `-o`).
- **`--split-panels`** writes **5** panels (6A–6E style): sequence heatmap, coverage histogram,
  three coverage tracks, two sequence logos. `generate_netmhc_figure_bundle.py` uses
  **`figures/fig6_ttn_as1_split.png`** → `fig6_ttn_as1_split_instances_a.png` … `_e.png` by default.
  Pass **`--also-write-unique`** on the bundle to also write `fig6_ttn_as1_split_unique_*.png` in the
  same folder; **`generate_canonical_manuscript_figures.py`** instead writes **unique** split panels
  under **`figures/supplementary/figure6_ttn_as1/`** only.

**Gating (canonical TTN script):** default is **`--gating iedb_sb`**: merge NetMHC wide rows to
**`--iedb-csv`** on **`stable_key`**, then apply the **same cohort-style SB bundle** as merged Fig 5
(immuno / processing / EL / IC50 via ``netmhc_sb_core``; **IC50 uses the merged IEDB IC50 column when
present**, otherwise the BA-derived fallback inside ``sb_mask_spec``). Use **`--gating netmhc`** for
**NetMHC wide XLS only** (BA_rank or IC50-from-``BA_score``; no IEDB immunogenicity/processing gates).

**Supplement grids (not the main multi-panel figure):** NetMHC-only sweeps →
``supplement/plot_figure6_ttn_as1_sb_sensitivity.py``. **IEDB+NetMHC 1D + LOO** on the merged TTN table →
``supplement/netmhc_ttn_merged_iedb_sb_sensitivity_robustness.py``. **IEDB+NetMHC Cartesian** grid on the
same merge → ``supplement/plot_fig6_ttn_merged_iedb_sb_combination_grid.py``.

**SB definition when ``--gating netmhc``:** **`--sb-criterion ba_rank`** (default BA_rank ≤ 0.5 %) or
`ic50` from `BA_score`; optional **`--require-el-rank`**.

### Figure 6 — sensitivity & Cartesian grids (supplement)

**NetMHC wide XLS only (no IEDB):** `supplement/plot_figure6_ttn_as1_sb_sensitivity.py` — BA_rank /
IC50-from-BA / optional EL rows. **Orchestrated:** **`figures/supplementary/netmhc_fig5_fig6_supplement/fig6_ttn_wide_netmhc_sb_sweeps/`**;
standalone default: **`data/netmhc/figures/fig6_ttn_as1_sensitivity/`**. **`fig6_ttn_as1_sb_sensitivity_overview.png`**.

**IEDB + NetMHC merged (TTN stable_key) — 1D + LOO:** `supplement/netmhc_ttn_merged_iedb_sb_sensitivity_robustness.py` —
same **baseline → leave-one-out → 1D sweeps** structure as Fig 5’s ``netmhc_sb_sensitivity_robustness.py``.
**Orchestrated:** **`figures/supplementary/netmhc_fig5_fig6_supplement/fig6_ttn_merged_iedb_1d_sensitivity_loo/`**;
standalone default: **`data/netmhc/figures/fig6_ttn_merged_iedb_1d_sensitivity/`**.

**IEDB + NetMHC merged (TTN stable_key) — Cartesian grid:** `supplement/plot_fig6_ttn_merged_iedb_sb_combination_grid.py` —
same **Cartesian** immuno × proc × EL × IC50 exploration style as Fig 5’s
``plot_fig5_netmhc_sb_combination_grid.py``, on the single-peptide long table. **Orchestrated:**
**`figures/supplementary/netmhc_fig5_fig6_supplement/fig6_ttn_merged_iedb_cartesian_sb_grid/`**;
standalone default: **`data/netmhc/figures/fig6_ttn_merged_iedb_sb_combinations/`**. **PNG suite (fold vs manuscript baseline):**
`fig6_ttn_merged_iedb_cartesian_sweep_a_heatmaps_fold_unique.png`, `..._b_imm.png` … `..._e_ic50.png`, plus
`..._abcde_combined.png` (pass ``--no-suite-pngs`` for CSV-only).

**Scope:** the NetMHC-only script does **not** sweep IEDB immunogenicity or processing. The TTN
merged IEDB **1D + LOO** and **Cartesian** scripts do. Fig 5 merged-cohort **1D + LOO** (sig vs coding)
lives in ``netmhc_sb_sensitivity_robustness.py`` (two-table design).

---

## Quick reference

| Goal | Entry point |
|------|----------------|
| Regenerate **everything** (catalog → NetMHC bundle → Fig 5–6 supplement → canonical extras → publication export) | `regenerate_all_figures.py` |
| Regenerate Fig 2–4A (SmProt catalog) | `generate_catalog_figures.py` |
| Regenerate Fig 5–6 NetMHC **canonical** | `generate_netmhc_figure_bundle.py` |
| Fig 5–6 NetMHC **supplement folder** (Fig 5 1D+LOO, Fig 5 Cartesian, Fig 6 NetMHC sweeps, Fig 6 merged IEDB 1D+LOO, Fig 6 merged IEDB Cartesian) | `generate_netmhc_fig5_fig6_supplement.py` |
| NetMHC **supplement** sensitivity bundle | `generate_netmhc_fig5_fig6_supplement.py` |
| NetMHC **legacy** wide-XLS 5A–5E (optional) | `generate_netmhc_supplement.py` (deprecated) |
| **Orchestrator roadmap** (canonical vs supplement vs NetMHC) | **`docs/figure_generation_overview.md`** |
| **Command cheat sheet** (orchestrators + one-off CLIs) | **`docs/netmhc_figure_commands.md`** |
| **Supplementary figure tree** (`figures/supplementary/`) | **`figures/supplementary/README.md`** |
| Merge NetMHC wide XLS + IEDB CSV | `scripts/merge_netmhcpan_xls_with_iedb.py` |
| SB filter logic shared by scripts | `scripts/netmhc_sb_core.py` |
| Full clean matrix under `figures/manuscript_netmhc/` | `supplement/regenerate_manuscript_netmhc_figures.py` |

Project skill for NetMHC manuscript figures: **`.cursor/skills/netmhc-manuscript-figures/SKILL.md`**.
