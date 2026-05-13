# Figure catalog

Structured notes for selected manuscript-style figures: what they show, how they are generated, inputs and outputs, inferential procedures (if any), and limitations.

**Output root:** all scripts below write under **`figures/<peptide_mode>/`** at the repository root (`UNDEFINED/figures/`), where **`peptide_mode`** is **`tcga_matrix`** or **`all_smprot_filtered`**. Regenerate both modes in one go: **`python generate_catalog_figures.py`**.

**Figure 2** — peptide-fraction bars (TCGA limma–z). **Figure 3** — composition vs `known_proteins.fasta` (**3A** 1-mer, **3B** volcano, **3C** log2FC dipeptide heatmaps in **separate** TCGA-matrix and all-filtered files, **3D** Tr log2FC heatmap). **Figure 4A** — TIS vs Ribo-seq p-values for Tr MPs (`plot_figure4a_tis_vs_ribo_tr_mps.py`). **3A–3B** default/alternate use **`figures/<mode>/`**; **3C–3D** and **4A** live under **`figures/`** at repo root.

---

## Figure 2 — Peptide-producing fraction among DE ∩ high-|z| genes, by cancer

**Paths:** `figures/<peptide_gene_set>/peptide_fraction/*.png` and `figures/<peptide_gene_set>/tr_de_peptide_fraction_by_cancer.csv`

| `peptide_gene_set` | Peptide TSV |
|--------------------|-------------|
| `tcga_matrix` (default) | `data/smprot_filtered_tcga_expr_genes.tsv` |
| `all_smprot_filtered` | `data/smprot_filtered.tsv` |

**Script:** `plot_tr_de_peptide_fractions_by_transition.py` — optional **`--figures-dir`** to override the `figures/` parent.

---

## Figure 3 — lncRNA peptide composition vs **`data/known_proteins.fasta`**

### Figure 3A — Pooled 1-mer frequency

**Paths:** `figures/<peptide_set>/aa_frequency_*_vs_known_proteins.{png,csv,txt}`

**Script:** `plot_aa_frequency_tcga_vs_proteome.py` — default **`tcga_matrix`**; **`--peptide-set all_smprot_filtered`** needs **`data/smprot_all_filtered_peptides.faa`**.

### Figure 3B — Dipeptide volcano

**Paths:** `figures/<peptide_set>/dipeptide_volcano_lnc_vs_proteome*.{png,csv,txt}`

**Script:** `plot_dipeptide_volcano_lnc_vs_proteome.py`

### Figure 3C — Dipeptide log2FC heatmaps (TCGA-matrix and all-filtered, **separate files**)

**Paths:**

- `figures/fig3c_dipeptide_log2fc_tcga_matrix_vs_proteome.png` — always (TCGA-matrix FASTA).
- `figures/fig3c_dipeptide_log2fc_all_smprot_filtered_vs_proteome.png` — when
  `data/smprot_all_filtered_peptides.faa` exists.

Each file has its own symmetric color scale (99th percentile of |log2 ratio|, floor 0.5).
**Script:** `plot_figure3cd_dipeptide_log2fc_heatmaps.py` — same +0.5 / 400-cell smoothing and
`log2(peptide/proteome)` as `plot_dipeptide_mp_figure2.py`.

### Figure 3D — Dipeptide log2FC: Tr-lncRNA MPs (TCGA-matrix FASTA) vs proteome

**Path:** `figures/fig3d_dipeptide_log2fc_tr_lncrna_tcga_vs_proteome.png`

**Script:** `plot_figure3cd_dipeptide_log2fc_heatmaps.py` (second output). Restricts
`data/smprot_tcga_filtered_peptides.faa` to **canonical** Tr gene symbols
(`canonical_significant_lncRNAs.txt` or `limma_z_intersection_genes.txt`).

**Legacy multi-panel figure:** `plot_dipeptide_mp_figure2.py` → `tr_lncrna_output/figures/fig2_dipeptide_mp_composition.png`

### Figure 4A — TIS vs Ribo-seq p-values (TCGA-matrix filtered lncRNA MPs)

**Path:** `figures/fig4a_tr_lncrna_mp_tis_vs_riboseq_pvalues.png` (filename legacy; plot is **not** Tr-only.)

**Script:** `plot_figure4a_tis_vs_ribo_tr_mps.py` — **all** rows from `data/smprot_filtered_tcga_expr_genes.tsv`
(TCGA expression–matrix gene filter only; **no** canonical Tr gene filter) with both **TISPvalue** and **RiboPvalue** ≤ 0.05; log-scaled
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
| 2 | Peptide fraction by cancer × transition | `plot_tr_de_peptide_fractions_by_transition.py` |
| 3A | 1-mer AA vs proteome | `plot_aa_frequency_tcga_vs_proteome.py` |
| 3B | Dipeptide volcano vs proteome | `plot_dipeptide_volcano_lnc_vs_proteome.py` |
| 3C–3D | Dipeptide log2FC heatmaps (3C split files / 3D Tr) | `plot_figure3cd_dipeptide_log2fc_heatmaps.py` |
| 4A | TIS vs Ribo-seq p scatter (TCGA filtered MPs) | `plot_figure4a_tis_vs_ribo_tr_mps.py` |
| — | Both modes + 3C–3D + 4A (orchestrator) | `generate_catalog_figures.py` |

---

## Figure 5 — NetMHCpan cohort (manuscript: merged IEDB+NetMHC, **SB default**)

**Manuscript stance:** Figure 5 in this catalog refers to the **merged** ``*_with_iedb.tsv`` pipeline
with **``sb_mode=full``** (default) and thresholds shared with Figure 6 IEDB gating:
**``FIG5_IEDB_*``** in ``scripts/netmhc_sb_core.py`` (immunogenicity **> 0.1**, processing **> 1.5**,
EL **< 1%** by default, IC50 **< 150 nM** using the IEDB BA-IC50 column when present, else local
``BA_score``). Merge inputs with ``scripts/merge_netmhcpan_xls_with_iedb.py``.

**Sig lnc cohort (all panels):** ``data/netmhc/netmhcpan_sig_lnc_with_iedb.tsv``.

**Coding cohorts (merged TSVs; two **5C** outputs):**

| Cohort | Merged TSV | Manuscript **5C** output stem (default ``--out-dir``) |
|--------|------------|--------------------------------------------------------|
| **Proportional whole** (length-matched parents) | ``data/netmhc/netmhcpan_coding_proportional_whole_with_iedb.tsv`` | ``fig5abc_sb_immuno_proc_el_ic50`` (same run as 5A–5B) |
| **Length-matched coding control** | ``data/netmhc/netmhcpan_coding_control_with_iedb.tsv`` | ``fig5abc_sb_immuno_proc_el_ic50_coding_control`` (``scripts/plot_fig5abc_netmhc_sb_triple.py --panels c``) |

**5B** is **significant lncRNA only** (one canonical stem per run); the second coding cohort does not re-emit 5B.

**Scripts:** ``scripts/plot_fig5abc_netmhc_sb_triple.py`` (5A–5C) and ``scripts/plot_fig5de_merged_iedb_sb_per_allele.py`` (5D–5E; run twice with ``--coding-tsv`` + ``--output-stem`` for proportional-whole vs random-fragment merged coding). Default plot metrics use **SB row instances** (``--fig5a-y-metric`` / ``--sharing-y-metric`` / ``--count-metric`` can switch to **unique**).

**Repo-root mirrors:** PNGs are written under ``data/netmhc/figures/`` and copied into **``figures/``** (same basename) unless ``--no-repo-mirror`` is passed.

| Panel | Uses coding cohort? | ``figures/`` PNG (manuscript merged **full** SB; default stems) |
|-------|---------------------|------------------------------------------------------------------|
| **5A** | No (sig lnc SB only) | ``fig5abc_sb_immuno_proc_el_ic50_5a_epitopes_vs_allele_frequency.png`` |
| **5B** | No | ``fig5abc_sb_immuno_proc_el_ic50_5b_epitope_sharing.png`` |
| **5C** | **Yes** — proportional whole (default ``--coding-tsv``) | ``fig5abc_sb_immuno_proc_el_ic50_5c_epitope_sharing.png`` |
| **5C** | **Yes** — coding control merged TSV (second run, ``--panels c``) | ``fig5abc_sb_immuno_proc_el_ic50_coding_control_5c_epitope_sharing.png`` |
| **5D** | No (sig lnc SB only; same plot in both merged runs) | ``fig5de_merged_iedb_sb_proportional_whole_5d_sig_per_allele.png`` (canonical; duplicate also under ``fig5de_merged_iedb_sb_random_fragments_5d_*``) |
| **5E** | **Yes** — proportional-whole merged coding | ``fig5de_merged_iedb_sb_proportional_whole_5e_coding_per_allele.png`` |
| **5E** | **Yes** — random-fragment merged coding | ``fig5de_merged_iedb_sb_random_fragments_5e_coding_per_allele.png`` |

**Wide XLS (legacy / supplement):** ``plot_figure5b_epitope_sharing_across_alleles.py`` cohort **5B / 5C** use **IC50-from-BA only** on ``*.xls``; they are **not** the manuscript default for those panels. ``generate_netmhc_figure_bundle.py`` runs them only with **`--include-wide-xls-fig5``** (also restores random-fragment wide **5C / 5E** companion paths under ``data/netmhc/figures/``).

CSV companions for 5A–5C live alongside the PNGs under ``data/netmhc/figures/`` (not mirrored to ``figures/`` by default).

**Orchestrator:** ``generate_netmhc_figure_bundle.py`` — default: wide **5A** only (IC50-from-BA); **merged** cohort **5B–5E** (canonical 5A–5C + second **5C**; **two** merged **5D–5E** stems for proportional-whole vs random-fragment coding); sensitivity / combination-grid helpers. Optional wide **5B–5E** via **`--include-wide-xls-fig5`**.

**Figure 5 — sensitivity & robustness (IEDB+NetMHC SB):**

- **`scripts/netmhc_sb_sensitivity_robustness.py`** — one-dimensional sweeps, leave-one-filter-out,
  main CSV + **`sb_threshold_sensitivity_robustness_fold_change_vs_baseline.csv`**.
  Default curves use **SB row instances**; ``--plot-metric unique`` restores unique-peptide y-axes.

**Figure 5 — supplement (combination grid, cohort / not Fig 6):**

- **`scripts/plot_fig5_netmhc_sb_combination_grid.py`** — Cartesian grid of SB thresholds,
  fold-change table, 3×3 sharing grids, heatmap slice. Default folder:
  **`data/netmhc/figures/fig5_netmhc_sb_combinations/`**.

> **Naming note:** an earlier draft put combination-grid outputs under `fig6_sb_combinations/`.
> Per this catalog, **Figure 6 is reserved for TTN-AS1**; use the **`fig5_netmhc_sb_combinations`**
> path going forward.

**Optional clean tree:** ``python scripts/regenerate_manuscript_netmhc_figures.py`` (see ``docs/netmhc_figure_commands.md``) → ``figures/manuscript_netmhc/``. **Wide XLS–only** cohort scripts (IC50-from-BA on ``*.xls``) remain documented in ``docs/netmhc_figure_commands.md`` for diagnostics; they are not the manuscript SB default in this table.

---

## Figure 6 — TTN-AS1 (smPEP 108065): single-peptide allele / epitope coverage

**Scope:** one parent peptide (**TTN-AS1**, 79 aa; default sequence embedded unless
`--parent-fasta` is passed) with NetMHCpan wide output **`data/netmhc/netmhcpan_ttn_as1_108065.xls`**.

**Scripts:** canonical CLI entry **`scripts/manuscript_figure6_ttn_as1.py`** (forwards to
`plot_figure6_ttn_as1_allele_coverage.py` at repo root).

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
  **`figures/fig6_ttn_as1_split.png`** → `fig6_ttn_as1_split_instances_a.png` … `_e.png` (and with
  **`--also-write-unique`**, `fig6_ttn_as1_split_unique_*.png`).

**Gating:** **`--gating netmhc`** (default) uses only the wide XLS. **`--gating iedb_sb`** merges
NetMHC rows to **`--iedb-csv`** on **`stable_key`** (with **`--iedb-parent-input-seq-id`**
matching IEDB `input_seq_id`) and applies the **same one-pass SB bundle** as merged Fig 5
(``--sb-mode`` + thresholds; **default IC50 cap 150 nM**). This is **not** a multi-threshold
grid on the figure itself—use ``plot_figure6_ttn_as1_sb_sensitivity.py`` for NetMHC-only sweeps.

**SB definition (netmhc):** **`--sb-criterion ba_rank`** (default BA_rank ≤ 0.5 %) or
`ic50` from `BA_score`; optional **`--require-el-rank`**.

### Figure 6 — sensitivity (TTN-AS1 SB sweeps)

**Script:** `plot_figure6_ttn_as1_sb_sensitivity.py`

- Tables + fold-change vs default (BA_rank 0.5 %, no EL requirement) under
  **`data/netmhc/figures/fig6_ttn_as1_sensitivity/`**.
- **`fig6_ttn_as1_sb_sensitivity_overview.png`** — small multi-panel diagnostic.

**Scope:** sensitivity mirrors the **NetMHCpan XLS** SB definition only (BA_rank / IC50 from
`BA_score` / optional EL). It does **not** sweep IEDB immunogenicity, processing, or merged
IEDB IC50—those belong to **Figure 5** merged-TSV scripts.

---

## Agent quick reference

| Goal | Entry point |
|------|----------------|
| Regenerate Fig 2–4A (SmProt catalog) | `generate_catalog_figures.py` |
| Regenerate Fig 5–6 NetMHC bundle | `generate_netmhc_figure_bundle.py` |
| **Command cheat sheet (all figures / tests)** | **`docs/netmhc_figure_commands.md`** (NetMHC + catalog); SmProt catalog details also in sections above |
| Merge NetMHC wide XLS + IEDB CSV | `scripts/merge_netmhcpan_xls_with_iedb.py` |
| SB filter logic shared by scripts | `scripts/netmhc_sb_core.py` |

Project skill for agents: **`.cursor/skills/netmhc-manuscript-figures/SKILL.md`**.
