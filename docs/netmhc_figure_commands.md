# NetMHC and catalog figure commands

**Copy-paste commands** for orchestrators and one-off figure runs (repo root). For **what each figure shows**, inputs, and **`figures/manuscript_netmhc/`** folder layout, see **`docs/figure_catalog.md`**. **Which driver to run when:** **`docs/figure_generation_overview.md`**. For NetMHC **files, merge, WSL**, see **`data/netmhc/README_netmhc.md`**.

---

## Orchestrators (batch)

| Goal | Command |
|------|---------|
| **Figures 2–4A** (SmProt / TCGA catalog): Fig **1B** + `figures/supplementary/tcga_matrix/`, `figures/supplementary/all_smprot_filtered/`, shared `figures/` for 3C (TCGA) + 3D + 4A | `python generate_catalog_figures.py` — add `--strict` to fail on missing optional inputs; `--only tcga_matrix` or `all_smprot_filtered` |
| **Main-text bundle (1B sklearn + OpenTSNE supplement, tcga-matrix 2–4A, NetMHC core; optional Fig 6 unique)** | `python generate_canonical_manuscript_figures.py` — **`--strict`**; **`--fig1b-embedding opentsne4`** uses OpenTSNE for canonical 1B; **`--write-fig6-unique-supplement`** emits Fig 6 unique under `figures/supplementary/figure6_ttn_as1/`; **`--skip-opentsne-supplement`** / **`--skip-netmhc`** as needed |
| **Figures 5–6** (NetMHC **canonical**: merged **5A–5E**, TTN Fig 6 coverage) | `python generate_netmhc_figure_bundle.py` — **`--strict`** for non-zero on failure; **`--skip-iedb-pipeline`** skips merged cohort steps; **`--canonical-main-text-only`** = proportional-whole coding cohort only + no random-fragment mirrors; **Fig 6** defaults to **instances** split panels only — add **`--also-write-unique`** to also emit `*_unique_*` next to `*_instances_*` under `figures/` |
| **Figures 5–6** (NetMHC **supplement bundle** under `figures/supplementary/netmhc_fig5_fig6_supplement/`) | `python generate_netmhc_fig5_fig6_supplement.py` — **`--strict`** |
| **Figures 5–6** (NetMHC **supplement** sensitivity bundle) | `python generate_netmhc_fig5_fig6_supplement.py` — **`--strict`** |
| **Figures 5–6** (legacy wide-XLS **5A–5E** + optional § bundle) | `python generate_netmhc_supplement.py` — **deprecated**; use only with **`--include-wide-xls-fig5`** for IC50-from-BA wide `*.xls` panels |
| **Clean matrix: Fig 5–6 + all SB modes + instances/unique + sensitivity** | `python supplement/regenerate_manuscript_netmhc_figures.py --clean --purge-repo-figures-netmhc --purge-data-netmhc-figures` — see script `--help` (`--dry-run`, skips, optional `--with-combination-grid`) |

---

## Catalog Figures 2–4A (individual scripts)

Same steps as `generate_catalog_figures.py`; override output roots with each script’s `--help`.

| Figure | Script | Typical invocation |
|--------|--------|--------------------|
| **1B** | `manuscript/plot_figure1b_tsne_stage_lncrna.py` | `python manuscript/plot_figure1b_tsne_stage_lncrna.py` — needs `data/primary_exp_stage_lncRNA.csv`; default embedding is sklearn (see `--embedding`; **openTSNE** only for `--embedding opentsne4`) |
| **2** | `manuscript/plot_tr_de_peptide_fractions_by_transition.py` | `python manuscript/plot_tr_de_peptide_fractions_by_transition.py --peptide-gene-set tcga_matrix` (or `all_smprot_filtered`) |
| **3A** | `manuscript/plot_aa_frequency_tcga_vs_proteome.py` | `python manuscript/plot_aa_frequency_tcga_vs_proteome.py --peptide-set tcga_matrix` |
| **3B** | `manuscript/plot_dipeptide_volcano_lnc_vs_proteome.py` | `python manuscript/plot_dipeptide_volcano_lnc_vs_proteome.py --peptide-set tcga_matrix` |
| **3C–3D** | `manuscript/plot_figure3cd_dipeptide_log2fc_heatmaps.py` | `python manuscript/plot_figure3cd_dipeptide_log2fc_heatmaps.py` — add `--only-tcga-matrix-3c` for TCGA-matrix 3C only (no all-filtered 3C) |
| **4A** | `manuscript/plot_figure4a_tis_vs_ribo_tr_mps.py` | `python manuscript/plot_figure4a_tis_vs_ribo_tr_mps.py` |

---

## Figure 5 — cohort NetMHCpan (wide XLS; **supplement** / legacy)

Default PNG/CSV under `data/netmhc/figures/` (some scripts also mirror copies into repo-root `figures/` or `figures/supplementary/...`; see each script and `figures/supplementary/README.md`). **Manuscript Figure 5** uses the **merged** pipeline (`plot_fig5abc_netmhc_sb_triple.py` / `plot_fig5de_merged_iedb_sb_per_allele.py`), not these wide-XLS cohort plots.

| Panel | Script | Notes |
|-------|--------|--------|
| **5A** | `manuscript/plot_netmhc_epitopes_vs_hla_frequency.py` | **Supplement** (orchestrated by `generate_netmhc_supplement.py --include-wide-xls-fig5`); IC50-from-BA on wide `*.xls`; default freq table `data/netmhc/hla_european27_allele_frequencies.csv` |
| **5B / 5C** | `supplement/plot_figure5b_epitope_sharing_across_alleles.py` | IC50-from-BA on wide XLS only |
| **5D / 5E** | `supplement/plot_figure5de_epitopes_per_allele.py` | IC50-from-BA on wide XLS |

---

## Figure 5 — merged `*_with_iedb.tsv` (**manuscript** default)

Requires merged long TSVs from `scripts/merge_netmhcpan_xls_with_iedb.py` (see `data/netmhc/README_netmhc.md`). **Fig 5A** in the merged bundle uses the same allele-frequency reference as wide 5A, but **SB counts come from merged rows** (IEDB + EL + IC50 gates), not from wide XLS alone.

| Product | Script | Default output dir |
|---------|--------|--------------------|
| **5A–5C** merged | `manuscript/plot_fig5abc_netmhc_sb_triple.py` | `--out-dir data/netmhc/figures` (default stem `fig5_merged` → `fig5_merged_5a_*` … `5c_*`); random-fragment **5C** only: `--panels c --coding-tsv data/netmhc/netmhcpan_coding_control_with_iedb.tsv --output-stem figS5c_random_fragments --repo-mirror-dir figures/supplementary/netmhc/coding_fragments_random_sample` |
| **5D–5E** merged | `manuscript/plot_fig5de_merged_iedb_sb_per_allele.py` | Default stem `fig5de_merged_whole` (mirrors to `figures/`); random-fragment cohort: `--coding-tsv data/netmhc/netmhcpan_coding_control_with_iedb.tsv --output-stem figS5de_random_fragments --repo-mirror-dir figures/supplementary/netmhc/coding_fragments_random_sample`; `--out-dir` default `data/netmhc/figures` |
| **Sensitivity / robustness** | `supplement/netmhc_sb_sensitivity_robustness.py` | Default `data/netmhc/figures/`; **`generate_netmhc_fig5_fig6_supplement.py`** → `figures/supplementary/netmhc_fig5_fig6_supplement/fig5_merged_cohort_1d_sensitivity_loo/` |
| **Fig 5 Cartesian SB grid** | `supplement/plot_fig5_netmhc_sb_combination_grid.py` | Default `data/netmhc/figures/fig5_netmhc_sb_combinations/`; orchestrator → `figures/supplementary/netmhc_fig5_fig6_supplement/fig5_merged_cohort_cartesian_sb_grid/` |
| **Fig 6 merged IEDB — 1D + LOO (TTN)** | `supplement/netmhc_ttn_merged_iedb_sb_sensitivity_robustness.py` | Default `data/netmhc/figures/fig6_ttn_merged_iedb_1d_sensitivity/`; orchestrator → `figures/supplementary/netmhc_fig5_fig6_supplement/fig6_ttn_merged_iedb_1d_sensitivity_loo/` |
| **Fig 6 Cartesian IEDB+NetMHC grid (TTN)** | `supplement/plot_fig6_ttn_merged_iedb_sb_combination_grid.py` | Default `data/netmhc/figures/fig6_ttn_merged_iedb_sb_combinations/`; orchestrator → `figures/supplementary/netmhc_fig5_fig6_supplement/fig6_ttn_merged_iedb_cartesian_sb_grid/` |

---

## Supplement — PCA (same expression matrix as Fig 1B)

| Product | Script | Default output |
|---------|--------|----------------|
| **Fig S — PCA panels** | `supplement/plot_supplement_pca_stage_samples.py` | `figures/supplementary/pca/figS_pca_stage_lncrna_samples_*.png` (override with `--out-dir`) |

`export_publication_figures.py` runs this step after Fig 1B unless you pass `--skip-supplement-pca`.

---

## Figure 6 — TTN-AS1 (smPEP 108065)

**Canonical manuscript entry point (one script per figure):** `manuscript/manuscript_figure6_ttn_as1.py` — forwards argv to `manuscript/plot_figure6_ttn_as1_allele_coverage.py`.

| Task | Command |
|------|---------|
| Default combined PNG (instances coverage) | `python manuscript/manuscript_figure6_ttn_as1.py` |
| Same, explicit path | `python manuscript/plot_figure6_ttn_as1_allele_coverage.py -o figures/fig6_ttn_as1_allele_coverage.png` |
| Split panels (6A–6E) | Instances: ``python manuscript/plot_figure6_ttn_as1_allele_coverage.py --split-panels -o figures/fig6_ttn_as1_split.png``. Unique (supplementary tree): add ``--coverage-output unique`` and ``-o figures/supplementary/figure6_ttn_as1/fig6_ttn_as1_split.png``. Optional: ``--also-write-unique`` with instances + ``-o figures/...`` writes both `*_instances_*` and `*_unique_*` in the **same** directory. |
| IEDB-gated SB | `python manuscript/plot_figure6_ttn_as1_allele_coverage.py --gating iedb_sb --iedb-csv <path> --iedb-parent-input-seq-id <id>` (see `--help`; on PowerShell avoid `$PID` as a variable name) |
| **TTN synthetic IEDB companion CSV** | `python supplement/build_ttn_iedb_companion_csv.py` |
| **Fig 6 NetMHC-only sensitivity** | `python supplement/plot_figure6_ttn_as1_sb_sensitivity.py` → default `data/netmhc/figures/fig6_ttn_as1_sensitivity/`; **`python generate_netmhc_fig5_fig6_supplement.py`** → `figures/supplementary/netmhc_fig5_fig6_supplement/fig6_ttn_wide_netmhc_sb_sweeps/` |
| **Fig 6 — local NetMHCpan (FASTA → XLS)** | WSL: `bash data/netmhc/run_netmhcpan_ttn_as1_108065.sh` → `data/netmhc/netmhcpan_ttn_as1_108065.xls`. Install, paths, flags: **`data/netmhc/README_netmhc.md`**; TTN one-pager: **`docs/figure6_ttn_as1_parameters.md`**. |
| **Fig 6 from IEDB Tools API (NetMHCpan 4.1 BA+EL)** | `python supplement/fetch_ttn_mhci_iedb_api_netmhc41.py` — request/CLI details: **`docs/iedb_tools_api.md`**. Then plot with `--netmhc-xls` pointing at the generated XLS (**`docs/figure6_ttn_as1_parameters.md`**). |
