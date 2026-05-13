# NetMHC and catalog figure commands

Single reference for **orchestrators**, **manuscript figure scripts**, **data prep**, and **sensitivity** jobs. Run everything from the repository root unless noted.

| Doc | Role |
|-----|------|
| This file | Copy-paste commands for Figures **2–4A**, **5–6**, merges, sensitivity |
| `docs/figure_catalog.md` | What each panel shows, inputs, limitations |
| `data/netmhc/README_netmhc.md` | NetMHC file layout, merge pipeline, column semantics |

---

## Orchestrators (batch)

| Goal | Command |
|------|---------|
| **Figures 2–4A** (SmProt / TCGA catalog): Fig **1B** + `figures/tcga_matrix/`, `figures/all_smprot_filtered/`, shared `figures/` for 3C–3D + 4A | `python generate_catalog_figures.py` — add `--strict` to fail on missing optional inputs; `--only tcga_matrix` or `all_smprot_filtered` |
| **Figures 5–6** (NetMHC **canonical**: merged **5A–5E**, TTN Fig 6 coverage) | `python generate_netmhc_figure_bundle.py` — **`--strict`** for non-zero on failure; **`--skip-iedb-pipeline`** skips merged cohort steps |
| **Figures 5–6** (NetMHC **supplement**: legacy wide **5A–5E**, cohort sensitivity, SB combination grid, Fig 6 SB sweeps) | `python generate_netmhc_supplement.py` — **`--include-wide-xls-fig5`** for wide-XLS **5A–5E**; **`--skip-sensitivity`** to skip sensitivity / combo grid / Fig 6 sweeps |
| **Clean matrix: Fig 5–6 + all SB modes + instances/unique + sensitivity** | `python supplement/regenerate_manuscript_netmhc_figures.py --clean --purge-repo-figures-netmhc --purge-data-netmhc-figures` — see script `--help` (`--dry-run`, skips, optional `--with-combination-grid`) |

### `figures/manuscript_netmhc/` layout (after regenerate script)

| Path | Contents |
|------|-----------|
| `fig5_wide_ic50_lt_150nm/{instances,unique}/` | Wide XLS 5A–5E (IC50-from-BA gate) |
| `fig5_merged/sb_{full,no_ic50,ic50_only}/{instances,unique}/` | Merged TSV 5A–5C (`full` = IEDB+EL+IC50; `no_ic50` = IEDB+EL without binding gate; `ic50_only` = binding/IC50 only). Two **5C** products per folder: proportional-whole coding + ``*_coding_control`` stem from ``--panels c`` |
| `fig5_merged_de/sb_*/{instances,unique}/` | Merged 5D–5E per-allele bars: two stems per folder, ``*_proportional_whole`` and ``*_random_fragments`` (same sig-lnc **5D** in both runs; **5E** differs by coding cohort) |
| `fig6/netmhc_default/` | TTN split panels, instances + unique companion PNGs |
| `fig6/iedb_sb_*/` | TTN with `--gating iedb_sb` per `--sb-mode` |
| `sensitivity/cohort/sb_*_{instances,unique}/` | `netmhc_sb_sensitivity_robustness.py` products |
| `sensitivity/fig6_netmhc_sb_sweeps/` | `plot_figure6_ttn_as1_sb_sensitivity.py` |

Wide Fig 5 scripts support **`--no-repo-mirror`** so outputs stay only under this tree (no extra copies in repo-root `figures/`).

---

## Catalog Figures 2–4A (individual scripts)

Same steps as `generate_catalog_figures.py`; override output roots with each script’s `--help`.

| Figure | Script | Typical invocation |
|--------|--------|--------------------|
| **1B** | `manuscript/plot_figure1b_tsne_stage_lncrna.py` | `python manuscript/plot_figure1b_tsne_stage_lncrna.py` — needs `data/primary_exp_stage_lncRNA.csv` and **openTSNE** |
| **2** | `manuscript/plot_tr_de_peptide_fractions_by_transition.py` | `python manuscript/plot_tr_de_peptide_fractions_by_transition.py --peptide-gene-set tcga_matrix` (or `all_smprot_filtered`) |
| **3A** | `manuscript/plot_aa_frequency_tcga_vs_proteome.py` | `python manuscript/plot_aa_frequency_tcga_vs_proteome.py --peptide-set tcga_matrix` |
| **3B** | `manuscript/plot_dipeptide_volcano_lnc_vs_proteome.py` | `python manuscript/plot_dipeptide_volcano_lnc_vs_proteome.py --peptide-set tcga_matrix` |
| **3C–3D** | `manuscript/plot_figure3cd_dipeptide_log2fc_heatmaps.py` | `python manuscript/plot_figure3cd_dipeptide_log2fc_heatmaps.py` |
| **4A** | `manuscript/plot_figure4a_tis_vs_ribo_tr_mps.py` | `python manuscript/plot_figure4a_tis_vs_ribo_tr_mps.py` |

---

## Figure 5 — cohort NetMHCpan (wide XLS; **supplement** / legacy)

Default PNG/CSV under `data/netmhc/figures/` (some scripts also mirror copies into repo-root `figures/`; see each script). **Manuscript Figure 5** uses the **merged** pipeline (`plot_fig5abc_netmhc_sb_triple.py` / `plot_fig5de_merged_iedb_sb_per_allele.py`), not these wide-XLS cohort plots.

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
| **5A–5C** merged | `manuscript/plot_fig5abc_netmhc_sb_triple.py` | `--out-dir data/netmhc/figures` (stem `fig5abc_sb_immuno_proc_el_ic50`); second coding cohort: `--panels c --coding-tsv data/netmhc/netmhcpan_coding_control_with_iedb.tsv --output-stem fig5abc_sb_immuno_proc_el_ic50_coding_control` |
| **5D–5E** merged | `manuscript/plot_fig5de_merged_iedb_sb_per_allele.py` | Default stem ``fig5de_merged_iedb_sb_proportional_whole``; second cohort: ``--coding-tsv data/netmhc/netmhcpan_coding_control_with_iedb.tsv --output-stem fig5de_merged_iedb_sb_random_fragments``; ``--out-dir`` default `data/netmhc/figures` |
| **Sensitivity / robustness** | `supplement/netmhc_sb_sensitivity_robustness.py` | Under `data/netmhc/figures/` (see `--help`) |
| **Combination grid** | `supplement/plot_fig5_netmhc_sb_combination_grid.py` | `data/netmhc/figures/fig5_netmhc_sb_combinations/` |

---

## Figure 6 — TTN-AS1 (smPEP 108065)

**Canonical manuscript entry point (one script per figure):** `manuscript/manuscript_figure6_ttn_as1.py` — forwards argv to `manuscript/plot_figure6_ttn_as1_allele_coverage.py`.

| Task | Command |
|------|---------|
| Default combined PNG (instances coverage) | `python manuscript/manuscript_figure6_ttn_as1.py` |
| Same, explicit path | `python manuscript/plot_figure6_ttn_as1_allele_coverage.py -o figures/fig6_ttn_as1_allele_coverage.png` |
| Split panels (6A–6E) + unique companion | `python manuscript/plot_figure6_ttn_as1_allele_coverage.py --split-panels --also-write-unique -o figures/fig6_ttn_as1_split.png` |
| IEDB-gated SB | `python manuscript/plot_figure6_ttn_as1_allele_coverage.py --gating iedb_sb --iedb-csv <path> --iedb-parent-input-seq-id <id>` (see `--help`; on PowerShell avoid `$PID` as a variable name) |
| **TTN synthetic IEDB companion CSV** | `python supplement/build_ttn_iedb_companion_csv.py` |
| **Fig 6 sensitivity** (NetMHC SB sweeps only) | `python supplement/plot_figure6_ttn_as1_sb_sensitivity.py` → `data/netmhc/figures/fig6_ttn_as1_sensitivity/` |
| **Fig 6 — local NetMHCpan (FASTA → XLS)** | WSL: `bash data/netmhc/run_netmhcpan_ttn_as1_108065.sh` → `data/netmhc/netmhcpan_ttn_as1_108065.xls`. Flags / install: **`data/netmhc/README_netmhc.md`** §5; one-page summary: **`docs/figure6_ttn_as1_parameters.md`** (*Local NetMHCpan*). Same layout under parent **`UNDEFINED/data/netmhc/`** if you run from the full tree. |
| **Fig 6 from IEDB Tools API (NetMHCpan 4.1 BA+EL)** | `python supplement/fetch_ttn_mhci_iedb_api_netmhc41.py` — **POST fields / fetch CLI:** `docs/iedb_tools_api.md`. Default output dir `data/netmhc/predictions/ttn_as1_smpep108065_iedb_api_netmhc41/` (see `README.md` there). On Windows add `--insecure` if TLS fails. Then plot with `--netmhc-xls …/netmhcpan_ttn_as1_iedb_api_netmhc41.xls` (same Fig 6 script; see **`docs/figure6_ttn_as1_parameters.md`**) |

## Tests

This repo may not ship a unified test suite for plotting. If `pytest` is configured, from the repo root run:

`python -m pytest`

If that reports “no tests collected”, rely on the orchestrators above and spot-check outputs under `figures/` and `data/netmhc/figures/`.
