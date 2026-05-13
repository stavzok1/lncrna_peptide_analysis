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
| **Figures 2–4A** (SmProt / TCGA catalog): `figures/tcga_matrix/`, `figures/all_smprot_filtered/`, shared `figures/` for 3C–3D + 4A | `python generate_catalog_figures.py` — add `--strict` to fail on missing optional inputs; `--only tcga_matrix` or `all_smprot_filtered` |
| **Figures 5–6** (NetMHC: wide **5A** only by default; **merged** cohort **5B–5E** including two merged **5D–5E** stems; optional legacy wide **5B–5E** + TTN Fig 6) | `python generate_netmhc_figure_bundle.py` — **`--include-wide-xls-fig5`** adds legacy wide-XLS epitope-sharing (5B/5C), wide **5D/5E**, and random-fragment wide controls; **`--strict`** for non-zero on failure; **`--skip-iedb-pipeline`** skips merged-TSV Fig 5ABC / merged 5DE / sensitivity / combo grid |
| **Clean matrix: Fig 5–6 + all SB modes + instances/unique + sensitivity** | `python scripts/regenerate_manuscript_netmhc_figures.py --clean --purge-repo-figures-netmhc --purge-data-netmhc-figures` — see script `--help` (`--dry-run`, skips, optional `--with-combination-grid`) |

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
| **2** | `plot_tr_de_peptide_fractions_by_transition.py` | `python plot_tr_de_peptide_fractions_by_transition.py --peptide-gene-set tcga_matrix` (or `all_smprot_filtered`) |
| **3A** | `plot_aa_frequency_tcga_vs_proteome.py` | `python plot_aa_frequency_tcga_vs_proteome.py --peptide-set tcga_matrix` |
| **3B** | `plot_dipeptide_volcano_lnc_vs_proteome.py` | `python plot_dipeptide_volcano_lnc_vs_proteome.py --peptide-set tcga_matrix` |
| **3C–3D** | `plot_figure3cd_dipeptide_log2fc_heatmaps.py` | `python plot_figure3cd_dipeptide_log2fc_heatmaps.py` |
| **4A** | `plot_figure4a_tis_vs_ribo_tr_mps.py` | `python plot_figure4a_tis_vs_ribo_tr_mps.py` |

---

## Figure 5 — cohort NetMHCpan (wide XLS; legacy / partial)

Default PNG/CSV under `data/netmhc/figures/` (some scripts also mirror copies into repo-root `figures/`; see each script). **Cohort epitope sharing (5B / 5C)** for the manuscript uses the **merged** pipeline in the next section, not these wide-XLS sharing plots.

| Panel | Script | Notes |
|-------|--------|--------|
| **5A** | `plot_netmhc_epitopes_vs_hla_frequency.py` | Defaults at repo root; `--help` for SB criterion and `--y-metric unique` |
| **5B / 5C** | `scripts/plot_fig5abc_netmhc_sb_triple.py` | **Manuscript default:** merged ``*_with_iedb.tsv``, ``sb_mode=full`` (see `netmhc_sb_core.py`). Two **5C** stems: proportional-whole + ``coding_control`` merged TSVs (second run: ``--panels c``). Legacy wide: ``plot_figure5b_epitope_sharing_across_alleles.py`` (IC50-from-BA only) |
| **5D / 5E** | `plot_figure5de_epitopes_per_allele.py` | **Legacy / supplement:** IC50-from-BA on wide XLS. Manuscript default: ``scripts/plot_fig5de_merged_iedb_sb_per_allele.py`` (two ``--output-stem`` values for proportional-whole vs fragment merged coding) |

---

## Figure 5 — merged `*_with_iedb.tsv` (IEDB + NetMHC SB)

Requires merged long TSVs from `scripts/merge_netmhcpan_xls_with_iedb.py` (see `data/netmhc/README_netmhc.md`).

| Product | Script | Default output dir |
|---------|--------|--------------------|
| **5A–5C** merged | `scripts/plot_fig5abc_netmhc_sb_triple.py` | `--out-dir data/netmhc/figures` (stem `fig5abc_sb_immuno_proc_el_ic50`); second coding cohort: `--panels c --coding-tsv data/netmhc/netmhcpan_coding_control_with_iedb.tsv --output-stem fig5abc_sb_immuno_proc_el_ic50_coding_control` |
| **5D–5E** merged | `scripts/plot_fig5de_merged_iedb_sb_per_allele.py` | Default stem ``fig5de_merged_iedb_sb_proportional_whole``; second cohort: ``--coding-tsv data/netmhc/netmhcpan_coding_control_with_iedb.tsv --output-stem fig5de_merged_iedb_sb_random_fragments``; ``--out-dir`` default `data/netmhc/figures` |
| **Sensitivity / robustness** | `scripts/netmhc_sb_sensitivity_robustness.py` | Under `data/netmhc/figures/` (see `--help`) |
| **Combination grid** | `scripts/plot_fig5_netmhc_sb_combination_grid.py` | `data/netmhc/figures/fig5_netmhc_sb_combinations/` |

---

## Figure 6 — TTN-AS1 (smPEP 108065)

**Canonical manuscript entry point (one script per figure):** `scripts/manuscript_figure6_ttn_as1.py` — forwards argv to `plot_figure6_ttn_as1_allele_coverage.py`.

| Task | Command |
|------|---------|
| Default combined PNG (instances coverage) | `python scripts/manuscript_figure6_ttn_as1.py` |
| Same, explicit path | `python plot_figure6_ttn_as1_allele_coverage.py -o figures/fig6_ttn_as1_allele_coverage.png` |
| Split panels (6A–6E) + unique companion | `python plot_figure6_ttn_as1_allele_coverage.py --split-panels --also-write-unique -o figures/fig6_ttn_as1_split.png` |
| IEDB-gated SB | `python plot_figure6_ttn_as1_allele_coverage.py --gating iedb_sb --iedb-csv <path> --iedb-parent-input-seq-id <id>` (see `--help`; on PowerShell avoid `$PID` as a variable name) |
| **TTN synthetic IEDB companion CSV** | `python scripts/build_ttn_iedb_companion_csv.py` |
| **Fig 6 sensitivity** (NetMHC SB sweeps only) | `python plot_figure6_ttn_as1_sb_sensitivity.py` → `data/netmhc/figures/fig6_ttn_as1_sensitivity/` |
| **Fig 6 from IEDB Tools API (NetMHCpan 4.1 BA+EL)** | `python scripts/fetch_ttn_mhci_iedb_api_netmhc41.py` (see `data/netmhc/predictions/ttn_as1_smpep108065_iedb_api_netmhc41/README.md`; add `--insecure` on Windows if TLS fails) then plot with `--netmhc-xls …/netmhcpan_ttn_as1_iedb_api_netmhc41.xls` |

## Tests

This repo may not ship a unified test suite for plotting. If `pytest` is configured, from the repo root run:

`python -m pytest`

If that reports “no tests collected”, rely on the orchestrators above and spot-check outputs under `figures/` and `data/netmhc/figures/`.
