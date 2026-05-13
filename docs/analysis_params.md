# Analysis parameters

Central notes for reproducible choices. This file is **not** a complete catalog of every matplotlib knob for every figure: it grew **ad hoc** with deep notes where logic is easy to misread from code alone (e.g. **Fig 2** limma ∩ z joins, **supplement dipeptide** Cohen’s h, **CD-HIT** word-size table).

## Where each figure’s parameters are documented

There is **no single markdown** that lists every argparse flag for every panel. Use this index plus each script’s **`python manuscript/<script>.py --help`**.

| Figure | Narrative (what it shows, paths, key thresholds) | Extra “parameter” notes | Full CLI |
|--------|--------------------------------------------------|---------------------------|----------|
| **1B** | `docs/figure_catalog.md` § Figure 1B | openTSNE / perplexity etc. only in script | `manuscript/plot_figure1b_tsne_stage_lncrna.py --help` |
| **2** (bars) | `docs/figure_catalog.md` § Figure 2 | **This file** — *Peptide fraction bar charts* (denominator/numerator, dashed lines, `MIN_SAMPLES_CANCER`) | `manuscript/plot_tr_de_peptide_fractions_by_transition.py --help` |
| **3A–3B** | `docs/figure_catalog.md` § Figure 3A–3B | Heatmap percentile note in catalog for 3C | `manuscript/plot_aa_frequency_tcga_vs_proteome.py --help`, `plot_dipeptide_volcano_lnc_vs_proteome.py --help` |
| **3C–3D** | `docs/figure_catalog.md` § Figure 3C–3D | **This file** — *Dipeptide Fig. 2* (`plot_dipeptide_mp_figure2.py`) for supplement dipeptide math | `manuscript/plot_figure3cd_dipeptide_log2fc_heatmaps.py --help` |
| **4A** | `docs/figure_catalog.md` § Figure 4A (p-axis range, shading, labels, `--top-extreme-labels`) | — | `manuscript/plot_figure4a_tis_vs_ribo_tr_mps.py --help` |
| **5** (merged SB) | `docs/figure_catalog.md` § Figure 5 | **`scripts/netmhc_sb_core.py`** (`FIG5_IEDB_*`); local wide XLS flags in **`docs/iedb_tools_api.md`** appendix; **`data/netmhc/README_netmhc.md`** §5 | `manuscript/plot_fig5abc_netmhc_sb_triple.py --help`, `plot_fig5de_merged_iedb_sb_per_allele.py --help` |
| **5** (supplement / sensitivity) | `docs/figure_catalog.md` (sensitivity, combo grid) | **`docs/netmhc_figure_commands.md`** | `supplement/netmhc_sb_sensitivity_robustness.py --help`, etc. |
| **6** | `docs/figure_catalog.md` § Figure 6 | **`docs/figure6_ttn_as1_parameters.md`** (API vs local XLS, SB flags) | `manuscript/plot_figure6_ttn_as1_allele_coverage.py --help` |
| **Orchestrators** | — | **`docs/netmhc_figure_commands.md`**, root **`README.md`** | `generate_catalog_figures.py --help`, `generate_netmhc_figure_bundle.py --help`, `generate_netmhc_supplement.py --help` |

**SmProt / export geometry:** **`docs/smprot_peptide_pipeline_stages.md`** and **This file** — *Transcript slice for MP FASTA*. **Artifact paths:** **`docs/pipeline_outputs_reference.md`**.

## TCGA expression CSVs (`data/primary_exp_*_lncRNA.csv`)

- **Provenance / filtering notebook:** `process_scratch.ipynb` — primary copy often at **`../process_scratch.ipynb`** (parent of `paper-github/`); bundle mirror: **`notebooks/process_scratch.ipynb`**. Column layout and log2 scale must match `pipeline/tr_lncrna_de_analysis.py`.
- **File checklist:** `data/README.md` → *Expression tables (TCGA lncRNA matrices)*.

## IEDB Tools API (MHC-I binding fetch — NetMHCpan 4.1 BA+EL)

- **Endpoint, POST fields, fetch CLI:** `docs/iedb_tools_api.md`.
- **Fig 6 after the XLS (SB, immuno/proc gates, `--help` flags):** `docs/figure6_ttn_as1_parameters.md` — these are **not** sent to the `mhci/` binding API; immuno/proc come from **peptide_table CSV** columns when `--gating iedb_sb`.
- **Fetch script:** `supplement/fetch_ttn_mhci_iedb_api_netmhc41.py`.

## NetMHCpan cohort wide XLS (local 4.2 shells)

- **Shell defaults** (`NETMHCpan`, `TMPDIR`, `-inptype 0`, `-l 9`, `-BA 1`, `-pathogen 1`, `-neo 1`, `-context 0`, `-xls 1`, comma-built `-a`, per-run `-f` / `-xlsfile`): **`docs/iedb_tools_api.md`** → *Appendix: cohort wide XLS — local NetMHCpan-4.2*.
- **Long-form** (install, WSL, proportional-whole one-liners): **`data/netmhc/README_netmhc.md`** §5.
- **Scripts:** `data/netmhc/run_netmhcpan42_example.sh` (sig lnc + coding control), `run_netmhcpan_ttn_as1_108065.sh` (TTN 108065).

## CD-HIT (`cluster_gene_peptides_cdhit.py`)

- **Program:** `cd-hit` (protein), invoked via **WSL** on Windows paths mapped under `/mnt/...`.
- **Identity:** `-c 0.6` and `-c 0.9` (global sequence identity; see CD-HIT manual).
- **Word size `-n`:** follows the official CD-HIT protein table ([User’s Guide — Choice of word size for `cd-hit`](https://github.com/weizhongli/cdhit/wiki/3.-User's-Guide)):

  | Identity threshold `-c` | Recommended `-n` |
  |-------------------------|------------------|
  | 0.7 ~ 1.0               | 5                |
  | 0.6 ~ 0.7               | 4                |
  | 0.5 ~ 0.6               | 3                |
  | 0.4 ~ 0.5               | 2                |

  So for this repo: **`-c 0.9` → `-n 5`**, **`-c 0.6` → `-n 4`**.
- **Minimum length `-l`:** **10** (cd-hit's default for `cd-hit`; sequences shorter than `-l` are skipped by CD-HIT). SmProt peptides used here are filtered to **≥10 aa** before FASTA export, so inputs align with this default.
- **Other flags in script:** `-d 0`, `-T 0`, `-M 16000` (see script `run_cdhit_wsl` for the exact command line).

Examples in the same upstream guide use e.g. `cd-hit ... -c 0.9 -n 5 ...` and `cd-hit ... -c 0.6 -n 4 ...`.

## Peptide fraction bar charts (`plot_tr_de_peptide_fractions_by_transition.py`) — **Figure 2** in `docs/figure_catalog.md`

- **Per-cancer bar (denominator):** genes that are **both** (1) limma FDR `< 0.05` for that **cancer type × transition** (`tr_lncrna_output/limma/limma_stage_FDR0.05.csv` or `limma_metastasis_FDR0.05.csv`) and (2) in **`tr_lncrna_output/tr_lncrnas_stage_detail.csv`** or **`tr_lncrnas_metastasis_detail.csv`** for the **same** cancer and transition with **|z| ≥ 3** (`tr_lncrna_de_analysis.py`). Rows are inner-joined on `gene`. **Up / down:** limma `logFC` and z-side `log2FC` are both `> 0` or both `< 0`. **Combined:** limma DE in either direction for that cancer × transition, still intersected with the z table for that stratum (|z| ≥ 3).
- **Per-cancer bar (numerator):** genes in that stratum with ≥1 peptide: unique `GeneSymbol` in the chosen SmProt TSV. **Default (`--peptide-gene-set tcga_matrix`):** `data/smprot_filtered_tcga_expr_genes.tsv` (TCGA-matrix genes; aligns with `smprot_tcga_filtered_peptides.faa`). **Alternate (`--peptide-gene-set all_smprot_filtered`):** `data/smprot_filtered.tsv` (FASTA-synced curated list: only peptides with a successful `smprot_all_filtered_peptides.faa` record).
- **Red dashed line (Overall Tr-lncRNAs):** `100 * |canonical ∩ peptide| / |canonical|` where canonical = `data/canonical_significant_lncRNAs.txt` (global limma gene union **∩** z union, ~1608), using the **same** peptide TSV as the bars.
- **Green dashed line (Overall TCGA lncRNAs):** `100 * |matrix lncRNA columns ∩ peptide| / |matrix lncRNA columns|` (`primary_exp_stage_lncRNA.csv` minus meta columns).
- **Cancers on the x-axis:** `> 100` samples (`MIN_SAMPLES_CANCER`); metastasis counts exclude `M1_s` when listing cohorts.
- **Outputs:** `figures/<peptide_gene_set>/peptide_fraction/*.png` and `figures/<peptide_gene_set>/tr_de_peptide_fraction_by_cancer.csv` (repo root `figures/`). Override parent with `--figures-dir`.

## Dipeptide Fig. 2 (`plot_dipeptide_mp_figure2.py`)

- **Reference:** `data/known_proteins.fasta` (overlapping dipeptide counts on standard 20 aa only).
- **MPs:** `data/smprot_tcga_filtered_peptides.faa` (filtered TCGA-panel MPs from the export script). **Tr-lncRNA-MPs:** same FASTA restricted to `GeneSymbol` in `canonical_significant_lncRNAs.txt`.
- **Panels A–B:** `log2` of frequency ratio (Dirichlet-style `+0.5` pseudocount per cell over 400 pairs before normalizing).
- **Panels C–D:** **Cohen's h** for two proportions per cell: `h = 2*arcsin(sqrt(p_mp)) - 2*arcsin(sqrt(p_ref))` with `p = (c+0.5)/(N+200)` (`N` = total overlapping dipeptides in that corpus). This is the standard effect for comparing proportions; it is not the same as Cohen's d for continuous data and is not a Wald z-statistic (those can exceed 1 arbitrarily for rare categories).
- **Output:** `tr_lncrna_output/figures/fig2_dipeptide_mp_composition.png` and `fig2_dipeptide_summary.txt`.

## Transcript slice for MP FASTA (`export_tcga_filtered_peptides_fasta.py`)

- Coordinates are **0-based, half-open** on the GENCODE transcript sequence as in the FASTA file; see module docstring. Translation is **forward** on that string (no reverse complement in this exporter).
