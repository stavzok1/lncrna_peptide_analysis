# Analysis parameters

Central notes for reproducible choices. Add new sections as workflows grow.

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
