# Pipeline outputs reference

This document lists **useful** CSV/TSV/TXT outputs from the UNDEFINED lncRNA / SmProt / limma workflow: what each file is, where it lives, and how it is produced. Large raw inputs (e.g. `SmProt2.txt`, full expression matrices) are summarized but not exhaustively column-documented.

---

## 1. Gene lists and clinical / expression inputs (`data/`)

| File | Description |
|------|-------------|
| **`canonical_significant_lncRNAs.txt`** | One gene symbol per line: **limma DE gene union ∩ z-score union** (~1608). Written by `tr_limma_de.R` when `limma_z_intersection_genes.txt` exists; also synced by `build_significant_lncs_smprot.py`. **Primary “Tr-lncRNA” gene list.** |
| **`significant_lncs.csv`** | Canonical genes with **`gene_id`** from `lncrna_genes_small.csv` mapping. Built by `build_significant_lncs_smprot.py`. |
| **`primary_exp_stage_lncRNA.csv`** | TCGA samples × meta + **lncRNA gene columns** (log2 expected count + 1) for **AJCC stage** analyses. |
| **`primary_exp_metastasis_lncRNA.csv`** | Same for **M_stage** (M1_s excluded in R scripts). |
| **`lncrna_genes_small.csv`** | GENCODE-style **`gene_id`**, **`gene_name`**, coordinates — used to map symbols ↔ ENSG. |
| **`significant_lnc_peptides_full.tsv`** | Rows from **FASTA-synced** `smprot_filtered.tsv` whose **GeneID / GeneSymbol** match **`significant_lncs.csv`**. Written by **`smprot_gene_match.write_significant_lnc_peptides_full`** when you run **`export_tcga_filtered_peptides_fasta.py --peptides-tsv data/smprot_filtered.tsv ...`** (not by `build_significant_lncs_smprot.py` alone). |
| **`significant_lnc_peptides.tsv`** | **Analyzed / exportable** subset: rows from **`significant_lnc_peptides_full.tsv`** whose `smPEP_ID` appears in **`significant_lnc_peptides.faa`**. The all-filtered export chain re-exports the FASTA and runs **`sync_significant_lnc_peptides_with_fasta.py`** automatically unless **`--no-chain-significant-fasta`** is passed. |
| **`significant_lnc_peptides.faa`** | Translated sequences for the significant peptide table (from **`export_tcga_filtered_peptides_fasta.py --peptides-tsv ..._full.tsv`**; optional **`--min-aa-length 9`** for short ORFs). |
| **`significant_lnc_peptide_gene_map_127.csv`** | **One row per** SmProt `GeneSymbol` present in the **analyzed** `significant_lnc_peptides.tsv` (historically **127** symbols from the full SmProt table; count drops if peptides are removed after FASTA sync): SmProt ENSG(s), transcript IDs, peptide counts/IDs, **canonical** `gene_name` / `gene_id` (via ENSG base match), and **TCGA matrix** symbol / ENSG (usually the canonical symbol as column name). Built by **`export_significant_peptide_gene_table.py`**. Resolves the **127 vs 70** issue: 70 counts genes whose **SmProt symbol string** is a matrix column; all **127** map to a **TCGA column symbol** when resolved through canonical ENSG. |
| **`smprot_filtered.tsv`** | **FASTA-synced** curated list: after **`build_significant_lncs_smprot.py`** (p &lt; 0.05, ≥10 aa, deduped `smPEP_ID` with sorted tie-break), **`export_tcga_filtered_peptides_fasta.py --peptides-tsv data/smprot_filtered.tsv ...`** **rewrites** this file to **only** rows that produced a successful in-silico translation. Downstream TCGA and significant tables assume this definition. |
| **`smprot_filtered_tcga_expr_genes.tsv`** | Subset of `smprot_filtered.tsv` where the peptide maps to a TCGA matrix column by **`GeneSymbol`** (exact header match) **or** by **`GeneID`** → `lncrna_genes_small.csv` **`gene_name`** that is a matrix column (**Ensembl-aware rescue**). Same **`smPEP_ID`** dedup in **`filter_peptides_tcga_expr_genes.py`** (sort by TIS, Ribo, RiboID then `keep="first"`). Column **`TCGA_match_via`** records how each row matched. Diff summaries: **`data/reports/tcga_expr_genes_filter_summary.txt`** (and removed/added TSVs when a prior output existed). |
| **`report_smprot_pipeline_stages.py`** → **`data/reports/smprot_pipeline_stages.{md,json}`** | Per-stage counts: all-filtered, TCGA-mapped, significant (gene-list) — rows kept, matrix exclusions, FASTA successes/failures. Narrative: **`docs/smprot_peptide_pipeline_stages.md`**. |
| **`smprot_gene_match.py`** | `write_significant_lnc_peptides_full()` — builds **`significant_lnc_peptides_full.tsv`** from FASTA-synced **`smprot_filtered.tsv`** + **`significant_lncs.csv`**. Invoked from **`export_tcga_filtered_peptides_fasta.py`** after restricting `smprot_filtered.tsv`. |
| **`smprot_filtered_no_tcga_summary.csv`** | Aggregate: unique peptide count and gene count from **`smprot_filtered.tsv`** (no TCGA filter). From **`export_filtered_peptide_stats.py`**. |
| **`smprot_filtered_no_tcga_peptides_by_gene.csv`** | Per-`GeneSymbol` peptide counts from **`smprot_filtered.tsv`**. |
| **`smprot_tcga_filtered_peptides.faa`** | Amino acid sequences exported from GENCODE transcript slice + translate for rows in **`smprot_filtered_tcga_expr_genes.tsv`**. See `export_tcga_filtered_peptides_fasta.py`. |
| **`compare_old_new/*.txt`**, **`comparison_summary.json`** | Old vs new gene/peptide set overlaps. Produced by **`compare_sig_old_new.py`**. |
| **`compare_old_new/significant_peptide_genes_tcga_old_new.csv`** | Rows comparing **new** vs **old** significant peptide sets with/without TCGA symbol filter; gene intersection lists. From **`export_filtered_peptide_stats.py`**. |

---

## 2. Z-score transition outputs (`tr_lncrna_output/`)

| File | Description |
|------|-------------|
| **`tr_lncrna_de_analysis.py`** → **`tr_lncrnas_stage_detail.csv`**, **`tr_lncrnas_metastasis_detail.csv`** | Per **cancer × transition × gene**: `log2FC`, **z** (|z|≥3 retained), sample counts. Feeds Tr-gene union and peptide-fraction plots (limma ∩ z per stratum). |
| **`tr_genes_union.txt`**, **`tr_genes_stage_unique.txt`**, **`tr_genes_metastasis_unique.txt`** | Z-union and stage/meta–specific Tr gene lists. |
| **`tr_lncrna_summary.json`** | Counts and cancer lists from z pipeline. |

---

## 3. Limma DE (`tr_lncrna_output/limma/`)

| File | Description |
|------|-------------|
| **`limma_stage_all_tests.csv`**, **`limma_metastasis_all_tests.csv`** | All moderated tests: `gene`, `logFC`, `adj.P.Val`, `cancer_type`, `transition`, **`n_early`**, **`n_late`**. |
| **`limma_stage_FDR0.05.csv`**, **`limma_metastasis_FDR0.05.csv`** | FDR **&lt; 0.05** hits only. |
| **`limma_z_intersection_genes.txt`** | Sorted intersection of z-union genes with limma DE union (feeds canonical list). |
| **`limma_summary.json`** | Run metadata and counts. |
| **`limma_z_genes_tcga_filtered_peptide_summary.csv`** | For **each canonical gene**: `n_unique_filtered_peptides`, `smPEP_IDs` (from **`smprot_filtered_tcga_expr_genes.tsv`**). From **`summarize_intersection_genes_peptides.py`**. |
| **`limma_z_genes_tcga_filtered_peptide_summary_min1_by_npeps.csv`** | Same, **only genes with ≥1** peptide, sorted by count. |

---

## 4. Figures and ancillary (`figures/`, `tr_lncrna_output/figures/`, scripts at repo root)

**Figure-level documentation:** **`docs/figure_catalog.md`**. **Catalog Fig. 2 & 3** primary outputs live under **`figures/tcga_matrix/`** and **`figures/all_smprot_filtered/`** (see **`generate_catalog_figures.py`**). Legacy/auxiliary plots may still use **`tr_lncrna_output/figures/`**.

| Output | Script |
|--------|--------|
| **`figures/<mode>/peptide_fraction/*.png`**, **`figures/<mode>/tr_de_peptide_fraction_by_cancer.csv`** | **Fig. 2** (`<mode>` = `tcga_matrix` or `all_smprot_filtered`). **`plot_tr_de_peptide_fractions_by_transition.py`**. |
| **`figures/<mode>/aa_frequency_*_vs_known_proteins.*`** | **Fig. 3A**. **`plot_aa_frequency_tcga_vs_proteome.py`**. |
| **`figures/<mode>/dipeptide_volcano_lnc_vs_proteome*.*`** | **Fig. 3B**. **`plot_dipeptide_volcano_lnc_vs_proteome.py`**. |
| **`figures/fig3c_dipeptide_log2fc_tcga_matrix_vs_proteome.png`**, **`figures/fig3c_dipeptide_log2fc_all_smprot_filtered_vs_proteome.png`** (if FASTA exists), **`figures/fig3d_dipeptide_log2fc_tr_lncrna_tcga_vs_proteome.png`** | **Fig. 3C–3D** log2FC 20×20 heatmaps. **`plot_figure3cd_dipeptide_log2fc_heatmaps.py`**. |
| **`figures/fig4a_tr_lncrna_mp_tis_vs_riboseq_pvalues.png`** | **Fig. 4A** TIS vs Ribo-seq p-values (Tr MPs). **`plot_figure4a_tis_vs_ribo_tr_mps.py`**. |
| **`generate_catalog_figures.py`** | Runs the three scripts per mode + **3C–3D** + **4A** (skips `all_smprot_filtered` mode if FASTA missing unless **`--strict`**). |
| **`fig2_dipeptide_mp_composition.png`**, **`fig2_dipeptide_summary.txt`** | Dipeptide heatmaps vs **`data/known_proteins.fasta`** (TCGA-filtered MPs, canonical Tr subset, and **significant exportable FASTA** when present). **`plot_dipeptide_mp_figure2.py`** → **`tr_lncrna_output/figures/`**. |
| **`starting_aa_filtered_vs_significant.png`**, **`starting_aa_distribution_counts.csv`**, **`starting_aa_distribution_report.txt`** | N-terminal AA distribution (filtered all, canonical subset, **SmProt significant exportable** when **`significant_lnc_peptides.faa`** exists). **`plot_starting_aa_distribution_mps.py`**. |
| **`data/cdhit_clustering/<GENE>/`** | Per-gene CD-HIT clustering. **`cluster_gene_peptides_cdhit.py`**. |

---

## 5. Parameter notes (`docs/`)

| File | Description |
|------|-------------|
| **`docs/analysis_params.md`** | CD-HIT `-n`/`-l`, peptide fraction plot definitions, dipeptide Fig. 2 formulas, transcript slice convention for FASTA export. |

---

## 6. NetMHCpan inputs (`data/netmhc/`)

| File | Description |
|------|-------------|
| **`hla_european27_class1.txt`** | **27** class I alleles (one per line, **LF**). Shells expand this to a **comma-separated** `-a` argument (do **not** pass the file path as a single allele). Replace with your paper’s list if needed. |
| **`README_netmhc.md`** | WSL install, FASTA prep, **§5** full `netMHCpan` examples (including proportional-whole). |
| **Shell parameter summary** | **`docs/iedb_tools_api.md`** → *Appendix: cohort wide XLS — local NetMHCpan-4.2* (same defaults as `run_netmhcpan42_example.sh` / `run_netmhcpan_ttn_as1_108065.sh`). |
| **`sig_peptides_missing_from_fasta.tsv`** | Optional (**`data/netmhc/`**): written by **`prepare_netmhc_tr_vs_coding_epitopes.py`** if the analyzed TSV still lists `smPEP_ID`s absent from the FASTA (should be empty after **`sync_significant_lnc_peptides_with_fasta.py`**). |
| **`sig_parent_micropeptides.csv`**, **`coding_control_parent_micropeptides.csv`**, **`ninemers_sig_lnc.fasta`**, **`ninemers_coding_control.fasta`**, **`ninemers_summary.csv`** | Produced by **`prepare_netmhc_tr_vs_coding_epitopes.py`**: all significant lncRNA peptides vs length-matched substrings from **`data/known_proteins.fasta`** (override with `--coding-fa`); sliding **9-mer** FASTA for local NetMHCpan-4.1. |

---

## 7. How to refresh key tables

```text
Rscript tr_limma_de.R <project_root>
python tr_lncrna_de_analysis.py
python build_significant_lncs_smprot.py
python export_tcga_filtered_peptides_fasta.py --peptides-tsv data/smprot_filtered.tsv --out-aa data/smprot_all_filtered_peptides.faa
python filter_peptides_tcga_expr_genes.py
python export_tcga_filtered_peptides_fasta.py
python summarize_intersection_genes_peptides.py
# significant_lnc_peptides.faa + .tsv: normally refreshed by the smprot_filtered export line above (omit the next two unless you used --no-chain-significant-fasta):
python export_tcga_filtered_peptides_fasta.py --peptides-tsv data/significant_lnc_peptides_full.tsv --transcripts-fa <gencode.transcripts.fa.gz> [--ensembl-fallback] [--min-aa-length 9]
python sync_significant_lnc_peptides_with_fasta.py
python export_significant_peptide_gene_table.py
python export_filtered_peptide_stats.py
python compare_sig_old_new.py
python prepare_netmhc_tr_vs_coding_epitopes.py
python plot_dipeptide_mp_figure2.py
python plot_starting_aa_distribution_mps.py
```

Order can vary slightly (e.g. limma before or after z script for intersection file). **`significant_lnc_peptides_full.tsv`** is written when you export from **`smprot_filtered.tsv`** (after **`significant_lncs.csv`** exists). **`significant_lnc_peptides.tsv`** is updated by the same export chain (FASTA + **`sync_significant_lnc_peptides_with_fasta.py`**) unless you skip it with **`--no-chain-significant-fasta`**.

---

## 8. 127 vs 70 genes (recap)

- **127** = distinct **SmProt `GeneSymbol`** strings in **`significant_lnc_peptides.tsv`**.
- **70** = genes whose **SmProt symbol alone** appears as a **column name** in the TCGA matrices (strict string match).
- Many peptides use **alternate symbols** for the same ENSG as the canonical list; **`significant_lnc_peptide_gene_map_127.csv`** adds **`canonical_gene_name`** / **`tcga_matrix_gene_symbol`** (usually identical) via **versionless ENSG** match to **`significant_lncs.csv`**, so all **127** rows can show the **TCGA** symbol used in expression data when the locus is in the panel.
