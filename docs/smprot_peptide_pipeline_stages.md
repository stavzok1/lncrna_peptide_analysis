# SmProt peptide processing (stages)

This page describes **three parallel peptide products** derived from SmProt lncRNA rows, and where counts are recorded.

## Where the numbers live

After any rebuild of `smprot_filtered.tsv`, TCGA filtering, or FASTA export/sync, refresh the quantitative report:

```bash
python report_smprot_pipeline_stages.py
```

That writes:

- `data/reports/smprot_pipeline_stages.md` — **Snapshot** table at the top with **live numbers** (rows, FASTA records, exclusions, `TCGA_match_via` counts), then per-stage detail.
- `data/reports/smprot_pipeline_stages.json` — same content for scripts

## The three branches (conceptual)

1. **All filtered** — `build_significant_lncs_smprot.py` then **`export_tcga_filtered_peptides_fasta.py --peptides-tsv data/smprot_filtered.tsv --out-aa data/smprot_all_filtered_peptides.faa`**  
   Build writes a **pre-export** `smprot_filtered.tsv` (p &lt; 0.05, length, deduped `smPEP_ID`). The export step **rewrites** that TSV to **only peptides that received a FASTA record** (translation failures are dropped from the table, not only from the FASTA). The same export run refreshes `significant_lnc_peptides_full.tsv` and, by default, re-exports `significant_lnc_peptides.faa` and runs `sync_significant_lnc_peptides_with_fasta.py`.

2. **TCGA matrix–mapped** — `data/smprot_filtered_tcga_expr_genes.tsv` → `data/smprot_tcga_filtered_peptides.faa`  
   Same starting pool, restricted to genes that appear as **columns** in the TCGA lncRNA expression matrices (`filter_peptides_tcga_expr_genes.py`). Matching uses **gene symbol** and/or **Ensembl `GeneID` → `gene_name`** via `data/lncrna_genes_small.csv`.

3. **Significant (canonical *gene* list)** — `data/significant_lnc_peptides_full.tsv` → `data/significant_lnc_peptides.faa` → `data/significant_lnc_peptides.tsv`  
   Peptides whose host gene matches **`data/significant_lncs.csv`** (built from the limma DE ∩ z-score gene union in `build_significant_lncs_smprot.py`). **“Canonical” here means that gene-list match**, not TCGA column presence. The analyzed TSV is the subset that **also** has a successful FASTA translation (`sync_significant_lnc_peptides_with_fasta.py`).

## Other diffs

- When re-running `filter_peptides_tcga_expr_genes.py`, see `data/reports/tcga_expr_genes_filter_summary.txt` and companion TSVs for **added/removed** `smPEP_ID` relative to the previous `smprot_filtered_tcga_expr_genes.tsv` on disk.
