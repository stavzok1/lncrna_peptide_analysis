# Data layout

Small **SmProt / TCGA** tables, **NetMHC helper** FASTAs/metadata, **`netmhcpan_ttn_as1_108065.xls`** (TTN-AS1 wide run), and **saved figure CSV sidecars** under `data/netmhc/figures/` ship in this repository. **Large** cohort wide `*.xls`, merged `*_with_iedb.tsv`, and **IEDB peptide_table** exports are gitignored: copy them from your full analysis tree or from an archive (e.g. Zenodo) into `data/netmhc/` as listed below.

## Required for `python generate_netmhc_figure_bundle.py`

Place these under `data/netmhc/` (paths match `generate_netmhc_figure_bundle.py` and `docs/netmhc_figure_commands.md`):

| File | Role |
|------|------|
| `netmhcpan_sig_lnc_with_iedb.tsv` | Merged cohort — significant lncRNA MPs |
| `netmhcpan_coding_proportional_whole_with_iedb.tsv` | Merged coding — proportional whole |
| `netmhcpan_coding_control_with_iedb.tsv` | Merged coding — length-matched control |
| `netmhcpan_ttn_as1_108065.xls` | Wide NetMHCpan output for TTN-AS1 smPEP 108065 (also bundled here if small) |

Optional **wide XLS** cohort panels: original `netmhcpan_sig_lnc.xls` and coding `*.xls` if you use **`python generate_netmhc_supplement.py --include-wide-xls-fig5`**.

**Allele frequencies** for Figure 5A: default `data/netmhc/figures/fig5a_epitopes_vs_allele_frequency_ic50_sb.csv` (or `epitopes_vs_allele_frequency_ic50_sb.csv`). Regenerate from wide 5A if missing.

## Required for `python generate_catalog_figures.py`

Repo includes `data/smprot_filtered_tcga_expr_genes.tsv`, `data/smprot_tcga_filtered_peptides.faa`, `data/known_proteins.fasta`, and gene lists. Optional: `data/smprot_all_filtered_peptides.faa` for `all_smprot_filtered` mode (build with `python pipeline/export_tcga_filtered_peptides_fasta.py`, see orchestrator docstring).

## NetMHC prep and merge

See `data/netmhc/README_netmhc.md` and `scripts/merge_netmhcpan_xls_with_iedb.py` to rebuild `*_with_iedb.tsv` from wide `*.xls` plus IEDB CSVs.
