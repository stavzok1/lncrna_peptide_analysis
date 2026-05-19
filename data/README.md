# Data layout

Small **SmProt / TCGA** tables, **NetMHC helper** FASTAs/metadata, **`netmhcpan_ttn_as1_108065.xls`** (TTN-AS1 wide run), and **canonical merged** CSV sidecars under `data/netmhc/figures/` (`fig5_merged_*`, `fig5de_merged_whole_*`, supplement stems) ship in this repository. **All PNGs** under `data/netmhc/figures/` stay on GitHub; **legacy wide-XLS / Cartesian-grid / alternate-stem CSVs** are gitignored (see root `.gitignore`). **Large** cohort wide `*.xls`, merged `*_with_iedb.tsv`, and **IEDB peptide_table** exports are gitignored: copy them from your full analysis tree or from an archive (e.g. Zenodo) into `data/netmhc/` as listed below. **Which tables go to Zenodo vs GitHub** is summarized in **`docs/ZENODO.md`** (*What goes where*).

A one-time copy log from the parent analysis tree (if you used it) may exist locally as **`data/MANIFEST_COPIED.txt`** (not tracked on GitHub).

---

## Expression tables (TCGA lncRNA matrices)

These drive **Fig 1B**, **Fig 2**, and **`generate_tr_lncrna_identification.py`** (Python z-screen + R limma). Values are **log2(expected count + 1)** per gene column (same scale as `pipeline/tr_lncrna_de_analysis.py`: prevalence filter **log2 ≥ 1** in ≥ **40%** of union samples per stratum).

| File | Meta columns (first columns) | Notes |
|------|------------------------------|--------|
| **`data/primary_exp_stage_lncRNA.csv`** | `sample_id`, `cancer_type`, `ajcc_t`, `ajcc_m`, **`stage`** | AJCC stage labels; used with `stage` for early/late transitions |
| **`data/primary_exp_metastasis_lncRNA.csv`** | `sample_id`, `cancer_type`, **`stage`**, **`M_stage`** | **`M1_s` excluded** in downstream code; metastasis transitions use `M_stage` |

Remaining columns are **lncRNA gene symbols** (GENCODE-style IDs). Cohort sizes per cancer follow **`MIN_SAMPLES_CANCER` > 100** in the DE scripts.

**Optional wider matrices** (not required by default figure orchestrators): e.g. `data/primary_exp_stage_final.csv`, `data/primary_exp_metastasis_final.csv` — full gene sets or alternate filters.

**How these CSVs were built:** notebook **`process_scratch.ipynb`**. In the full **UNDEFINED** checkout it usually lives **next to** `paper-github/` as **`../process_scratch.ipynb`**. This bundle also keeps **`notebooks/process_scratch.ipynb`** — copy from the parent when you update the scratch workflow. It covers TCGA matrix import, **sample** and **lncRNA** filtering, and **cancer type / stage / M_stage** annotation to match the column contracts above (paths to raw downloads are local to your machine).

---

## Required for `python generate_catalog_figures.py`

Orchestrator: Fig **1B** (t-SNE) + **2–4A**. Run from repo root; paths are under `data/` and `tr_lncrna_output/` unless overridden. **Limma / z tables for Fig 2:** run **`python generate_tr_lncrna_identification.py`** first if `tr_lncrna_output/limma/` is empty (see Tr-lncRNA identification section below).

| Path | Role |
|------|------|
| **`data/primary_exp_stage_lncRNA.csv`** | Fig **1B** t-SNE + same sample set as DE pipeline |
| **`data/primary_exp_metastasis_lncRNA.csv`** | Metastasis panels in Fig **2** (peptide fractions) |
| **`data/smprot_filtered_tcga_expr_genes.tsv`** | Peptide gene table (**tcga_matrix** mode; Fig 2 / 3A–3B) |
| **`data/significant_lnc_peptides.tsv`** | **Fig. 4A** default (~501 exportable MPs; NetMHC cohort) |
| **`data/smprot_tcga_filtered_peptides.faa`** | FASTA for **3A–3B**, **3C** (TCGA-matrix branch), **3D** |
| **`data/known_proteins.fasta`** | Proteome reference for composition / volcano / heatmaps |
| **`tr_lncrna_output/limma/*.csv`** | Limma DE tables (e.g. `limma_stage_FDR0.05.csv`, `limma_metastasis_FDR0.05.csv`) |
| **`tr_lncrna_output/tr_lncrnas_stage_detail.csv`**, **`tr_lncrna_output/tr_lncrnas_metastasis_detail.csv`** | z-detail tables merged with limma for Fig **2** |
| **`data/canonical_significant_lncRNAs.txt`** (or **`tr_lncrna_output/limma/limma_z_intersection_genes.txt`**) | Canonical Tr gene list for dashed-line / Fig **3D** logic |

**Optional — `all_smprot_filtered` mode (`--only` or full run):** `data/smprot_all_filtered_peptides.faa` (see orchestrator docstring for build command).

**Optional / not used by default figure scripts:** `data/primary_exp_stage_final.csv`, `data/primary_exp_metastasis_final.csv`, etc. (full expression matrices for other analyses).

---

## Tr-lncRNA identification (z-scores + limma)

Run from repo root:

```bash
python generate_tr_lncrna_identification.py
```

**Inputs:** same **`data/primary_exp_stage_lncRNA.csv`** and **`data/primary_exp_metastasis_lncRNA.csv`** as Fig 1B / Fig 2.

**Step 1 — Python z-screen** (`pipeline/tr_lncrna_de_analysis.py`): per cancer and transition, expression prevalence filter (log2(expr+1) ≥ 1 in ≥ 40% of union samples), mean log2FC early vs late, **z-score across genes** in that stratum, keep **|z| ≥ 3**. Writes:

| Output | Content |
|--------|---------|
| `tr_lncrna_output/tr_lncrnas_stage_detail.csv` | Stage hits (gene, cancer, transition, log2FC, z, n) |
| `tr_lncrna_output/tr_lncrnas_metastasis_detail.csv` | M_stage hits (M1_s excluded) |
| `tr_lncrna_output/tr_genes_union.txt` | Union of z-hit gene symbols |
| `tr_lncrna_output/tr_lncrna_summary.json` | Counts / parameters |
| `tr_lncrna_output/*.png` | Quick diagnostic bar charts |

**Step 2 — R limma** (`tr_limma_de.R`): same matrices, same prevalence filter, **moderated t-tests** (BH **FDR < 0.05**). Always writes under **`tr_lncrna_output/limma/`** (`limma_stage_FDR0.05.csv`, `limma_metastasis_FDR0.05.csv`, `limma_de_genes_union_FDR0.05.txt`, all-tests tables, `limma_summary.json`, histograms).

**If `tr_lncrna_output/tr_genes_union.txt` exists** (from step 1): R intersects that z-gene union with the limma DE union (FDR &lt; 0.05), writes **`tr_lncrna_output/limma/limma_z_intersection_genes.txt`**, and if that intersection is non-empty overwrites **`data/canonical_significant_lncRNAs.txt`**.

**If `tr_genes_union.txt` is missing** (e.g. you ran limma before the z-step, or deleted the file): limma outputs above are still produced; **no** `limma_z_intersection_genes.txt` is written; **`canonical_significant_lncRNAs.txt` is not updated** by R (downstream scripts may still use an older canonical list if present). Run **`python pipeline/tr_lncrna_de_analysis.py`** first, then R again.

**R packages:** `limma` (Bioconductor), `jsonlite` (CRAN). From repo root: **`Rscript install_r_dependencies.R`** (see root **`requirements.txt`** comments).

**Next (SmProt / MPs, not in the orchestrator):** `python pipeline/build_significant_lncs_smprot.py` — needs **`data/SmProt2.txt`**, `data/lncrna_genes_small.csv`, etc. **`SmProt2.txt`** and **`primary_exp_*_final.csv`** are **not** tracked on GitHub (see `.gitignore`); download them from the **Zenodo dataset** (see **`docs/ZENODO.md`**) into your local `data/` before running that pipeline.

---

## Required for `python generate_netmhc_figure_bundle.py`

Place these under `data/netmhc/` (paths match `generate_netmhc_figure_bundle.py` and `docs/netmhc_figure_commands.md`):

| File | Role |
|------|------|
| `netmhcpan_sig_lnc_with_iedb.tsv` | Merged cohort — significant lncRNA MPs |
| `netmhcpan_coding_proportional_whole_with_iedb.tsv` | Merged coding — proportional whole |
| `netmhcpan_coding_control_with_iedb.tsv` | Merged coding — length-matched control |
| `netmhcpan_ttn_as1_108065.xls` | Wide NetMHCpan output for TTN-AS1 smPEP 108065 (also bundled here if small) |
| `hla_european27_allele_frequencies.csv` | Reference **allele × population frequency** for merged Fig **5A** (and optional wide 5A in supplement) |

Optional **wide XLS** cohort panels: `netmhcpan_sig_lnc.xls` and coding `*.xls` if you use **`python generate_netmhc_supplement.py --include-wide-xls-fig5`**.

**Supplement sensitivity / grids** (`python generate_netmhc_fig5_fig6_supplement.py`, or `regenerate_all_figures.py`): needs the **same merged cohort TSVs** as the canonical NetMHC Fig 5 scripts, plus TTN wide XLS + IEDB companion CSV for Fig 6 supplement steps. **Rebuild merged TSVs** from wide XLS + IEDB CSVs: `python rebuild_netmhc_merged_tsvs.py` (see `--help`; IEDB CSVs often live only on your machine under a sibling `data/netmhc/`).

**Allele frequencies (legacy sidecar):** older workflows used `data/netmhc/figures/fig5a_epitopes_vs_allele_frequency_ic50_sb.csv` from wide 5A; the **default merged pipeline** now uses the bundled **`hla_european27_allele_frequencies.csv`**.

## NetMHC prep and merge

See `data/netmhc/README_netmhc.md` and `scripts/merge_netmhcpan_xls_with_iedb.py` to rebuild `*_with_iedb.tsv` from wide `*.xls` plus IEDB CSVs.

---

## Git (Windows): CRLF warnings on `*.tsv`

If you see **`CRLF will be replaced by LF the next time Git touches it`** when staging NetMHC TSVs, that is **Git + `.gitattributes`** (`*.tsv` → `eol=lf`) normalizing line endings. It is harmless, but noisy.

1. **Do not commit** large merged TSVs or **`merge_rebuild_verify/`** outputs — they belong under **`data/netmhc/**/*.tsv`** / **`merge_rebuild_verify/`** in `.gitignore`. If you already staged them, remove from the index only:

   ```bash
   git rm -r --cached data/netmhc/merge_rebuild_verify
   ```

   (Add any other wrongly tracked paths the same way.)

2. **Clear the warning** for files you *do* keep: save those TSVs as **LF** (e.g. open in an editor and set line ending to LF, or use `dos2unix` in Git Bash), then `git add` again.

3. Optional — stop Git from scanning huge ignored trees: after `git rm --cached`, run `git status` and confirm those paths show as **untracked** / ignored, not staged.
