# Supplementary dataset (large inputs & NetMHC / IEDB artifacts)

This archive accompanies the **code** repository (GitHub, and optionally the GitHub→Zenodo software record). It holds files that are **gitignored** from the public repo because of size or policy (see `docs/ZENODO.md` in the code bundle).

**Suggested Zenodo metadata**

- **Resource type:** Dataset  
- **Title (example):** Supplementary data: TCGA lncRNA expression matrices, SmProt raw export, NetMHCpan wide outputs, IEDB inputs, merged epitope tables  
- **License:** Choose an appropriate **data license** (e.g. CC BY 4.0) in Zenodo; it may differ from the code repo (often MIT).  
- **Related identifiers:** Link the Zenodo **software** DOI (GitHub release) as “is supplemented by this upload” or vice versa.

### One zip vs many separate files on Zenodo

**Yes — you do not have to use a single zip.** One Zenodo **dataset** record can list **many files**; each file gets its own download link, so people can grab only what they need.

**Zenodo limits (standard deposit):** up to **100 files** and **50 GB** total per record ([Manage files](https://help.zenodo.org/docs/deposit/manage-files/)). If your bundle has **more than 100 files** (for example a large `predictions/` tree), use **one or more zips**, apply for a higher quota, or split into a second record.

**Folders:** the web uploader is built around **files**, not preserving a full directory tree as nested Zenodo “folders.” The numbered directories here are for **your** layout before upload. If you upload **loose files**, use **clear, unique names** (or prefixes like `01_primary_exp_stage_lncRNA.csv`) so the file list on Zenodo stays readable. A **zip** is still useful when you want one download that recreates the same folder tree after unpack.

**Practical mix:** ship the big matrices and cohort tables as **individual files**; put deeply nested optional content in **a small extra zip** if file count approaches 100.

---

## Folder layout (after running `scripts/build_zenodo_dataset_bundle.ps1`)

| Folder | Contents | Typical destination in a clone |
|--------|------------|----------------------------------|
| **`01_tcga_lncrna_expression/`** | `primary_exp_stage_lncRNA.csv`, `primary_exp_metastasis_lncRNA.csv`, optional `primary_exp_*_final.csv` | `paper-github/data/` |
| **`02_smprot_raw/`** | `SmProt2.txt` (if present at source) | `paper-github/data/` |
| **`03_netmhc_merged/`** | `*_with_iedb.tsv` merged cohort tables | `paper-github/data/netmhc/` |
| **`04_netmhc_wide_xls/`** | Cohort wide NetMHCpan `*.xls` (pre-merge) | `paper-github/data/netmhc/` |
| **`05_iedb_inputs/`** | `iedb_*.csv`, `iedb_*.json` used by merge scripts | `paper-github/data/netmhc/` |
| **`06_predictions_optional/`** | Optional TTN / API prediction subtrees (`predictions/…`) | `paper-github/data/netmhc/predictions/` (preserve inner folders) |
| **`07_provenance/`** | `MANIFEST_COPIED.txt` and `FILE_MANIFEST.txt` from the build | `paper-github/data/` (or keep as documentation only) |
| **`08_curated_lists_optional/`** | e.g. `lncrna_genes_small.csv` and other small lists not on GitHub | `paper-github/data/` |

Exact filenames may vary; see **`FILE_MANIFEST.txt`** in the bundle root for what was actually copied.

---

## How to use this archive for reproduction

1. Clone or download the **code** repository (`paper-github`). Install Python/R dependencies per root `README.md`.  
2. Obtain this dataset: **unzip** a bundle **or** download only the Zenodo files you need.  
3. Copy files into the clone so paths match what `repo_paths.py` expects (see table above). On Windows, you can drag folders or use:

   ```powershell
   $clone = "C:\path\to\paper-github"
   $zen   = "C:\path\to\unzipped_bundle"
   Copy-Item "$zen\01_tcga_lncrna_expression\*"   "$clone\data\" -Force
   Copy-Item "$zen\02_smprot_raw\*"               "$clone\data\" -Force -ErrorAction SilentlyContinue
   Copy-Item "$zen\03_netmhc_merged\*"            "$clone\data\netmhc\" -Force -ErrorAction SilentlyContinue
   Copy-Item "$zen\04_netmhc_wide_xls\*"          "$clone\data\netmhc\" -Force -ErrorAction SilentlyContinue
   Copy-Item "$zen\05_iedb_inputs\*"              "$clone\data\netmhc\" -Force -ErrorAction SilentlyContinue
   if (Test-Path "$zen\06_predictions_optional") {
     Copy-Item "$zen\06_predictions_optional\*" "$clone\data\netmhc\predictions\" -Recurse -Force -ErrorAction SilentlyContinue
   }
   ```

4. Re-run orchestrators with `--strict` as in `docs/CODE_REVIEW.md` / `README.md` once inputs are in place.

---

## Citation

Cite this dataset’s **Zenodo DOI** in the manuscript **Data availability** section, together with the **software** DOI (GitHub / code Zenodo record) if you split them.

---

## Build notes (for maintainers)

The bundle is generated locally by:

**Windows (PowerShell or cmd):**

```powershell
cd paper-github
.\scripts\build_zenodo_dataset_bundle.ps1 -DataParent "C:\path\to\UNDEFINED" -Zip
```

**Git Bash** (use Unix-style paths, not `c:\...`):

```bash
cd /c/Users/you/Desktop/masters/UNDEFINED/paper-github
./scripts/build_zenodo_dataset_bundle.sh /c/Users/you/Desktop/masters/UNDEFINED --zip
```

**WSL:**

```bash
cd /mnt/c/Users/you/Desktop/masters/UNDEFINED/paper-github
./scripts/build_zenodo_dataset_bundle.sh /mnt/c/Users/you/Desktop/masters/UNDEFINED --zip
```

`-DataParent` must be a directory that contains a **`data`** folder (e.g. your full project root where `data\primary_exp_*.csv` lives). The script skips missing files and lists what was copied in `FILE_MANIFEST.txt`.
