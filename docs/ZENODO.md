# Archiving this repository on Zenodo

Zenodo ([zenodo.org](https://zenodo.org)) is a long-term archive operated by CERN. It gives you a **citable DOI** for a specific snapshot of your code + data (often linked to a **GitHub release**).

---

## What goes where (recommended split)

Think in **two deposits**: (A) **GitHub release → Zenodo “software”** = everything **git-committed** in `paper-github/`; (B) **Zenodo “dataset”** (manual upload) = large tables that are **gitignored** so the repo stays clone-friendly.

### A — GitHub (and therefore the auto-archived GitHub→Zenodo zip)

**Commit and push** anything needed to **re-run scripts** and **reproduce figures** that you are happy to version in git:

| Area | Typical contents |
|------|-------------------|
| **Code** | `manuscript/`, `supplement/`, `pipeline/`, `scripts/`, `*.py`, `*.R`, orchestrators |
| **Docs** | `docs/`, `data/README.md`, root `README.md` |
| **Small / bundled inputs** | `data/known_proteins.fasta`, `data/smprot_*.tsv` / `.faa` (if sizes OK), `data/netmhc/hla_european27_*.csv` / `.txt`, `data/canonical_significant_lncRNAs.txt`, TTN helper FASTAs, bundled **`netmhcpan_ttn_as1_108065.xls`** when small |
| **Expression + DE** | `tr_lncrna_output/**` — OK on GitHub if totals stay under limits. **`data/primary_exp_stage_lncRNA.csv`** and **`data/primary_exp_metastasis_lncRNA.csv`** are usually **>100 MiB** → **do not commit** (see `.gitignore`); ship on the **Zenodo dataset** instead (or use **Git LFS** if you insist on GitHub). |
| **Figures** | Canonical main-text PNGs at `figures/` root + `figures/supplementary/`; `data/netmhc/figures/*.png`. **Not** on GitHub: `figures/external/`, `figures/publication/`, `figures/paper_tifs/`, ECCB source deck, combined manifest, Fig 1B single-panel PNGs (see `.gitignore`) |

**Do not expect** the following on GitHub: merged NetMHC TSVs, cohort wide XLS, IEDB exports, **`SmProt2.txt`**, **`primary_exp_*_final.csv`**, **`primary_exp_*_lncRNA.csv`** (when large) — they are **gitignored** by design; put them on the **Zenodo dataset** (see below). Smaller **`data/`** files that remain under the GitHub size limit can stay on GitHub per `data/README.md`.

### B — Zenodo **dataset** record (manual; second DOI) — **your chosen bundle**

**Folder layout + README for upload:** see **`zenodo_dataset_bundle/README.md`** (human-facing text for the zip root) and run **`scripts/build_zenodo_dataset_bundle.ps1`** to assemble numbered folders (`01_tcga_lncrna_expression/`, `02_smprot_raw/`, …) from your machine’s **`data/`** tree before zipping.

Upload one **zip** (or multiple files) with everything you want archived **but not on GitHub**. This repo’s policy for your manuscript is:

| Path (relative to project `data/`) | Include on Zenodo? | Role |
|------------------------------------|--------------------|------|
| **`SmProt2.txt`** | **Yes** | Raw SmProt export for curation / `build_significant_lncs_smprot.py` |
| **`primary_exp_stage_lncRNA.csv`**, **`primary_exp_metastasis_lncRNA.csv`** | **Yes** | TCGA lncRNA matrices for Fig 1B / `tr_lncrna_de_analysis.py` / `tr_limma_de.R` (too large for GitHub) |
| **`primary_exp_stage_final.csv`**, **`primary_exp_metastasis_final.csv`** (and any other **`primary_exp_*_final.csv`**) | **Yes** | Full / wide expression matrices after staging |
| **`netmhc/netmhcpan_sig_lnc_with_iedb.tsv`** | **Yes** | Merged cohort — sig lnc |
| **`netmhc/netmhcpan_coding_proportional_whole_with_iedb.tsv`** | **Yes** | Merged coding — proportional whole |
| **`netmhc/netmhcpan_coding_control_with_iedb.tsv`** | **Yes** | Merged coding — control |
| **`netmhc/netmhcpan_sig_lnc.xls`**, **`netmhc/netmhcpan_coding_proportional_whole.xls`**, **`netmhc/netmhcpan_coding_control.xls`** (and any other cohort **`.xls`**) | **Yes** | Wide NetMHCpan runs before merge |
| **`netmhc/iedb_*.csv`**, **`netmhc/iedb_*.json`** | **Yes** | IEDB peptide_table / side inputs used in merge |
| **`netmhc/predictions/**/netmhcpan*.xls`** (TTN API / local API runs) | **Yes** (if you use them) | Optional TTN / API wide tables |
| **`lncrna_genes_small.csv`**, other SmProt-side TSVs **not** on GitHub | **As needed** | Curated gene lists, etc. |
| **`MANIFEST_COPIED.txt`** or **`SHA256SUMS.txt`** | **Yes** | Provenance / checksums |

**Not** required on Zenodo if they remain **small and committed** on GitHub: filtered **`smprot_*.tsv`/`*.faa`**, **`known_proteins.fasta`**, allele frequency CSV, canonical gene list — unless you later choose to move those to Zenodo too. The two **`primary_exp_*_lncRNA.csv`** matrices are **required on Zenodo** when they are gitignored from GitHub.

**License on the dataset record:** use Zenodo’s picker (e.g. **CC BY 4.0** for data-only) if it should differ from code **MIT**.

---

## What you need (accounts & metadata)

1. **A Zenodo account** (free; “Sign in with GitHub” is easiest).
2. **A coherent snapshot** of what you want archived:
   - **GitHub path:** tag a **release** on the public repo (Zenodo’s GitHub integration zips that release automatically).
   - **Large files:** anything **gitignored** in this repo (merged NetMHC TSVs, wide XLS, IEDB exports) is **not** in the GitHub zip. Upload those as a **separate Zenodo record** or a **Zenodo “version”** with multiple files, and cite **both** DOIs in the paper (code + data).
3. **Metadata** Zenodo will ask for: title, authors, description, license, keywords, related paper (optional). Match the **`LICENSE`** file in the repo (MIT here unless you change it).

## GitHub → Zenodo (recommended)

1. Push **`paper-github/`** to GitHub as its own repository (or your chosen org repo).
2. On Zenodo: **Account → GitHub** → enable the repository → flip the Zenodo archive switch.
3. On GitHub: **Releases → Create a new release**, pick a tag (e.g. `v1.0.0`), add release notes, publish.
4. Zenodo ingests the release and assigns a **DOI** (reserve DOI in advance if you need it for the manuscript before public release).

## Manual upload (no GitHub integration)

1. **Zenodo → New upload**.
2. Upload a **zip** you built locally (include `figures/` and `figures/manuscript_netmhc/` if you want figures in the archive even when not on GitHub).
3. Fill metadata, choose **license**, publish → get DOI.

## Optional niceties

- **`.zenodo.json`** in the repo root can pre-fill Zenodo metadata (title, creators, license) for GitHub-triggered uploads.
- **`CITATION.cff`** helps GitHub show a “Cite this repository” button; Zenodo can also read it in some workflows.

## Updating “a folder” later

- **In Git:** change files under e.g. `figures/manuscript_netmhc/`, commit, tag a **new** release (`v1.0.1`) → Zenodo creates a **new version** with a new DOI (old DOI remains valid).
- **From terminal into the repo:** use `robocopy` / `rsync` / `cp -r` from your full analysis tree into `paper-github/figures/...`, then `git add` and commit (same as we did for `figures/manuscript_netmhc/`).

---

## GitHub: what to add, commit, and push (checklist)

Run from **`paper-github/`** (or your clone root). Adjust branch name (`main` vs `master`).

```bash
git status
# Track everything you intend to ship (code, docs, LICENSE, figures, data that is NOT gitignored):
git add LICENSE README.md requirements.txt install_r_dependencies.R tr_limma_de.R repo_paths.py
git add manuscript/ supplement/ pipeline/ scripts/ notebooks/ docs/
git add generate_*.py rebuild_*.py
git add data/ figures/
git add tr_lncrna_output/   # only if you want DE/z outputs in the repo and they are not huge
git status   # confirm no accidental huge gitignored files listed as "deleted" etc.
git commit -m "Release bundle: code, docs, figures, bundled data"
git remote add origin https://github.com/<you>/<repo>.git   # skip if already set
git push -u origin main
```

Then **GitHub → Releases → Draft a new release → choose tag `v1.0.0` → Publish**.

**If GitHub rejects a large file:** use **[Git LFS](https://git-lfs.com/)** for that path, or **remove** it from git and put it only on the **Zenodo dataset** record instead.

## Zenodo dataset (second DOI): build the zip from your full machine

Preferred: from **`paper-github/`** run the bundler.

**PowerShell (Windows):**

```powershell
cd C:\path\to\paper-github
.\scripts\build_zenodo_dataset_bundle.ps1 -DataParent "C:\path\to\UNDEFINED" -Zip
```

**Git Bash or WSL** (bash does not run `.ps1` files; use the wrapper and Unix paths):

```bash
cd /c/Users/you/Desktop/masters/UNDEFINED/paper-github   # Git Bash — note /c/ not c:\
./scripts/build_zenodo_dataset_bundle.sh /c/Users/you/Desktop/masters/UNDEFINED --zip
```

```bash
cd /mnt/c/Users/you/Desktop/masters/UNDEFINED/paper-github   # WSL
./scripts/build_zenodo_dataset_bundle.sh /mnt/c/Users/you/Desktop/masters/UNDEFINED --zip
```

This writes **`zenodo_dataset_staging/`** under `paper-github/` (gitignored) with **`README.md`**, **`FILE_MANIFEST.txt`**, and numbered folders, then **`zenodo_dataset_staging.zip`** when `-Zip` / `--zip` is used. Adjust `-DataParent` if your `data/` lives directly under `paper-github` (then pass the `paper-github` directory as `-DataParent`).

Upload the **zip** on Zenodo → **New upload** → resource type **Dataset** → publish.

**Manual alternative (flat list):** if you prefer not to use the script, you can still `Compress-Archive` a hand-picked file list; the table above remains the checklist.

**If a file was ever committed to GitHub before adding `.gitignore`:** remove it from git tracking (keep local copy):  
`git rm --cached data/SmProt2.txt` then commit, so the public repo no longer grows with that blob.

In the paper **Data availability**, cite **both** DOIs: **software** (GitHub release) + **dataset** (this zip).

