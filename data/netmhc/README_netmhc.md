# NetMHCpan-4.1 (analyzed significant lncRNA peptides vs coding control)

This folder holds **HLA allele lists** and (after you run the prep script) **9-mer FASTA** inputs for an in-silico class I screen: sliding **9-mers** (step size 1) from each parent micropeptide, **NetMHCpan-4.1** in **9-mer** mode, over a **27-allele** European-oriented panel.

The **primary set** is all **`smPEP_ID`** in **`data/significant_lnc_peptides.tsv`** (the **analyzed** table: rows whose sequences appear in **`significant_lnc_peptides.faa`**; typically **~506** after excluding untranslatable SmProt coordinates). The **control** is the same number of parent peptides **per amino-acid length**, each sampled as a contiguous substring from a **human coding proteome FASTA** (default **`data/known_proteins.fasta`**, UniProt-style sequences in this repo).

---

## 1. Install NetMHCpan-4.1 on WSL2 (Ubuntu)

These steps assume **Windows 11** (or 10) with **WSL2** and **Ubuntu** from the Microsoft Store.

### 1a. Enable WSL and install Ubuntu

- Open **PowerShell as Administrator** and run: `wsl --install` (or install **Ubuntu** from the Store, then `wsl --set-default-version 2`).
- Reboot if prompted. Launch **Ubuntu**, create a UNIX user, then run `sudo apt update`.

### 1b. Packages NetMHCpan usually needs

```bash
sudo apt install -y build-essential gawk tcsh unzip wget
```

Some DTU tarballs also expect a C compiler and **Perl** (often already present). If `INSTALL` mentions extra tools, install those too.

### 1c. Download and unpack NetMHCpan-4.1

1. Open [DTU Health Tech — NetMHCpan-4.1](https://services.healthtech.dtu.dk/service.php?NetMHCpan-4.1) in a browser, complete **academic registration**, and download the **Linux** tarball.
2. In Ubuntu, copy the file into your home directory (e.g. from `/mnt/c/Users/.../Downloads/`).

```bash
cd ~
tar xzf netMHCpan-4.1*.tar.gz
cd netMHCpan-4.1*
```

### 1d. Compile and test

Follow **`README`** and **`INSTALL`** inside the unpacked folder (often `make` or a provided install script). Then:

```bash
./netMHCpan -h
```

If the program cannot find data files, set paths as described in the package README (some installs require editing a small settings file or exporting an environment variable).

### 1e. Run from your project (Windows files under `/mnt/c/...`)

```bash
cd /mnt/c/Users/<you>/Desktop/masters/UNDEFINED
```

Use Linux paths in all `netMHCpan` arguments. Large outputs are fine on the NTFS mount; if I/O is slow, copy the `data/netmhc/` inputs into `~/work` inside WSL and run there.

---

## 2. Allele list

- **`hla_european27_class1.txt`**: one allele per line (example panel). **Replace** with the exact list from your paper or supplement if needed.
- **NetMHCpan-4.2** expects names as in **`$NETMHCpan/data/MHC_pseudo.dat`** (first column). Human alleles use the **no-asterisk** form, e.g. **`HLA-A02:01`**, **`HLA-B07:02`** — **not** `HLA-A*02:01`. If you see `cannot be found in hla_pseudo list`, fix the spelling/format of each line (or pick names from **`$NETMHCpan/data/allelenames`** / **`netMHCpan -listMHC`** if your build supports it).
- **Line endings:** the allele file must be **Unix LF** only. If it was saved on Windows with **CRLF**, each line becomes `HLA-A01:01\r` and **no longer matches** the pseudo list. Run **`dos2unix data/netmhc/hla_european27_class1.txt`** in WSL, or rewrite the file with LF-only newlines.
- **`-a` syntax (4.2):** **`-a` takes comma-separated allele names**, not a path to a text file. If you pass **`data/netmhc/hla_european27_class1.txt`** (or **`file:...`**) literally, NetMHCpan treats that whole string as **one** allele ID and prints **`cannot be found in hla_pseudo list`**. Build the argument from your allele file, for example: **`ALLELE_ARG=$(tr -d '\r' < data/netmhc/hla_european27_class1.txt | grep -v '^[[:space:]]*$' | paste -sd, -)`** then **`-a "$ALLELE_ARG"`** (see **`run_netmhcpan42_example.sh`**).

---

## 3. Build the 549-peptide FASTA (GENCODE + translate)

`prepare_netmhc_tr_vs_coding_epitopes.py` expects **`data/significant_lnc_peptides.faa`**. Generate it with the same transcript-slice exporter used for TCGA-filtered peptides, pointing at **`data/significant_lnc_peptides.tsv`**:

```bash
# From repo root; set --transcripts-fa to your GENCODE transcript FASTA (.fa or .fa.gz)
python export_tcga_filtered_peptides_fasta.py \
  --peptides-tsv data/significant_lnc_peptides_full.tsv \
  --transcripts-fa /path/to/gencode.v49.transcripts.fa.gz \
  --ensembl-fallback \
  --min-aa-length 9

python sync_significant_lnc_peptides_with_fasta.py
```

`--min-aa-length 9` recovers more rows than the default **10** (short ORFs). **`data/significant_lnc_peptides.faa`** may still omit some full-list IDs; **`sync_significant_lnc_peptides_with_fasta.py`** rewrites **`data/significant_lnc_peptides.tsv`** to match the FASTA exactly.

By default this writes **`data/significant_lnc_peptides.faa`**. If many ENST versions are missing locally, keep **`--ensembl-fallback`** (uses Ensembl REST when the local slice is invalid or missing).

---

## 4. Prepare 9-mer FASTA (Python, repo root)

```bash
python prepare_netmhc_tr_vs_coding_epitopes.py --seed 1
```

### Coding controls: fragments (default) vs whole proteins

The pipeline matches **counts per parent length** to the significant set. How those coding sequences are taken from **`--coding-fa`** is controlled by **`--coding-control-mode`**:

| Mode | Meaning |
|------|--------|
| **`fragments`** (default) | For each required length **L**, draw a **contiguous substring of length L** from a **longer** proteome sequence (random start, uniform over valid starts). This is the usual choice for **short** lnc micropeptides: a canonical human SwissProt-style FASTA has almost **no** entries whose **entire** length is e.g. 14 aa. |
| **`whole_protein`** | Use only **complete FASTA records** whose cleaned length is **exactly L**. Statistically closest to “pick proteins of the same length,” but often **infeasible** for short MPs unless you supply a peptide/sORF/peptidome FASTA with many short full entries. Requires **`--max-proteins 0`** (full scan) in practice so every length bucket can fill. |

```bash
# Default: length-matched fragments from longer proteins
python prepare_netmhc_tr_vs_coding_epitopes.py --seed 1

# Optional: whole-FASTA-entry controls (strict exact length; may raise if pools are too small)
python prepare_netmhc_tr_vs_coding_epitopes.py --seed 1 --coding-control-mode whole_protein --max-proteins 0
```

After each run, **`data/netmhc/coding_control_sampling_mode.txt`** records **`fragments`** or **`whole_protein`**.

Optional:

- **`--coding-fa path/to/proteome.fasta`**: human coding reference (e.g. UniProt reviewed + isoforms). Default: **`data/known_proteins.fasta`**.
- **`--max-proteins N`**: for **fragment** mode, stop after scanning **N** proteome entries (default **0** = scan until every length bucket is filled). For **whole_protein** mode, **`0`** is recommended so the scan can find enough exact-length entries.

Writes under **`data/netmhc/`**:

| File | Description |
|------|-------------|
| `coding_control_sampling_mode.txt` | `fragments`, `whole_protein`, or `proportional_whole` (+ parameters when applicable). |
| `sig_mp_length_distribution.csv` | Length histogram of the analyzed parent peptides (see `significant_lnc_peptides.tsv`). |
| `sig_parent_micropeptides.csv` | Parent sequences aligned to the analyzed significant TSV. |
| `coding_control_parent_micropeptides.csv` | Coding parents for **fragments** / **whole_protein** only (not written for `proportional_whole`; use `coding_proportional_whole_parent_micropeptides.csv` instead). |
| `ninemers_sig_lnc.fasta` | One 9-mer per record (headers encode parent + offset). |
| `ninemers_coding_control.fasta` | Same for the **fragment** or **whole_protein** control set (not used for `proportional_whole`). |
| `ninemers_coding_proportional_whole.fasta` | **Only** when `--coding-control-mode proportional_whole`: 9-mers from whole-proteome parents (default path; override with `--coding-ninemers-out`). |
| `coding_proportional_whole_parent_micropeptides.csv` | Parent sequences for proportional-whole controls (separate from `coding_control_parent_micropeptides.csv`). |
| `coding_proportional_bin_summary.csv` | Per length-bin: reference counts/fractions, targets, sampled counts, `sample_fraction_total`, `delta_sample_minus_ref_fraction`. |
| `coding_proportional_vs_reference_report.md` | Overall TVD / L2 / Pearson r between reference and sampled bin distributions. |
| `ninemers_summary.csv` | Parent and 9-mer counts; includes `ninemers_fasta` column with output paths. |

### 4b. **Proportional whole proteins** (length-bin matched to reference)

Matches the **fraction of parent micropeptides per length bin** (equal-width bins, default **5 aa**) to the **reference** FASTA, sampling **entire** proteome records whose length falls in each bin. Does **not** overwrite **`ninemers_coding_control.fasta`**; it writes **`ninemers_coding_proportional_whole.fasta`** (or **`--coding-ninemers-out`**) plus **`coding_proportional_whole_parent_micropeptides.csv`**, **`coding_proportional_bin_summary.csv`**, and **`coding_proportional_vs_reference_report.md`**.

**Manuscript / significant cohort (~501 analyzed MPs):** use **defaults** — do **not** pass `--peptide-faa` / `--peptide-tsv` (they default to `data/significant_lnc_peptides.faa` + `data/significant_lnc_peptides.tsv`).

```bash
# Reference = significant lncRNA MPs (~501) vs proportional whole-proteome coding parents
python prepare_netmhc_tr_vs_coding_epitopes.py --coding-control-mode proportional_whole --max-proteins 0 --seed 1
```

**Optional — larger non-significant reference:** the TCGA-expression–filtered SmProt export (`data/smprot_tcga_filtered_peptides.faa` + `data/smprot_filtered_tcga_expr_genes.tsv`) is **~2.6k parents**, not the 501 significant set. Only use it when you explicitly want that broader length reference:

```bash
python prepare_netmhc_tr_vs_coding_epitopes.py --coding-control-mode proportional_whole --max-proteins 0 \
  --peptide-faa data/smprot_tcga_filtered_peptides.faa \
  --peptide-tsv data/smprot_filtered_tcga_expr_genes.tsv
```

Optional: **`--no-proportional-top-up`**, **`--bin-width-aa`**, **`--proportional-target-n`**, **`--coding-ninemers-out path/to/out.fasta`**.

### 4c. Histogram: **proteome sampling** vs significant micropeptides

After **`prepare_netmhc_tr_vs_coding_epitopes.py`**, this compares:

- **Significant** parent MPs — **`data/netmhc/sig_parent_micropeptides.csv`**
- **Proteome sampling** — the empirical distribution of the **length-matched coding control** parents in **`data/netmhc/coding_control_parent_micropeptides.csv`** (contiguous fragments actually drawn from the coding FASTA for NetMHC controls), *not* the length distribution of every full SwissProt entry.

```bash
python plot_proteome_vs_sig_mp_length_histogram.py
```

Writes **`data/netmhc/figures/proteome_sampling_vs_sig_mp_histograms.png`** (overlay length histogram + pooled amino-acid frequency bars), merged length-count CSV, AA frequency table, and a short report.

Optional **reference** panel — whole-proteome FASTA **entry** lengths (can be slow on large files):

```bash
python plot_proteome_vs_sig_mp_length_histogram.py --full-proteome-protein-lengths --max-proteome-records 50000
```

Other flags: **`--sig-csv`**, **`--control-csv`**, **`--proteome-fa`**, **`--out-dir`**, **`--proteome-xmax`**.

### 4d. **Figure 5 — significant vs proportional-whole coding** (separate folder)

After you have ``netmhcpan_coding_proportional_whole.xls`` and ``ninemers_coding_proportional_whole.fasta``:

```bash
python scripts/run_fig5_sig_vs_proportional_coding.py
```

Writes merged coding TSV ``data/netmhc/netmhcpan_coding_proportional_whole_with_iedb.tsv`` (synthetic IEDB pass columns; see ``merge_netmhcpan_xls_with_iedb.py --synthetic-iedb-pass``) and Fig 5 PNG/CSV under ``figures/manuscript_netmhc/fig5_sig_vs_proportional_whole/``. Options: ``--skip-merge``, ``--out-dir``, paths to XLS/FASTA/TSV.

---

## 5. Run NetMHCpan-4.1 / 4.2 (CLI)

**Quick reference (shell defaults used in `run_netmhcpan42_example.sh` / `run_netmhcpan_ttn_as1_108065.sh`):** `docs/iedb_tools_api.md` in the **`paper-github/`** bundle → section *Appendix: cohort wide XLS — local NetMHCpan-4.2*.

### 5a. Why the examples used `-BA` only

In NetMHCpan-4.x, the **default** output is centred on **antigen presentation (EL)** style scores (`Score_EL`, `%Rank_EL`, binder labels, etc.). The **`-BA`** switch does **not** replace that: it **adds** predicted **binding affinity** columns (IC50-style) on top of the EL output. So “BA only” in a command line was never “BA instead of everything”—it was “turn on affinity as well.”

For the **full** set of optional columns your tarball exposes, turn on the extra heads your build documents (4.2 help excerpt):

| Flag | Meaning (when set to `1`) |
|------|---------------------------|
| **`-BA 1`** | Include binding-affinity prediction. |
| **`-pathogen 1`** | Include pathogen-epitope (IEDB-fine-tuned) predictions. |
| **`-neo 1`** | Include neoepitope (CEDAR-fine-tuned) predictions. |
| **`-context 1`** | Use context encoding (for FASTA, context is derived from the parent sequence). |

Always confirm exact spelling with:

```bash
/path/to/netMHCpan-4.2 -h
```

### 5b. Set **`NETMHCpan`** (install root) — fixes `Unable to open $NETMHCpan/data/version`

The binary looks up data files under **`$NETMHCpan/data/`** (see `netMHCpan-4.2 -h`: **`-rdir`**, **`-version`**). The environment variable name is **`NETMHCpan`** (capital `NETMHC`, lowercase `pan`) — not `NETMHCPAN`.

From WSL, point it at the **top-level unpacked directory** (the folder that contains `data/version`):

```bash
export NETMHCpan="$HOME/netMHCpan-4.2"   # or: export NETMHCpan=/path/to/netMHCpan-4.2
```

Then invoke the binary (often under `Linux_x86_64/bin/`). Optional: **`export PATH="$NETMHCpan/Linux_x86_64/bin:$PATH"`** if you call `netMHCpan-4.2` by name.

### 5b2. **`sh: 1: .../netMHCpan-4.2/bin/nnalign_...: not found`** (and **`estimate_PCC: not found`**)

The core program calls helper binaries as **`$NETMHCpan/bin/...`**. In many Linux tarballs the executables actually live under **`$NETMHCpan/Linux_x86_64/bin/`**, and nothing creates **`$NETMHCpan/bin/`** unless you ran the full install step from the package **`README`/`INSTALL`**.

**One-time fix** (from WSL, adjust the path if your install root differs):

```bash
export NETMHCpan="$HOME/netMHCpan-4.2"
ln -sfn "$NETMHCpan/Linux_x86_64/bin" "$NETMHCpan/bin"
```

After that, **`ls "$NETMHCpan/bin/nnalign_gaps_pan_play_MA_wgt_MN_v2_two_outputs_context_allelelist_pboth_inv_allmhc_v2"`** should resolve. The repo script **`run_netmhcpan42_example.sh`** creates this symlink automatically when **`bin/`** is missing.

### 5c. **`Cannot make tmpdir. Exit`** (WSL + `/mnt/c/`)

NetMHCpan creates scratch space under **`$TMPDIR`** (see **`-tdir`** in `netMHCpan-4.2 -h`). If **`TMPDIR` is unset**, empty, or points somewhere creation fails, you get **`Cannot make tmpdir. Exit`**.

Working from **`/mnt/c/...`** (Windows drives) can also break temp-dir creation depending on WSL / DrvFs behaviour. Fix by forcing a temp root on the **Linux (ext4) side**:

```bash
mkdir -p "$HOME/tmp"
export TMPDIR="$HOME/tmp"
```

Then rerun NetMHCpan from the same shell. (Optional: add those two lines to `~/.bashrc`.) You can also pass an explicit **`-tdir /home/you/some_empty_dir`** after `mkdir -p` that directory, if you prefer not to touch global `TMPDIR`.

### 5d. NetMHCpan-4.2: **`-out` is invalid** — use **`-xlsfile`**

Standalone **4.2** rejects **`-out`**. For Excel-style output, use **`-xls 1`** and **`-xlsfile <path>`** (defaults to `NetMHCpan_out.xls` if you omit the filename).

Example (WSL, from repo root on `/mnt/c/.../UNDEFINED`):

```bash
export NETMHCpan="$HOME/netMHCpan-4.2"
mkdir -p "$HOME/tmp"
export TMPDIR="$HOME/tmp"
NETMHCPAN="$NETMHCpan/Linux_x86_64/bin/netMHCpan-4.2"
ALLELES=data/netmhc/hla_european27_class1.txt
ALLELE_ARG="$(tr -d '\r' < "$ALLELES" | grep -v '^[[:space:]]*$' | paste -sd, -)"

$NETMHCPAN -f data/netmhc/ninemers_sig_lnc.fasta -inptype 0 -a "$ALLELE_ARG" -l 9 \
  -BA 1 -pathogen 1 -neo 1 -context 0 \
  -xls 1 -xlsfile data/netmhc/netmhcpan_sig_lnc.xls

$NETMHCPAN -f data/netmhc/ninemers_coding_control.fasta -inptype 0 -a "$ALLELE_ARG" -l 9 \
  -BA 1 -pathogen 1 -neo 1 -context 0 \
  -xls 1 -xlsfile data/netmhc/netmhcpan_coding_control.xls

# Proportional-whole control (separate FASTA; does not replace ninemers_coding_control.fasta)
$NETMHCPAN -f data/netmhc/ninemers_coding_proportional_whole.fasta -inptype 0 -a "$ALLELE_ARG" -l 9 \
  -BA 1 -pathogen 1 -neo 1 -context 0 \
  -xls 1 -xlsfile data/netmhc/netmhcpan_coding_proportional_whole.xls
```

If you prefer **tabular stdout** instead of XLS, drop `-xls 1 -xlsfile ...` and redirect:  
`$NETMHCPAN ... > data/netmhc/netmhcpan_sig_lnc.tsv` (column layout is documented in the package).

### 5e. Full copy-paste (WSL): one block to run NetMHCpan-4.2 from this repo

From **Ubuntu (WSL)**, paste the whole block after **`prepare_netmhc_tr_vs_coding_epitopes.py`** has created the FASTA files under **`data/netmhc/`**. Change **`cd`** to your actual repo path (Windows Desktop clone is usually under **`/mnt/c/Users/<you>/.../UNDEFINED`**).

```bash
cd /mnt/c/Users/<you>/Desktop/masters/UNDEFINED

export NETMHCpan="${NETMHCpan:-$HOME/netMHCpan-4.2}"
mkdir -p "$HOME/tmp"
export TMPDIR="$HOME/tmp"

# One-time per install: NetMHCpan shells out to $NETMHCpan/bin/* (see §5b2).
# The example script creates this symlink if missing; including it here makes the block self-contained.
if [[ -d "$NETMHCpan/Linux_x86_64/bin" && ! -d "$NETMHCpan/bin" ]]; then
  ln -sfn "$NETMHCpan/Linux_x86_64/bin" "$NETMHCpan/bin"
fi

bash data/netmhc/run_netmhcpan42_example.sh
```

**What `bash data/netmhc/run_netmhcpan42_example.sh` does (read this for methods text):**

- It runs **NetMHCpan-4.2 twice**, not once.
- **First run — significant lncRNA-derived 9-mers:** input **`data/netmhc/ninemers_sig_lnc.fasta`** (sliding 9-mers from the **analyzed** significant lnc micropeptide set prepared by **`prepare_netmhc_tr_vs_coding_epitopes.py`**, aligned to **`data/significant_lnc_peptides.tsv`** / **`significant_lnc_peptides.faa`**). Output **`data/netmhc/netmhcpan_sig_lnc.xls`**.
- **Second run — length-matched coding control 9-mers:** input **`data/netmhc/ninemers_coding_control.fasta`**. Output **`data/netmhc/netmhcpan_coding_control.xls`**.

So the script is **not** “significant epitopes only”: it is **significant lnc 9-mers plus the paired coding control** in one go, same allele string and same flags for both.

**Equivalent explicit 4.2 invocations** (same as the script; use only if you prefer not to call the shell script):

```bash
cd /mnt/c/Users/<you>/Desktop/masters/UNDEFINED

export NETMHCpan="${NETMHCpan:-$HOME/netMHCpan-4.2}"
mkdir -p "$HOME/tmp"
export TMPDIR="$HOME/tmp"
if [[ -d "$NETMHCpan/Linux_x86_64/bin" && ! -d "$NETMHCpan/bin" ]]; then
  ln -sfn "$NETMHCpan/Linux_x86_64/bin" "$NETMHCpan/bin"
fi

NETMHCPAN="$NETMHCpan/Linux_x86_64/bin/netMHCpan-4.2"
ALLELES="data/netmhc/hla_european27_class1.txt"
ALLELE_ARG="$(tr -d '\r' < "$ALLELES" | grep -v '^[[:space:]]*$' | paste -sd, -)"

"$NETMHCPAN" -f data/netmhc/ninemers_sig_lnc.fasta -inptype 0 -a "$ALLELE_ARG" -l 9 \
  -BA 1 -pathogen 1 -neo 1 -context 0 -xls 1 -xlsfile data/netmhc/netmhcpan_sig_lnc.xls

"$NETMHCPAN" -f data/netmhc/ninemers_coding_control.fasta -inptype 0 -a "$ALLELE_ARG" -l 9 \
  -BA 1 -pathogen 1 -neo 1 -context 0 -xls 1 -xlsfile data/netmhc/netmhcpan_coding_control.xls
```

### 5f. NetMHCpan-4.1 wrapper (if your install still uses `-out`)

Older wrappers sometimes accepted **`-out`**. If **your** 4.1 binary accepts it, keep your local command; otherwise use **`-xlsfile`** the same way as 4.2 after checking **`netMHCpan -h`**.

**Notes**

- `-inptype 0` = **FASTA** input in 4.2 help (`0` FASTA, `1` PEPTIDE, …).
- `-l 9` = peptide length 9 for FASTA digestion (here your FASTA is already 9-mers; length still matches).
- Very large FASTA files may need splitting; see the package docs for threads / memory (`-t` threshold etc.).

### 5g. Monitoring progress during a long NetMHCpan run

NetMHCpan **does not** expose a standard **percentage-complete** progress bar for huge FASTA jobs. Practical options:

| Approach | What to do |
|----------|------------|
| **Verbose flag** | Add **`-v 1`** to the `netMHCpan-4.2` command line (see **`netMHCpan-4.2 -h`**: *Verbose mode*). Output volume depends on the build; it may still be sparse per peptide. |
| **Log file** | Run with **`tee`**: `bash data/netmhc/run_netmhcpan42_example.sh 2>&1 \| tee data/netmhc/netmhcpan_run.log` so you can scroll the log and confirm the process is still writing. |
| **Output / temp** | In a second terminal, **`ls -lh data/netmhc/netmhcpan_*.xls`** or watch **`$TMPDIR`** for new **`netMHCpan_*`** scratch dirs while the job runs. Some builds buffer XLS until the end, so file size may jump late. |
| **`-dirty 1`** | Keeps temporary files/dirs for inspection (see **`-h`**); useful for debugging, not for routine huge runs. |
| **Chunked runs** | Split **`ninemers_*.fasta`** into smaller files, run NetMHCpan per chunk, and concatenate outputs if you need **incremental** results. |

---

## 5b. Figure catalog, cohort vs TTN-AS1, and bundles

Manuscript-style numbering for **NetMHC-related** panels is documented in
**`docs/figure_catalog.md`** (sections **Figure 5** and **Figure 6**).

- **Figure 5** = **cohort** (many 9-mers: significant lncRNA MPs vs coding control): allele–epitope
  interplay, wide XLS scripts at repo root, optional **IEDB-merged** pipeline under `scripts/`.
- **Figure 6** = **single peptide** **TTN-AS1** (smPEP 108065): `manuscript/manuscript_figure6_ttn_as1.py`
  (→ `plot_figure6_ttn_as1_allele_coverage.py`; default PNG under repo `figures/`), NetMHC-only sweeps
  `plot_figure6_ttn_as1_sb_sensitivity.py`, merged IEDB **1D+LOO** `netmhc_ttn_merged_iedb_sb_sensitivity_robustness.py`,
  and merged IEDB **Cartesian** grid `plot_fig6_ttn_merged_iedb_sb_combination_grid.py` (defaults under
  `data/netmhc/figures/...`; orchestrated copies under `figures/supplementary/netmhc_fig5_fig6_supplement/...`).

To regenerate the **canonical** NetMHC figure set (merged 5–6 + optional random-fragment mirrors):

```bash
python generate_netmhc_figure_bundle.py
```

For **Fig 5–6 supplement** (1D + LOO cohort, Cartesian Fig 5, TTN NetMHC sweeps, TTN merged IEDB 1D+LOO, TTN merged IEDB Cartesian):

```bash
python generate_netmhc_fig5_fig6_supplement.py
```

Full orchestration map: **`docs/figure_generation_overview.md`**.

Agents should read **`.cursor/skills/netmhc-manuscript-figures/SKILL.md`** for the same map.

---

## 6. Provenance

Document: **NetMHCpan version**, **command line**, **allele file source**, **`--seed`**, **GENCODE (or Ensembl) release** used for **`significant_lnc_peptides.faa`**, and **proteome file** used for **`--coding-fa`**.
