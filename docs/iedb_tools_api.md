# IEDB Tools API — TTN-AS1 NetMHCpan 4.1 fetch

This repo calls the **public IEDB MHC class I binding prediction** endpoint from `supplement/fetch_ttn_mhci_iedb_api_netmhc41.py` (no local NetMHCpan binary).

**Local NetMHCpan instead (FASTA → real `*.xls`):** that is a **separate** workflow (WSL binary, `-xls 1 -xlsfile …`). See **`docs/figure6_ttn_as1_parameters.md`** (section *Local NetMHCpan*) and **`data/netmhc/README_netmhc.md`** §5 — the checked-in example is **`data/netmhc/run_netmhcpan_ttn_as1_108065.sh`**. In a full **UNDEFINED** checkout the same paths usually live under **`data/netmhc/`** at the tree root as well as inside **`paper-github/`**.

## Endpoint

Default URL (overridable with `--url`):

`https://tools-cluster-interface.iedb.org/tools_api/mhci/`

## HTTP method and body

- **Method:** `POST`
- **Content-Type:** `application/x-www-form-urlencoded`
- **Form fields** (as implemented in the script):

| Field | Value in this workflow |
|-------|-------------------------|
| `method` | `netmhcpan_ba-4.1` (binding / IC50) or `netmhcpan_el-4.1` (eluted ligand score) |
| `allele` | IEDB token, e.g. `HLA-A*01:01` (script maps from `HLA-A01:01` lines in `data/netmhc/hla_european27_class1.txt`) |
| `length` | `9` (9-mer peptides) |
| `sequence_text` | Peptide **FASTA** string: headers `>id` and one sequence line per peptide (see `--fasta`) |

The script parses the returned **tab-separated** table (`csv.DictReader` with `delimiter="\t"`). Column names used downstream include `ic50`, `percentile_rank` (as BA/EL rank), and EL `score`.

**Not sent in this POST:** immunogenicity scores, antigen-processing scores, or EL/IC50 *thresholds*. Those are **not** knobs on the `mhci/` binding call in our code; they show up later only if you merge to an **IEDB peptide_table** export and run Fig 6 with `--gating iedb_sb` (see `docs/figure6_ttn_as1_parameters.md`).

## CLI parameters (`python supplement/fetch_ttn_mhci_iedb_api_netmhc41.py --help`)

| Flag | Default | Role |
|------|---------|------|
| `--out-dir` | `data/netmhc/predictions/ttn_as1_smpep108065_iedb_api_netmhc41` | Raw per-allele TSV logs + synthetic wide `*.xls` |
| `--fasta` | `data/netmhc/ttn_as1_108065_ninemers.fasta` | Input peptides |
| `--alleles` | `data/netmhc/hla_european27_class1.txt` | Alleles to query (one per line) |
| `--url` | see above | API base URL |
| `--sleep` | `1.0` | Seconds between alleles / batches (politeness) |
| `--batch-size` | `20` | Peptides per POST (`0` = single POST with all peptides; large payloads can trigger HTTP 500) |
| `--max-retries` | `3` | Retries on HTTP 5xx |
| `--retry-backoff` | `2.0` | Initial backoff seconds before first retry |
| `--insecure` | off | Disable TLS verification (Windows revocation issues) |

## Synthetic wide XLS

The script merges BA + EL rows per allele into a **NetMHCpan-style wide workbook** (see module docstring for BA_score reconstruction from IC50). Plotting, SB thresholds, and optional peptide-table merge for Fig 6 are documented separately in **`docs/figure6_ttn_as1_parameters.md`**.

---

## Appendix: cohort wide XLS — local NetMHCpan-4.2 (shell defaults)

**Not IEDB:** the rows below are the **checked-in bash** defaults used to generate **cohort** wide `*.xls` files (Fig 5 wide / merge inputs). Same **prediction flags** are used in **`data/netmhc/run_netmhcpan42_example.sh`** and **`data/netmhc/run_netmhcpan_ttn_as1_108065.sh`** (TTN only differs by FASTA and `-xlsfile` path).

### Environment and paths (both scripts)

| Item | Default / value |
|------|------------------|
| `NETMHCpan` | `$HOME/netMHCpan-4.2` (override by exporting before run) |
| Binary | `$NETMHCpan/Linux_x86_64/bin/netMHCpan-4.2` (override with `NETMHCPAN=…`) |
| Allele list file | `data/netmhc/hla_european27_class1.txt` (override with `ALLELES=…`) |
| Alleles on CLI | Comma-separated contents of that file (`tr -d '\r' \| … \| paste -sd, -`) — **not** a path passed to `-a` |
| `TMPDIR` | `$HOME/tmp` (created if missing; avoids WSL `/mnt/c` tmp failures) |
| `bin` → `Linux_x86_64/bin` symlink | Created if missing so helper binaries resolve |

### NetMHCpan prediction flags (`FLAGS` in shell)

| Flag | Value | Role |
|------|-------|------|
| `-inptype` | `0` | Peptide input (FASTA) |
| `-a` | `"$ALLELE_ARG"` | All alleles in one run (comma-separated) |
| `-l` | `9` | Peptide length 9 |
| `-BA` | `1` | Include binding affinity / IC50-style BA output |
| `-pathogen` | `1` | Pathogen / IEDB-fine-tuned head (4.2; omit on 4.1 if unsupported) |
| `-neo` | `1` | Neoepitope / CEDAR-fine-tuned head (4.2; omit on 4.1 if unsupported) |
| `-context` | `0` | Context encoding off for these FASTA runs |
| `-xls` | `1` | Wide XLS-style tab-separated output |
| `-f` | FASTA path | Input ninemers |
| `-xlsfile` | output `.xls` | Destination wide table |

### Cohort example (`run_netmhcpan42_example.sh`)

| Step | Path |
|------|------|
| Input (sig lnc) | `data/netmhc/ninemers_sig_lnc.fasta` |
| Output | `data/netmhc/netmhcpan_sig_lnc.xls` |
| Input (coding control) | `data/netmhc/ninemers_coding_control.fasta` |
| Output | `data/netmhc/netmhcpan_coding_control.xls` |

**Proportional-whole coding cohort** uses the **same flag block** in **`data/netmhc/README_netmhc.md`** (§5d / §5 examples): point `-f` at `data/netmhc/ninemers_coding_proportional_whole.fasta` and `-xlsfile` at `data/netmhc/netmhcpan_coding_proportional_whole.xls`.

### TTN single-peptide example (`run_netmhcpan_ttn_as1_108065.sh`)

| Input | `data/netmhc/ttn_as1_108065_ninemers.fasta` |
| Output | `data/netmhc/netmhcpan_ttn_as1_108065.xls` |

For install, `-xlsfile` vs deprecated `-out`, and 4.1 vs 4.2 caveats, see **`data/netmhc/README_netmhc.md`** §5.

## Official IEDB documentation

For field semantics and rate limits, refer to the **IEDB Tools API** documentation on [iedb.org](https://www.iedb.org/) (MHC-I binding prediction tools cluster). This file only documents what **this repository** sends and expects.
