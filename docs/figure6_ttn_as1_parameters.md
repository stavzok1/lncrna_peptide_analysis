# Figure 6 — TTN-AS1 allele coverage (NetMHCpan wide XLS + optional IEDB peptide table)

Manuscript script: `manuscript/plot_figure6_ttn_as1_allele_coverage.py` (CLI entry `manuscript/manuscript_figure6_ttn_as1.py`). Full flag list: `python manuscript/plot_figure6_ttn_as1_allele_coverage.py --help`.

## IEDB MHC-I binding API vs immunogenicity / processing

The helper `supplement/fetch_ttn_mhci_iedb_api_netmhc41.py` calls only the **MHC class I binding prediction** endpoint (`tools_api/mhci/`). In this repository the POST body contains **only** `method`, `allele`, `length`, and `sequence_text` (see `docs/iedb_tools_api.md`). There is **no** parameter in our client that tells that endpoint to “run processing” or to tune **immunogenicity** thresholds: NetMHCpan **BA** and **EL** methods return binding / elution-related columns only.

**Immunogenicity** (`iedb_score`) and **antigen processing** (`iedb_processing_score`) appear as columns after you **merge** NetMHC rows with an **IEDB peptide_table CSV** (`--iedb-csv` when `--gating iedb_sb`). Those scores are produced by IEDB’s **analysis / peptide resource** workflows (exported table), not by the same minimal POST used for BA+EL prediction in the fetch script. Defaults `--imm-min` and `--proc-min` are **downstream filters** on those CSV columns (see `scripts/netmhc_sb_core.py`), aligned with Fig 5 merged cohort logic.

**Summary:** API fetch → wide XLS (BA/EL ranks, IC50, etc.). Optional immuno/proc gates → **local** thresholds on merged peptide-table columns.

## Local NetMHCpan: FASTA → wide XLS (real binary output)

Use this when you want **NetMHCpan** itself (not the IEDB REST client) to score the same 9-mers and write a **wide XLS** that `plot_figure6_ttn_as1_allele_coverage.py` can read via `--netmhc-xls`. **Tabulated shell defaults** (`NETMHCpan`, `TMPDIR`, every `-` flag in `FLAGS`): **`docs/iedb_tools_api.md`** → *Appendix: cohort wide XLS — local NetMHCpan-4.2*. Full install, allele comma-list gotchas, and cohort FASTA prep: **`data/netmhc/README_netmhc.md`** (especially **§5**).

### TTN-AS1 smPEP 108065 (bundled example)

From the **`paper-github/`** repo root inside **WSL** (Linux paths):

```bash
bash data/netmhc/run_netmhcpan_ttn_as1_108065.sh
```

That script exports `NETMHCpan`, builds **`ALLELE_ARG`** from `data/netmhc/hla_european27_class1.txt`, then runs the binary on **`data/netmhc/ttn_as1_108065_ninemers.fasta`** with:

`-inptype 0 -a "$ALLELE_ARG" -l 9 -BA 1 -pathogen 1 -neo 1 -context 0 -xls 1 -xlsfile data/netmhc/netmhcpan_ttn_as1_108065.xls`

Output path matches the default **`--netmhc-xls`** in the Fig 6 script. If your **NetMHCpan-4.1** build does not accept `-pathogen`/`-neo`, omit flags your `netMHCpan -h` does not list (see **`data/netmhc/predictions/ttn_as1_smpep108065_netmhc41/README.md`**).

### Cohort wide XLS (Fig 5 inputs)

Example for significant lnc + coding-control ninemers:

```bash
bash data/netmhc/run_netmhcpan42_example.sh
```

(Writes `data/netmhc/netmhcpan_sig_lnc.xls` and `data/netmhc/netmhcpan_coding_control.xls` — see script headers.)

### Python wrapper (optional)

`supplement/plot_figure6_ttn_netmhc41_bundle.py` shells out to plot Fig 6 after you point **`--netmhc-xls`** at a local XLS.

### UNDEFINED vs `paper-github/`

If you work in the **parent** **`UNDEFINED`** tree (`…/masters/UNDEFINED`), you typically have **`data/netmhc/`** at that root too. The **same** `run_netmhcpan_*.sh` examples are maintained under **`paper-github/data/netmhc/`** for the GitHub bundle—run whichever copy sits next to the `data/netmhc/*.fasta` files you are using, from that tree’s repo root so relative paths resolve.

## Default numeric gates (IEDB SB path, Fig 5–aligned)

Constants in `scripts/netmhc_sb_core.py`; argparse defaults in `plot_figure6_ttn_as1_allele_coverage.py`:

| CLI flag | Default | Role |
|----------|---------|------|
| `--imm-min` | **0.1** | Require `iedb_score` above this when SB mode uses immuno |
| `--proc-min` | **1.5** | Require `iedb_processing_score` above this when SB mode uses processing |
| `--el-rank-max` | **1.0** | Compare to NetMHC **EL_rank** (strict less-than unless `--el-rank-lte`) |
| `--ic50-max-nm` | **150** | IC50 upper bound (nM) when the IC50 gate is on |

## Main CLI flags (coverage + SB stack)

| Flag | Default | Role |
|------|---------|------|
| `--gating` | `iedb_sb` | `netmhc`: SB from XLS ranks/IC50 only (no peptide table). `iedb_sb`: merge XLS to `--iedb-csv` then cohort-style SB via `netmhc_sb_core`. |
| `--netmhc-xls` | bundled TTN wide XLS | Local NetMHCpan output or synthetic XLS from `fetch_ttn_mhci_iedb_api_netmhc41.py`. |
| `--iedb-csv` / `--iedb-parent-input-seq-id` | synthetic companion | Required for `iedb_sb`. |
| `--sb-mode` | `full` | `full` — immuno + processing + EL + IC50/BA; `no_ic50` — drop binding gate; `ic50_only` — turn off immuno/proc/EL gates. |
| `--sb-criterion` | `ba_rank` | `ba_rank` vs `ic50` for binder definition. |
| `--ba-rank-pct` | `0.5` | SB if BA_rank ≤ this percentile (top 0.5%). |
| `--ic50-nm` | `150` | Used when `--sb-criterion ic50`. |
| `--require-el-rank` | off | Also require EL_rank cutoff. |
| `--coverage-output` | `instances` | `instances` / `unique` / `both`. |
| `--split-panels`, `--also-write-unique` | off | Multi-PNG layout; ``--also-write-unique`` adds ``*_unique_*`` next to ``*_instances_*`` in the **same** ``-o`` directory (orchestrators do **not** enable this by default). |

## Related

- **Fetch BA+EL via IEDB API (no local binary):** `docs/iedb_tools_api.md`, `supplement/fetch_ttn_mhci_iedb_api_netmhc41.py`.
- **Local NetMHCpan install + cohort prep:** `data/netmhc/README_netmhc.md`.
- **Commands cheat sheet:** `docs/netmhc_figure_commands.md` (Figure 6 block).
