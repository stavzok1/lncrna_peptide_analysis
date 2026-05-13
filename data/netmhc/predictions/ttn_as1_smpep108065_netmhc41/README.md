# TTN-AS1 (smPEP 108065) — NetMHCpan **4.1** predictions (BA + EL)

## Version vs IEDB

- **This folder is for local NetMHCpan-4.1 outputs**, not files downloaded from the IEDB web UI.
- The **IEDB MHC Class I binding / processing tools** that call NetMHCpan historically track a specific **NetMHCpan major version** (often **4.1** in older documentation; **4.2** appears in newer bundles). **IEDB does not stamp “4.1 vs 4.2” into your XLS** — you choose that by **which standalone tarball you install** and which binary you run.
- The canonical Fig 6 XLS in this repo was produced with **`run_netmhcpan_ttn_as1_108065.sh` → NetMHCpan-4.2**. For **4.1**, use **`run_netmhcpan41_ttn_as1_108065.sh`** in this directory and point Fig 6 at the new XLS.

## What gets predicted

- **Input:** `../../ttn_as1_108065_ninemers.fasta` (sliding 9-mers over the 79 aa TTN-AS1 parent).
- **Alleles:** `../../hla_european27_class1.txt` (27-allele panel; same as cohort scripts).
- **Outputs:** wide **`-xls 1`** table with **EL_score / EL_rank** and, with **`-BA 1`**, **BA_score / BA_rank** — the same column layout `plot_figure6_ttn_as1_allele_coverage.py` expects.

## Flags (“IEDB-style” where supported)

- **Always (4.1-safe):** `-inptype 0 -l 9 -BA 1 -xls 1` (FASTA in, 9-mer length, BA columns, XLS out).
- **Optional (4.2-style extras):** if you set `IEDB_STYLE_EXTRA_FLAGS="-pathogen 1 -neo 1 -context 0"` **and** your 4.1 build accepts them, export that before running. Many **4.1** installs **do not** implement `-pathogen`/`-neo`; if `netMHCpan -h` does not list them, **omit** them — **BA + EL still come from `-BA 1` + standard XLS**.

## Run (WSL / Linux)

From repo root:

```bash
export NETMHCpan="$HOME/netMHCpan-4.1"   # your unpacked 4.1 directory
bash data/netmhc/predictions/ttn_as1_smpep108065_netmhc41/run_netmhcpan41_ttn_as1_108065.sh
```

Writes:

- `netmhcpan_ttn_as1_108065_netmhc41.xls`
- `RUN_INFO.txt` (command, flags, timestamp)

## Fig 6 A–E from this XLS

From repo root (Windows PowerShell or WSL with same Python env):

```bash
python scripts/plot_figure6_ttn_netmhc41_bundle.py
```

Or call `plot_figure6_ttn_as1_allele_coverage.py` yourself with `--netmhc-xls` pointing at the XLS above and `-o` under `figures/...`.
