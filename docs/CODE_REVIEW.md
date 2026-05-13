# Code review — `paper-github` bundle

Living audit for **reproducibility, safety, and maintainability**. Last pass used automated scans (paths, subprocess, bare `except`, `shell=True`, TODOs) plus spot-checks of orchestrators and shared libs.

## Inventory

| Area | Count | Notes |
|------|-------|--------|
| **`manuscript/`** | 13 × `.py` | Default figure entrypoints |
| **`supplement/`** | 19 × `.py` | Sensitivity, legacy wide cohort, IEDB helpers, shims |
| **`pipeline/`** | 15 × `.py` | SmProt / TCGA / NetMHC prep |
| **`scripts/`** | 2 × `.py` | Merge + SB core (no `repo_paths`; CLI-oriented) |
| **Root** | 4 × `.py` | Orchestrators + `repo_paths.py` + `figure_palettes.py` |
| **R** | `tr_limma_de.R` | `ROOT` from argv; mirrors Python filter constants (see drift risk) |
| **Shell** | `data/netmhc/*.sh` | WSL NetMHCpan examples |

## Automated scan (summary)

| Check | Result |
|-------|--------|
| Hardcoded user paths (`C:\\Users`, `/home/.../UNDEFINED`, etc.) in `*.py` | **None** found |
| `subprocess` with `shell=True` | **None** |
| Bare `except:` | **None** (line-only grep) |
| `eval(` / `pickle.load` | **None** in scanned `*.py` |
| `TODO` / `FIXME` | **None** in `*.py` |

## Path hygiene

- **Good:** Almost every runnable script bootstraps `sys.path` and imports **`repo_paths`** (`REPO_ROOT`, `DATA`, `NETMHC_DATA`, …). `repo_paths.py` anchors on `Path(__file__).parent`.
- **Exceptions (acceptable):** `figure_palettes.py` (constants only); `scripts/merge_netmhcpan_xls_with_iedb.py` / `scripts/netmhc_sb_core.py` (library/CLI, paths via arguments); `supplement/plot_stage_tsne_lncrna.py` (deprecated shim → subprocess to `manuscript/`).

## Subprocess / orchestration

| Pattern | Files | Risk |
|---------|-------|------|
| `subprocess.call` return value propagated to `SystemExit` (via `orchestrate_subprocess.call_echo`) | `generate_catalog_figures.py`, `generate_netmhc_figure_bundle.py`, `generate_netmhc_supplement.py`, `rebuild_netmhc_merged_tsvs.py`, `generate_tr_lncrna_identification.py` | **Low** — exit code surfaces |
| `subprocess.run` + explicit `returncode` check | `regenerate_manuscript_netmhc_figures.py`, `run_fig5_sig_vs_proportional_coding.py`, `plot_figure6_ttn_netmhc41_bundle.py`, `cluster_gene_peptides_cdhit.py` | **Low** |
| `subprocess.run` + **warn** on non-zero, continue (unless `--strict`) | `export_tcga_filtered_peptides_fasta.py` (chained significant FASTA + sync) | **Low** — default warns; `--strict` fails CI on partial chain |

## R ↔ Python duplication (maintainability)

`tr_limma_de.R` and `pipeline/tr_lncrna_de_analysis.py` both define **`LOG2_THRESH`**, **`EXPR_FRAC_MIN`**, **`MIN_SAMPLES_CANCER`**, transition lists, etc. **Drift risk:** change one without the other → paper inconsistency.

**Recommendation:** single JSON or small shared constants file generated from Python once, consumed by R; or a CI check that greps both files for the same literals. **Implemented:** `tests/test_tr_lncrna_constants_parity.py` compares Python module values to `tr_limma_de.R` assignments and transition names.

## Security / network

- **`supplement/fetch_ttn_mhci_iedb_api_netmhc41.py`:** public IEDB POST; `--insecure` documented for Windows TLS. No secrets in repo. **Rate / batching** already documented in `docs/iedb_tools_api.md`.

## Suggested follow-ups (priority)

| Priority | Action | Status |
|----------|--------|--------|
| **P0** | Before publication: run **`python generate_catalog_figures.py --strict`** and **`python generate_netmhc_figure_bundle.py --strict`** on a clean clone (+ Zenodo data unpacked) and capture logs as release evidence. | Manual (release checklist). |
| **P1** | **Smoke tests** (`pytest`): `manuscript/*.py`, `supplement/*.py`, and merge/SB scripts run **`--help`** — see `tests/test_smoke_cli.py`. | Done. |
| **P1** | **`--strict`** on `pipeline/export_tcga_filtered_peptides_fasta.py` chain (significant FASTA + sync) for non-zero exit on failure. | Done. |
| **P2** | Shared **`call_echo` / `run_echo`** in `orchestrate_subprocess.py` for orchestrators. | Done. |
| **P2** | **`mypy`** on `scripts/netmhc_sb_core.py` + `scripts/merge_netmhcpan_xls_with_iedb.py` — `pip install -r requirements-dev.txt` then `mypy scripts/netmhc_sb_core.py scripts/merge_netmhcpan_xls_with_iedb.py`. | Config in `pyproject.toml`. |
| **Maintainability** | R ↔ Python filter literals: `tests/test_tr_lncrna_constants_parity.py` asserts shared constants and transition names stay aligned. | Done. |

## Per-area reviewer checklist (manual)

Use this when touching a file in a PR:

1. **Paths:** only `repo_paths` / `argparse` paths / `Path(__file__).resolve().parent` — no new string literals under `data/` at repo root without `REPO_ROOT`.
2. **CLI:** defaults point at bundled small data; large inputs remain gitignored (see `.gitignore` + `docs/ZENODO.md`).
3. **Figures:** write under `FIGURES` or `NETMHC_FIGURES`; document if mirroring to repo `figures/`.
4. **Subprocess:** exit codes must propagate or be checked; document if warning-and-continue.
5. **Randomness:** any `--seed` documented in `analysis_params.md` / script docstring.

---

*Update this file when you complete a review pass or change global conventions.*
