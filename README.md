# lncRNA micropeptides — analysis code (GitHub bundle)

This folder is a **trimmed copy** of the parent project: **orchestrators**, **figure scripts**, **NetMHC + IEDB merge utilities**, and **documentation**, plus **small tabular/FASTA inputs**. Large NetMHC / IEDB tables are **not** committed (see `data/README.md` and `.gitignore`).

## Setup

```bash
python -m venv .venv
.venv\Scripts\activate   # Windows
pip install -r requirements.txt
```

Run commands from this directory (repository root).

## Regenerate figures

| Command | Output |
|---------|--------|
| `python generate_catalog_figures.py` | Figures 2–4A (`figures/tcga_matrix/`, `figures/all_smprot_filtered/`, shared `figures/`) |
| `python generate_netmhc_figure_bundle.py` | Figures 5–6 NetMHC bundle (needs merged `*_with_iedb.tsv`; see `data/README.md`) |
| `python scripts/regenerate_manuscript_netmhc_figures.py --help` | Full clean matrix under `figures/manuscript_netmhc/` |

Documentation: `docs/figure_catalog.md`, `docs/netmhc_figure_commands.md`, `data/netmhc/README_netmhc.md`.

## R

`tr_limma_de.R` is included for the limma leg of the DE workflow; run in R as in your original pipeline.

## License

Use the same terms as your thesis / manuscript; add a `LICENSE` file when you publish the GitHub repo.
