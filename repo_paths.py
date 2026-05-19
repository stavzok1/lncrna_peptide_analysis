"""
Single source for repository-relative paths. All scripts derive locations from here
(no absolute paths to user directories).

Import after putting the repo root on ``sys.path`` (see any script under
``manuscript/``, ``supplement/``, or ``pipeline/``).
"""
from __future__ import annotations

from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent
DATA = REPO_ROOT / "data"
FIGURES = REPO_ROOT / "figures"
# Non–main-text figure mirrors (PCA supplement, NetMHC coding-fragment cohort mirrors, etc.).
FIGURES_SUPPLEMENTARY = FIGURES / "supplementary"
FIGURES_SUPPLEMENTARY_PCA = FIGURES_SUPPLEMENTARY / "pca"
FIGURES_SUPPLEMENTARY_NETMHC = FIGURES_SUPPLEMENTARY / "netmhc"
FIGURES_SUPPLEMENTARY_NETMHC_CODING_FRAGMENTS = FIGURES_SUPPLEMENTARY_NETMHC / "coding_fragments_random_sample"
# Fig 1 OpenTSNE supplement (dims 1–2 and 3–4) and other alternate embeddings live here.
FIGURES_SUPPLEMENTARY_EMBEDDING = FIGURES_SUPPLEMENTARY / "embedding"
# Fig 6 TTN-AS1 split panels (e.g. unique coverage) when not part of main-text bundle.
FIGURES_SUPPLEMENTARY_FIG6_TTN = FIGURES_SUPPLEMENTARY / "figure6_ttn_as1"
# Catalog Fig 2–3 side outputs (TCGA-matrix vs all-filtered peptide modes): data tables +
# mode-specific PNGs live under supplementary/; canonical copies (fig2b, fig3a–b, fig3c TCGA, fig3d)
# remain at ``figures/`` root when applicable.
FIGURES_SUPPLEMENTARY_TCGA_MATRIX = FIGURES_SUPPLEMENTARY / "tcga_matrix"
FIGURES_SUPPLEMENTARY_ALL_SMPROT_FILTERED = FIGURES_SUPPLEMENTARY / "all_smprot_filtered"
# Exploratory limma–z stratum log2FC histograms (supplement/plot_z_stratum_logfc_histograms.py).
FIGURES_SUPPLEMENTARY_Z_STRATUM_LOGFC = FIGURES_SUPPLEMENTARY / "z_stratum_logfc_histograms"
# Fig 5 cohort SB sensitivity + SB combination grid (merged TSV) + Fig 6 TTN NetMHC-only SB sweeps.
# Populated by ``generate_netmhc_fig5_fig6_supplement.py`` (see ``docs/figure_generation_overview.md``).
FIGURES_SUPPLEMENTARY_NETMHC_FIG5_FIG6 = FIGURES_SUPPLEMENTARY / "netmhc_fig5_fig6_supplement"
# Canonical figures: PDF + TIFF exports (see export_publication_figures.py).
FIGURES_PUBLICATION = FIGURES / "publication"
NETMHC_DATA = DATA / "netmhc"
NETMHC_FIGURES = NETMHC_DATA / "figures"
# European 27 class-I panel population frequencies (reference table for merged Fig 5A).
NETMHC_HLA27_ALLELE_FREQ_CSV = NETMHC_DATA / "hla_european27_allele_frequencies.csv"
SCRIPTS_DIR = REPO_ROOT / "scripts"
MANUSCRIPT_DIR = REPO_ROOT / "manuscript"
SUPPLEMENT_DIR = REPO_ROOT / "supplement"
PIPELINE_DIR = REPO_ROOT / "pipeline"
