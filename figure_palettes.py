"""
Color-blind-friendly plotting defaults (Okabe–Ito + CVD-safe colormaps).

References: Wong, Nature Methods 8, 441 (2011); see also Paul Tol’s notes.
Import as ``import figure_palettes as pal`` from the repository root on ``sys.path``.
"""
from __future__ import annotations

from functools import lru_cache

from matplotlib.colors import LinearSegmentedColormap

# --- Okabe–Ito (color-blind safe) ---
OI_ORANGE = "#E69F00"
OI_SKY_BLUE = "#56B4E9"
OI_BLUISH_GREEN = "#009E73"
OI_YELLOW = "#F0E442"
OI_BLUE = "#0072B2"
OI_VERMILLION = "#D55E00"
OI_REDDISH_PURPLE = "#CC79A7"
OI_BLACK = "#000000"

# --- Manuscript roles (two-category + accents) ---
SIG_LNC = OI_VERMILLION
CODING_CONTROL = OI_BLUE
BAR_SINGLE_SERIES = OI_BLUE
# Fig 5D/5E: same allele → same edge color across lncRNA vs coding bars (rank-shift callouts)
FIG5DE_TRACK_ALLELE_EDGES: dict[str, str] = {
    "HLA-A*31:01": OI_ORANGE,
    "HLA-B*15:01": OI_BLUISH_GREEN,
    "HLA-B*08:01": OI_REDDISH_PURPLE,
    "HLA-A*30:01": OI_SKY_BLUE,
}
SCATTER_DEFAULT = OI_SKY_BLUE
LINE_PRIMARY = OI_BLUE
LINE_SECONDARY = OI_VERMILLION
HIST_FILL = OI_SKY_BLUE

VOLCANO_ENRICHED = OI_VERMILLION
VOLCANO_DEPLETED = OI_BLUE
VOLCANO_SMALL_EFFECT = OI_YELLOW
VOLCANO_NS = "#999999"

AA_FREQ_TCGA = OI_BLUE
AA_FREQ_PROTEOME = OI_ORANGE

LENGTH_HIST_SINGLE = OI_BLUE
LENGTH_HIST_ALT = OI_REDDISH_PURPLE

# Fig 4A (dark theme): regions still distinguishable without red–green reliance
F4A_SHADE_TIS_STRONG = OI_BLUISH_GREEN
F4A_SHADE_RIBO_STRONG = OI_SKY_BLUE
F4A_SHADE_BOTH_STRONG = OI_REDDISH_PURPLE
F4A_POINTS_BULK = OI_ORANGE
F4A_POINT_EXTREME = OI_VERMILLION
F4A_GRID = "#cccccc"
F4A_FRAME = "#ffffff"

# Gene callout frames (dark fill + bright edge)
F4A_FRAME_STYLES: dict[str, dict[str, str]] = {
    "PTPRG-AS1": {"edgecolor": OI_VERMILLION, "facecolor": "#1a0a06"},
    "LINC00326": {"edgecolor": OI_REDDISH_PURPLE, "facecolor": "#120818"},
    "LINC00958": {"edgecolor": OI_SKY_BLUE, "facecolor": "#050a18"},
}

# Fig 4A light theme (white figure/axes; Okabe–Ito–safe contrast)
F4A_LIGHT_AX = "#ffffff"
F4A_LIGHT_TEXT = "#1a1a1a"
F4A_LIGHT_TICK = "#333333"
F4A_LIGHT_GRID_MAJOR = "#9a9a9a"
F4A_LIGHT_GRID_MINOR = "#c8c8c8"
F4A_LIGHT_FRAME = "#222222"
F4A_LIGHT_POINTS_BULK = OI_BLUE
F4A_LIGHT_SHADE_ALPHA = 0.35
F4A_LIGHT_FRAME_STYLES: dict[str, dict[str, str]] = {
    "PTPRG-AS1": {"edgecolor": OI_VERMILLION, "facecolor": "#fff5f0"},
    "LINC00326": {"edgecolor": OI_REDDISH_PURPLE, "facecolor": "#f8f0ff"},
    "LINC00958": {"edgecolor": OI_SKY_BLUE, "facecolor": "#f0f8ff"},
}

# Extra framed callouts for multi-MP genes (Fig 4A); cycled after the three highlights above.
F4A_EXTRA_FRAME_STYLES_LIGHT: tuple[dict[str, str], ...] = (
    {"edgecolor": OI_ORANGE, "facecolor": "#fff8f0"},
    {"edgecolor": OI_BLUISH_GREEN, "facecolor": "#f0faf6"},
    {"edgecolor": OI_BLUE, "facecolor": "#f0f4ff"},
    {"edgecolor": OI_BLACK, "facecolor": "#f4f4f4"},
    {"edgecolor": OI_YELLOW, "facecolor": "#fffef0"},
)
F4A_EXTRA_FRAME_STYLES_DARK: tuple[dict[str, str], ...] = (
    {"edgecolor": OI_ORANGE, "facecolor": "#1a1206"},
    {"edgecolor": OI_BLUISH_GREEN, "facecolor": "#061410"},
    {"edgecolor": OI_BLUE, "facecolor": "#060a14"},
    {"edgecolor": "#aaaaaa", "facecolor": "#141414"},
    {"edgecolor": OI_YELLOW, "facecolor": "#141206"},
)


@lru_cache(maxsize=1)
def diverging_log2fc_cmap() -> LinearSegmentedColormap:
    """Blue ↔ white ↔ orange (no red–green)."""
    return LinearSegmentedColormap.from_list(
        "cvd_blue_orange",
        [OI_BLUE, "#f6f6f6", OI_ORANGE],
        N=256,
    )


def sequential_heatmap() -> str:
    """Matplotlib registered name; perceptually uniform, CVD-friendly."""
    return "cividis"


def sequential_count_heatmap() -> str:
    """For non-negative intensity grids (e.g. SB combination heatmap)."""
    return "cividis"
