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
NETMHC_DATA = DATA / "netmhc"
NETMHC_FIGURES = NETMHC_DATA / "figures"
SCRIPTS_DIR = REPO_ROOT / "scripts"
MANUSCRIPT_DIR = REPO_ROOT / "manuscript"
SUPPLEMENT_DIR = REPO_ROOT / "supplement"
PIPELINE_DIR = REPO_ROOT / "pipeline"
