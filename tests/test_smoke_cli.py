"""Smoke tests: CLI --help and import hygiene for manuscript and supplement scripts."""
from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent

_UTF8_ENV = {**os.environ, "PYTHONUTF8": "1", "PYTHONIOENCODING": "utf-8"}


def _py_help(script: Path) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, str(script), "--help"],
        cwd=str(REPO),
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        timeout=120,
        check=False,
        env=_UTF8_ENV,
    )


@pytest.mark.parametrize(
    "rel",
    sorted(p.relative_to(REPO).as_posix() for p in (REPO / "manuscript").glob("*.py")),
)
def test_manuscript_help(rel: str) -> None:
    r = _py_help(REPO / rel)
    assert r.returncode == 0, f"{rel} --help failed: {r.stderr[:2000]}"


@pytest.mark.parametrize(
    "rel",
    sorted(p.relative_to(REPO).as_posix() for p in (REPO / "supplement").glob("*.py")),
)
def test_supplement_help(rel: str) -> None:
    r = _py_help(REPO / rel)
    assert r.returncode == 0, f"{rel} --help failed: {r.stderr[:2000]}"


def test_merge_script_help() -> None:
    r = _py_help(REPO / "scripts" / "merge_netmhcpan_xls_with_iedb.py")
    assert r.returncode == 0


def test_netmhc_sb_core_help() -> None:
    r = _py_help(REPO / "scripts" / "netmhc_sb_core.py")
    assert r.returncode == 0
