"""Shared subprocess helpers for repo orchestrators (echo command + propagate exit code)."""
from __future__ import annotations

import subprocess
from pathlib import Path
from typing import Sequence, Union

PathLike = Union[Path, str]


def call_echo(cmd: Sequence[str], *, cwd: PathLike) -> int:
    print("+", " ".join(cmd), flush=True)
    return subprocess.call(list(cmd), cwd=str(cwd))


def run_echo(cmd: Sequence[str], *, cwd: PathLike) -> subprocess.CompletedProcess:
    print("+", " ".join(cmd), flush=True)
    return subprocess.run(list(cmd), cwd=str(cwd))
