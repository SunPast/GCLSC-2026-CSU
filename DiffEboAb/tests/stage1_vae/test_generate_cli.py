"""CLI contract tests for scripts/stage1_vae/generate.py."""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
GEN = ROOT / "scripts" / "stage1_vae" / "generate.py"


def test_missing_out_exits_nonzero():
    proc = subprocess.run(
        [sys.executable, str(GEN), "--n", "1"],
        cwd=ROOT,
        capture_output=True,
        text=True,
    )
    assert proc.returncode != 0
    err = (proc.stderr + proc.stdout).lower()
    assert "-o" in err or "--out" in err or "required" in err


def test_help_lists_out():
    proc = subprocess.run(
        [sys.executable, str(GEN), "-h"],
        cwd=ROOT,
        capture_output=True,
        text=True,
    )
    assert proc.returncode == 0
    assert "--out" in proc.stdout
    assert "-o" in proc.stdout
