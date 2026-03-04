"""Tests for path-related configuration behavior."""

from pathlib import Path
import mmap
import re
import shutil
import time

import pytest


def _parse_orca_version_with_mmap(orca_path: Path) -> str:
    version_pattern = rb"Program Version (\d+\.\d+\.\d+)"
    with open(orca_path, "rb") as f:
        with mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ) as mm:
            match = re.search(version_pattern, mm)
    if match is None:
        raise ValueError(f"Could not parse ORCA version from {orca_path}")
    return match.group(1).decode("utf-8")


def _parse_orca_version_with_read(orca_path: Path) -> str:
    version_pattern = rb"Program Version (\d+\.\d+\.\d+)"
    with open(orca_path, "rb") as f:
        content = f.read()
    match = re.search(version_pattern, content)
    if match is None:
        raise ValueError(f"Could not parse ORCA version from {orca_path}")
    return match.group(1).decode("utf-8")


def test_orca_version_parsing_mmap_matches_previous_and_reports_speedup():
    """Compare mmap ORCA version parsing with previous full-read approach."""
    orca = shutil.which("orca")
    if orca is None:
        pytest.skip("ORCA is not present in PATH.")

    orca_path = Path(orca)

    version_from_mmap = _parse_orca_version_with_mmap(orca_path)
    version_from_read = _parse_orca_version_with_read(orca_path)
    assert version_from_mmap == version_from_read

    repeats = 25

    start = time.perf_counter()
    for _ in range(repeats):
        _parse_orca_version_with_read(orca_path)
    read_time = time.perf_counter() - start

    start = time.perf_counter()
    for _ in range(repeats):
        _parse_orca_version_with_mmap(orca_path)
    mmap_time = time.perf_counter() - start

    speedup = read_time / mmap_time if mmap_time > 0 else float("inf")
    print(
        "ORCA version parsing benchmark "
        f"(path={orca_path}, repeats={repeats}): "
        f"read={read_time:.6f}s, mmap={mmap_time:.6f}s, speedup={speedup:.2f}x"
    )
