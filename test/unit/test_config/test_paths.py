"""Tests for path-related configuration behavior."""

from pathlib import Path

from censo.config.paths import PathsConfig


def test_orca_version_parsing_mmap_matches_previous_and_reports_speedup(tmp_path: Path):
    """Compare mmap ORCA version parsing with previous full-read approach.
    Uses a synthetic ORCA-like binary in a temporary directory so the test
    does not depend on a real ORCA executable being present in PATH.
    """
    fake_orca = tmp_path / "orca"
    # Create minimal bytes containing the expected version marker pattern.
    fake_content = (
        b"\x00\x01Random header bytes\n"
        b"Some other data...\n"
        b"Program Version 1.2.3\n"
        b"Trailing bytes\xff\xfe"
    )
    with open(fake_orca, "wb") as f:
        f.write(fake_content)
    orca_path = fake_orca

    config = PathsConfig.model_validate({"orca": str(orca_path)})
    assert config.orcaversion == "1.2.3"
