#!/usr/bin/env python3
import shutil
from pathlib import Path
import sys

def main():
    # Define source and destination directories
    src_dir = Path("4_NMR")
    dest_dir = Path("anmr")

    # Create destination directory
    dest_dir.mkdir(parents=True, exist_ok=True)

    # Copy all files matching anmr_* into dest_dir
    for file in Path.cwd().glob("anmr_*"):
        if file.is_file():
            shutil.copy(file, dest_dir / file.name)

    # Loop through CONF* subdirectories in src_dir
    for conf_dir in src_dir.glob("CONF*"):
        if not conf_dir.is_dir():
            continue

        # Build target NMR directory under dest_dir/CONF#/NMR
        new_nmr_dir = dest_dir / conf_dir.name / "NMR"
        new_nmr_dir.mkdir(parents=True, exist_ok=True)

        # Copy the two specific files if they exist
        for fname in ("nmrprop.dat", "coord"):
            src_file = conf_dir / fname
            if src_file.exists():
                shutil.copy(src_file, new_nmr_dir / fname)
            else:
                print(f"Warning: {src_file} not found", file=sys.stderr)

if __name__ == "__main__":
    main()

