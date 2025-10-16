#!/usr/bin/env python3
import argparse
from pathlib import Path
import sys
import json
import math
from typing import Any, cast
import shutil

# A type alias for the conformer data structure
ConformerData = dict[str, int | float | str]


def parse() -> argparse.Namespace:
    """
    Parses command-line arguments for ANMR setup.

    :return: Parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(
        "c2anmr",
        description="Generates files and sets up directory for calculation of NMR spectra using ANMR from a CENSO run.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Positional arguments
    parser.add_argument(
        "nmr_file",
        nargs="?",
        default="4_NMR.json",
        help="Path to the JSON file containing conformer energies and NMR parameters.",
    )

    # General settings
    parser.add_argument(
        "--mf",
        type=int,
        default=300,
        help="Carrier wave frequency of the NMR spectrometer in MHz.",
    )
    parser.add_argument("--lw", type=float, default=1.0, help="Linewidth for plotting.")
    parser.add_argument(
        "--j",
        choices=["on", "off"],
        default="on",
        help="Evaluate couplings ('on' or 'off').",
    )
    parser.add_argument(
        "--s",
        choices=["on", "off"],
        default="on",
        help="Evaluate shifts ('on' or 'off').",
    )
    parser.add_argument(
        "--T",
        type=float,
        default=298.15,
        help="Temperature in K for Boltzmann weighting.",
    )

    # Reference shieldings
    parser.add_argument(
        "--ref",
        nargs=4,
        action="append",
        metavar=("ATOM_NUM", "SHIELDING", "EXP_SHIFT", "ACTIVE"),
        help="Reference shielding data. Format: ATOM_NUM SHIELDING EXP_SHIFT ACTIVE. Can be specified multiple times for different elements."
        + " Reference values can also be specified as a file anmr.ref in the same format.",
    )

    args: argparse.Namespace = parser.parse_args()
    return args


def create_conformers_list(
    json_file_path: Path | str, temperature: float
) -> list[ConformerData] | None:
    """
    Reads conformer data from a JSON file, calculates Boltzmann weights,
    and returns a list of dictionaries for the anmr_enso file.

    :param json_file_path: Path to the JSON file containing conformer data.
    :param temperature: Temperature in K for Boltzmann weighting.
    :return: List of conformer data dictionaries or None if error.
    """
    K_HARTREE: float = 3.1668114e-6  # Boltzmann constant in Hartree/K

    try:
        with open(json_file_path) as f:
            json_data: dict[str, Any] = json.load(f)
    except FileNotFoundError:
        print(f"Error: The file {json_file_path} was not found.", file=sys.stderr)
        return None
    except json.JSONDecodeError:
        print(f"Error: Could not decode JSON from {json_file_path}.", file=sys.stderr)
        return None

    conformer_raw_data: dict[str, Any] | None = json_data.get("data")
    if not conformer_raw_data:
        print(f"Error: No 'data' key found in {json_file_path}.", file=sys.stderr)
        return None

    for d in conformer_raw_data.values():
        d["gtot"] = d["energy"] + d["gsolv"] + d["grrho"]

    try:
        min_gtot: float = min(d["gtot"] for d in conformer_raw_data.values())
    except (KeyError, TypeError):
        print(f"Error: Malformed 'gtot' data in {json_file_path}.", file=sys.stderr)
        return None

    kt: float = K_HARTREE * temperature
    partition_q: float = sum(
        math.exp(-(d["gtot"] - min_gtot) / kt) for d in conformer_raw_data.values()
    )

    if partition_q == 0:
        print(
            "Error: Partition function is zero, cannot calculate Boltzmann weights.",
            file=sys.stderr,
        )
        return None

    conformers_data: list[ConformerData] = []
    for name, values in conformer_raw_data.items():
        try:
            conf_index: int = int(name.replace("CONF", ""))
            boltzmann_weight: float = (
                math.exp(-(values["gtot"] - min_gtot) / kt) / partition_q
            )

            conformers_data.append(
                {
                    "ONOFF": 1,
                    "NMR": conf_index,
                    "CONF": conf_index,
                    "BW": boltzmann_weight,
                    "Energy": values["energy"],
                    "Gsolv": values["gsolv"],
                    "RRHO": values["grrho"],
                    "shieldings": values["shieldings"],
                    "couplings": values["couplings"],
                    "nat": int(values["nat"]),
                }
            )
        except (KeyError, ValueError) as e:
            print(
                f"Warning: Skipping malformed conformer entry '{name}'. Reason: {e}",
                file=sys.stderr,
            )
            continue

    conformers_data.sort(key=lambda x: x["CONF"])  # type: ignore

    print(
        f"Info: Processed {len(conformers_data)} conformers from {json_file_path}",
        file=sys.stderr,
    )
    return conformers_data


def write_nmrprop(
    directory: Path, nat: int, shieldings: list[Any], couplings: list[Any]
) -> None:
    """
    Writes the nmrprop.dat file for a single conformer.

    :param directory: Directory to write the file in.
    :param nat: Number of atoms.
    :param shieldings: List of shieldings.
    :param couplings: List of couplings.
    :return: None
    """
    lines: list[str] = []

    # Process shieldings
    if shieldings:
        shielding_map: dict[int, float] = {s[0]: s[1] for s in shieldings}
        for i, shielding in shielding_map.items():
            lines.append(f"{i + 1:4} {shielding_map[i]:.3f}\n")

        for i in range(nat - len(shielding_map)):
            lines.append("\n")

    # Process couplings
    if couplings:
        for c in couplings:
            i, j = c[0]
            coupling_val = c[1]
            lines.append(f"{i + 1:4} {j + 1:4} {coupling_val:.3f}\n")

    # Write the file
    file_path = directory / "nmrprop.dat"
    try:
        with open(file_path, "w") as f:
            f.writelines(lines)
    except OSError as e:
        print(f"Error writing to file {file_path}: {e}", file=sys.stderr)


def write_anmrrc(args: argparse.Namespace, directory: Path | str) -> None:
    """
    Writes the .anmrrc file to the specified directory.

    :param args: Parsed command-line arguments.
    :param directory: Directory to write the file in.
    :return: None
    """
    lines: list[str] = []
    lines.append("7 8 XH acid atoms")
    lines.append(
        f"ENSO qm= TM mf= {args.mf} lw= {args.lw}  J= {args.j} S= {args.s} T= {args.T}"
    )
    lines.append("TMS[chcl3] pbe0[COSMO]/def2-TZVP//pbeh-3c[DCOSMO-RS]/def2-mSVP")
    for ref_data in args.ref:
        lines.append(
            f"{ref_data[0]:<2} {ref_data[1]:<9} {ref_data[2]:<6} {ref_data[3]}"
        )

    file_path: Path = Path(directory) / ".anmrrc"
    try:
        file_path.write_text("\n".join(lines) + "\n")
        print(f"Info: Written ANMR config to {file_path}", file=sys.stderr)
    except OSError as e:
        print(f"Error: Could not write to {file_path}. {e}", file=sys.stderr)
        sys.exit(1)


def write_anmr_enso(
    conformers_data: list[ConformerData], output_filename: Path | str
) -> None:
    """
    Writes the anmr_enso file from a list of conformer data.

    :param conformers_data: List of conformer data dictionaries.
    :param output_filename: Path to the output file.
    :return: None
    """
    header: str = (
        f"{'ONOFF':<6}{'NMR':<5}{'CONF':<5}{'BW':>10}{'Energy':>12}{'Gsolv':>10}{'RRHO':>10}\n"
    )
    try:
        with open(output_filename, "w") as f:
            f.write(header)
            for conf in conformers_data:
                line: str = (
                    f"{conf['ONOFF']:<6d}"
                    f"{conf['NMR']:<5d}"
                    f"{conf['CONF']:<5d}"
                    f"{conf['BW']:>10.5f}"
                    f"{conf['Energy']:>12.5f}"
                    f"{conf['Gsolv']:>10.5f}"
                    f"{conf['RRHO']:>10.5f}\n"
                )
                f.write(line)
        print(f"Info: Written anmr_enso file to {output_filename}", file=sys.stderr)
    except OSError as e:
        print(f"Error writing to file {output_filename}: {e}", file=sys.stderr)
    except KeyError as e:
        print(f"Error: Missing key {e} in conformers_data.", file=sys.stderr)


def load_references_from_config(
    config_path: Path | str,
) -> list[list[str]] | None:
    """
    Loads reference shieldings from a local config file.

    :param config_path: Path to the config file.
    :return: List of reference shieldings or None if error.
    """
    refs: list[list[str]] = []
    try:
        with open(config_path) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts: list[str] = line.split()
                if len(parts) == 4:
                    refs.append(parts)
                else:
                    print(
                        f"Warning: Skipping malformed line in {config_path}: {line}",
                        file=sys.stderr,
                    )
        print(f"Info: Loaded reference shieldings from {config_path}", file=sys.stderr)
        return refs
    except OSError:
        return None


def main() -> None:
    """
    Main execution function.

    :return: None
    """
    args: argparse.Namespace = parse()

    dest_dir: Path = Path("anmr")

    if args.ref is None:
        config_file: Path = Path.cwd() / "anmr.ref"
        refs_from_file: list[list[str]] | None = load_references_from_config(
            config_file
        )
        if refs_from_file is not None:
            args.ref = refs_from_file
        else:
            print(
                "Info: No --ref argument or anmr.ref file found. Using default shieldings.",
                file=sys.stderr,
            )
            args.ref = [
                ["1", "31.786", "0.0", "1"],
                ["6", "189.674", "0.0", "0"],
                ["9", "182.57", "0.0", "0"],
                ["15", "291.9", "0.0", "0"],
            ]

    dest_dir.mkdir(parents=True, exist_ok=True)

    # Process JSON and write anmr_enso file
    conformers: list[ConformerData] | None = create_conformers_list(
        args.nmr_file, args.T
    )
    if conformers:
        write_anmr_enso(conformers, output_filename=dest_dir / "anmr_enso")
    else:
        print(f"No conformers could be read from {args.nmr_file}.")
        dest_dir.rmdir()
        sys.exit(2)

    # Write the .anmrrc file
    write_anmrrc(args, dest_dir)

    # Copy anmr_nucinfo and anmr_rotamer
    for file in ("anmr_rotamer", "anmr_nucinfo"):
        file_path = Path(file)
        if not file_path.is_file():
            print(f"Warning: {file} not found. ANMR cannot run without it!")
        else:
            shutil.copy(file_path, dest_dir)

    # Loop through conformers
    for conf in conformers:
        conf_dir = dest_dir / ("CONF" + str(conf["CONF"]))

        new_nmr_dir: Path = conf_dir / "NMR"
        new_nmr_dir.mkdir(parents=True, exist_ok=True)

        write_nmrprop(
            new_nmr_dir,
            cast(int, conf["nat"]),
            cast(list[Any], conf["shieldings"]),
            cast(list[Any], conf["couplings"]),
        )

    print("\nFinished setting up anmr directory. Please check .anmrrc for correctness.")
    print("Run anmr -plain")


if __name__ == "__main__":
    main()
