# CENSO - Commandline and Python Package for Conformer-Rotamer Ensembles

![CENSO Logo](https://github.com/user-attachments/assets/0e6bac6a-2637-4207-8eca-3122ab90112a)

**CENSO** is a Python package designed to automate the refinement of **Conformer-Rotamer Ensembles** at the **DFT level**. It supports both **command-line usage** and usage via a fully implemented **Python API**. CENSO enables calculations of ensemble properties like **NMR parameters**, **optical rotation**, and **UV-Vis spectra**.

---

## **Key Features**
- **Prescreening, Screening, Optimization, and Refinement** of conformer-rotamer ensembles.
  - **Geometry Optimization** using:
    - ORCA (native optimizer)
    - xtb (ANCOPT)
  - **Thermal contributions** using xtb single-point Hessians (`--bhess`)
- **Solvent Models**:
  - ORCA: **CPCM**, **SMD**
  - Turbomole: **COSMO**, **DCOSMORS**, **COSMORS**
  - xtb: **ALPB**, **GBSA**
- **Property Calculations**:
  - **NMR**: ORCA, Turbomole
  - **Optical Rotation**: Turbomole
  - **UV/Vis**: ORCA
- **Agnostic to global optimizer**, input files just need to be in xyz-format
- CLI and Python API

---

## **Update: CENSO 3.0**
Version **3.0** introduces significant architectural improvements, particularly enhancing its usability as a **Python package**. Overall usage is also improved, e.g.:
- automatic system path exploration (finding program binaries, program versions, etc.),
- automatic system recognition,
- improved clarity for printout,
- improved logging,
- more configuration options for command-line interface.

Configuration is now managed using **Pydantic V2**, and several keywords have been updated. Refer to `example.censo2rc` for details.

Resources and tasks are now managed using Dask. This is in preparation of a multi-node implementation.

New auxiliary command-line scripts (`nmrplot`, `uvvisplot`, `c2anmr`) have been added for improved workflow support.
`nmrplot` and `uvvisplot` are utilities for immediate visualization of NMR and UV/Vis data, respectively.
`c2anmr` is a tool to ensure backwards compatibility with ANMR, due to the new directory structure since CENSO 2.0.

The documentation has also been updated to reflect the new version and includes more guidance on usage.
It now also includes an example for the calculation of the NMR spectrum of glycerol.

---

## **Installation**

### **Using `pip`**
From pypi:
```bash
pip install censo
```
or after cloning the repository:
```bash
pip install .
```
To install additional command-line scripts and dependencies:
```bash
pip install .[scripts]
```

### **Manual Installation**
1. Add `CENSO/src` to your `$PYTHONPATH`.
2. Add `CENSO` to your `$PATH`.

---

## **Usage**

### **Command-Line Interface (CLI)**
After installation, use the `censo` command. For help, run:
```bash
censo -h
```

> **Note**: Opposed to previous versions, default values are now assumed for `-i` (`crest_conformers.xyz`) and `--maxcores` (total CPU cores).

In the CLI version, CENSO will try to automatically detect binary locations for necessary programs, as well as versions. If you want override, use the `[paths]` section in your `.censo2rc` file.
Furthermore, it will try to automatically detect the number of available CPU cores. To override these settings, as well as general settings, you can use the command-line arguments.

```
usage: censo [-h] [--prescreening] [--screening] [--optimization] [--refinement] [--nmr] [--rot] [--uvvis] [-i INP] [-n NCONF] [-c CHARGE] [-u UNPAIRED] [-v] [--cleanup] [--cleanup-all]
             [--new-config] [--inprc INPRCPATH] [--constraints CONSTRAINTS] [--maxcores MAXCORES] [--omp-min OMPMIN] [--loglevel {DEBUG,INFO,WARNING,ERROR,CRITICAL}] [--logpath LOGPATH]
             [--reload RELOAD [RELOAD ...]] [--keep-all] [--ignore-failed] [-T TEMPERATURE] [--trange start end step] [--bhess] [--consider-sym] [--sm-rrho SM_RRHO] [--evaluate-rrho]
             [--solvent SOLVENT] [--gas-phase] [--imagthr IMAGTHR] [--sthr STHR] [--scale SCALE]

Energetic sorting of Conformer Rotamer Ensembles (command line version).

options:
  -h, --help            show this help message and exit

RUN SETTINGS:
  --prescreening, -P    Enable the prescreening part.
  --screening, -S       Enable the screening part.
  --optimization, -O    Enable the optimization part.
  --refinement, -R      Enable the refinement part.
  --nmr                 Enable the NMR calculation part.
  --rot                 Enable the optical rotation calculation part.
  --uvvis               Enable the UV/Vis calculation part.
  -i, --input INP       Relative path to ensemble file, e.g. crest_conformers.xyz (default).
  -n, --nconf NCONF     The first 'nconf' conformers will be considered.
  -c, --charge CHARGE   Integer charge of the investigated molecule.
  -u, --unpaired UNPAIRED
                        Integer number of unpaired electrons of the investigated molecule.
  -v, --version         Print CENSO version and exit.
  --cleanup             Delete unneeded files from current working directory.
  --cleanup-all         Delete all CENSO files from previous runs from current working directory. Stronger than --cleanup !
  --new-config          Write new configuration file, which is placed into the current directory.
  --inprc INPRCPATH     Use to provide a path to the CENSO configuration file if you want to use a different one than the default (~/.censo2rc).
  --constraints CONSTRAINTS
                        Path to a file containing constraints in xtb format for use in geometry optimizations.
  --maxcores MAXCORES   Number of cores that should be used for CENSO on the machine. If this is not provided CENSO will use the maximum number available (also checks for slurm
                        environment variables).
  --omp-min OMPMIN      Minimum number of OpenMP threads per process, default is 4. This is mostly important if load balancing is enabled.
  --loglevel {DEBUG,INFO,WARNING,ERROR,CRITICAL}
                        Set the loglevel for all modules to a specified level.
  --logpath LOGPATH     Relative/absolute path to the logfile. If no path is provided, censo.log will be written in the output directory.
  --reload RELOAD [RELOAD ...]
                        Reload data from json output files. List all file names separated by spaces. Will be loaded in order (first to last). Note that all conformers from the current
                        ensemble need to be included in the output data keys.
  --keep-all            Do not cut down the ensemble, keep all conformers from start to end.
  --ignore-failed       Ignore failed calculations. If this is not set, failed calculations will raise RuntimeError.

GENERAL SETTINGS:
  -T, --temperature TEMPERATURE
                        Temperature in Kelvin for thermostatistical evaluation.
  --trange start end step
                        specify a temperature range [start, end, step] e.g.: 250.0 300.0 10.0 resulting in the range [250.0, 260.0, 270.0, 280.0, 290.0, 300.0].
  --bhess               Uses SPH and applies structure constraint to input/DFT geometry for mRRHO calcuation.
  --consider-sym        Consider symmetry in mRRHO calcuation (based on desy xtb threshold).
  --sm-rrho SM_RRHO     Solvation model used in xTB GmRRHO calculation. Applied if not in gas-phase. Options are 'gbsa' or 'alpb'.
  --evaluate-rrho       Evaluate mRRHO contribution.
  --solvent SOLVENT     Solvent to be used for Gsolv calculation.
  --gas-phase           Run calculation in gas-phase, overriding all solvation settings.
  --imagthr IMAGTHR     threshold for inverting imaginary frequencies for thermo in cm-1, e.g. -30.0.
  --sthr STHR           Rotor cut-off for thermo in cm-1, e.g. 50.0.
  --scale SCALE         Scaling factor for frequencies, e.g. 1.0.
```

---

### **Python API Usage**
CENSO can be integrated into Python scripts for custom workflows. Below is a basic setup example:

```python
from censo.ensembledata import EnsembleData
from censo.configuration import configure
from censo.ensembleopt import prescreening, screening, optimization
from censo.properties import nmr
from censo.config import GeneralConfig
from censo.parallel import get_cluster

# CENSO outputs files in the current working directory (os.getcwd())
input_path = "rel/path/to/your/inputfile"  # Relative to working directory
ensemble = EnsembleData(input_file=input_path)

# For charged/open-shell systems:
# ensemble = EnsembleData()
# ensemble.read_input(input_path, charge=-1, unpaired=1)

# Load a custom rcfile (optional)
config = configure(rcpath="/path/to/rcfile")

# Configure dask client
cluster = get_cluster()
client = cluster.get_client()

# Ensure valid configuration
config.general.solvent = "dmso"
config = config.model_validate(config)

# Execute workflow steps
results = [
    part(ensemble, config, client)
    for part in [prescreening, screening, optimization, nmr]
]
```

> **Note**: Results are also stored in `<part>.json` files. For multiple runs, rename or move output folders to avoid overwriting.

---

## Development Environment Setup

- **Python Version**: Requires Python >= 3.12. Set up a virtual environment using your preferred method:
  - Conda: `conda env create -f environment.yaml`
  - Pipenv, venv, or others.

- **Install the Package**: Install the package and its dependencies in development mode:
  ```bash
  pip install -e '.[dev]'
  ```

- **Pre-commit Hooks**: Install pre-commit hooks:
  ```bash
  pre-commit install
  ```

- **Testing**:
  - **Pytest Markers**: Custom markers are available for conditional test execution. Optional tests (e.g., integration tests) require third-party programs like **xtb**, **ORCA**, **Turbomole**, or **COSMOtherm**.
    ```python
    markers = [
        "optional: mark test as optional (e.g., integration tests)",
        "requires_xtb: skip test if xtb is not available",
        "requires_orca: skip test if ORCA is not available",
        "requires_turbomole: skip test if Turbomole is not available",
        "requires_cosmotherm: skip test if COSMOtherm is not available",
    ]
    ```
  - **Run Tests**:
    - Default (skip optional tests):
      ```bash
      pytest
      ```
    - Include optional tests:
      ```bash
      pytest --run-optional
      ```
  - **Run Tox**: Execute tests across multiple environments:
    ```bash
    tox
    ```

---

## **License**
CENSO is **free software** under the **GNU Lesser General Public License (LGPLv3)**. You may redistribute and modify it under the license terms. **No warranty** is provided. See the [LGPLv3](https://www.gnu.org/licenses/lgpl-3.0.html) for details.
