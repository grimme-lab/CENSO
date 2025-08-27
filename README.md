# CENSO - Commandline and Python Package for Conformer-Rotamer Ensembles

![CENSO Logo](https://github.com/user-attachments/assets/0e6bac6a-2637-4207-8eca-3122ab90112a)

**CENSO** is a Python package designed to automate the refinement of **Conformer-Rotamer Ensembles** at the **DFT level**. It supports both **command-line usage** and integration as a **fully customizable Python package**. CENSO enables calculations of ensemble properties like **NMR parameters**, **optical rotation**, and **UV-Vis spectra**.

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

---

## **Update: CENSO 3.0**
Version **3.0** introduces significant architectural improvements, particularly enhancing its usability as a **Python package**. Overall usage is also improved, e.g.:
- automatic system path exploration (finding program binaries, program versions, etc.),
- improved clarity for printout,
- improved logging,
- more configuration options for command-line interface.

Configuration is now managed using **Pydantic V2**, and several keywords have been updated. Refer to `example.censo2rc` for details.

New auxiliary command-line scripts (`nmrplot`, `uvvisplot`, `c2anmr`) have been added for improved workflow support.

---

## **Installation**

### **Using `pip`**
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

---

### **Python Package Usage**
CENSO can be integrated into Python scripts for custom workflows. Below is a basic setup example:

```python
from censo.ensembledata import EnsembleData
from censo.configuration import configure
from censo.ensembleopt import prescreening, screening, optimization
from censo.properties import nmr
from censo.config import GeneralConfig
from censo.config.parallel_config import ParallelConfig
from censo.parallel import Config

# CENSO outputs files in the current working directory (os.getcwd())
input_path = "rel/path/to/your/inputfile"  # Relative to working directory
ensemble = EnsembleData(input_file=input_path)

# For charged/open-shell systems:
# ensemble = EnsembleData()
# ensemble.read_input(input_path, charge=-1, unpaired=1)

# Load a custom rcfile (optional)
config = configure(rcpath="/path/to/rcfile")

# Configure parallelization
parallel_config = ParallelConfig(ncores=os.cpu_count(), omp=1)

# Ensure valid configuration
config.general.solvent = "dmso"
config = config.model_validate(config)

# Execute workflow steps
timings = [
    part(ensemble, config, parallel_config=parallel_config)
    for part in [prescreening, screening, optimization, nmr]
]
```

> **Note**: Results are stored in `<part>.json` files. For multiple runs, rename or move output folders to avoid overwriting.

---

## **License**
CENSO is **free software** under the **GNU Lesser General Public License (LGPLv3)**. You may redistribute and modify it under the license terms. **No warranty** is provided. See the [LGPLv3](https://www.gnu.org/licenses/lgpl-3.0.html) for details.
