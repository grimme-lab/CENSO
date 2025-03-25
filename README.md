# NEW: CENSO 2.0
This is the updated version of the former CENSO 1.3 program. New features include the possibility to use CENSO as a package from within Python, template files, json outputs, and more! For more information about the use and the capabilities of CENSO 2.0 visit the documentation [here](https://xtb-docs.readthedocs.io/en/latest/CENSO_docs/censo.html).

## Update: CENSO 2.2
The new update 2.2 improves upon version 2.0 by significantly upgrading the program's architecture, which is especially important for usage as a Python package.

# Installation
Can be installed using `pip` by running

    pip install .

If you want to install and run CENSO without `pip` you can add the `CENSO/src` directory to your `$PYTHONPATH` and add `CENSO/bin` to your `$PATH`.

# Usage
An entry point `censo` will be installed by pip. Furthermore the scripts under the `bin` directory will be added to your path.

For information about command line options use the `-h` option.

The helper script is also called via `censo`.

Previous versions required the ``-i`` and ``--maxcores`` keywords. Since version 2.1.3 these are no longer required 
and instead assume default values, ``crest_conformers.xyz`` and the total number of CPU cores on the machine 
running CENSO, respectively.

CENSO can also be used as a package within Python. A basic setup for a CENSO run in a Python file could look like this:
```python
from censo.ensembledata import EnsembleData
from censo.configuration import configure
from censo.ensembleopt import prescreening, screening, optimization
from censo.properties import nmr
from censo.params import NCORES, OMP
from censo.config import GeneralConfig

# CENSO will put all files in the current working directory (os.getcwd())
input_path = "rel/path/to/your/inputfile" # path relative to the working directory
ensemble = EnsembleData(input_file=input_path) 
# the above can be used if you molecule is neutral and closed shell, otherwise
# it is necessary to proceed with e.g.
# ensemble = EnsembleData()
# ensemble.read_input(input_path, charge=-1, unpaired=1)

# If the user wants to use a specific rcfile:
configure("/abs/path/to/rcfile")

# Get the number of available cpu cores on this machine
# This is also the default value that CENSO uses
# This number can also be set to any other integer value and automatically checked for validity
NCORES = os.cpu_count()

# Another possibly important setting is OMP, which will get used if you disabled the automatic 
# load balancing in the settings
OMP = 4

# The user can also choose to change specific settings of the parts
# Please take note of the following:

# It is also possible to use a dict to set multiple values in one step
settings = {
    "threshold": 3.5,
    "func": "pbeh-3c",
    "implicit": True,
}

# To temporarily disable assignment validation for a specific config:
GeneralConfig.model_config["validate_assignment"] = False
config.general.solvent = "dmso"

# Setup and run all the parts that the user wants to run
# Running the parts in order here, while it is also possible to use a custom order or run some parts multiple times
# Running a part will return an instance of the respective type
# References to the resulting part instances will be appended to a list in the EnsembleData object (ensemble.results)
# Note though, that currently this will lead to results being overwritten in your working directory
# (you could circumvent this by moving/renaming the folders)
timings = [part(ensemble, config) for part in [prescreening, screening, optimization, nmr]]

# You access the results using the ensemble object
# You can also find all the results the <part>.json output files
```

# License

``CENSO`` is free software: you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

``CENSO`` is distributed in the hope that it will be useful,
but without any warranty; without even the implied warranty of
merchantability or fitness for a particular purpose. See the
GNU Lesser General Public License for more details.
