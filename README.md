# CENSO - Commandline ENergetic SOrting of Conformer-Rotamer Ensembles
![censo_logo_300dpi](https://github.com/user-attachments/assets/0e6bac6a-2637-4207-8eca-3122ab90112a)
CENSO is a Python package meant to automate refinement of Conformer-Rotamer ensembles on DFT level, as well as calculation of some ensemble properties, e.g. NMR parameters.
It can be used from the command line as well as to define custom pipelines.

## Update: CENSO 2.2
The new update 2.2 improves upon version 2.0 by significantly upgrading the program's architecture, which is especially important for usage as a Python package.
Configuration is handled using Pydantic V2.
Some keywords were modified, please refer to the example rcfile for reference.
Also, new auxiliary command line scripts were added (`nmrplot`, `uvvisplot`, `c2anmr`).

# Installation
Can be installed using `pip` by running

    pip install .

If you want to install and run CENSO without `pip` you can add the `CENSO/src` directory to your `$PYTHONPATH` and add `CENSO` to your `$PATH`.

To install further entry points (new command line auxiliary scripts) and their dependencies run
```
    pip install .[scripts]
```

# Usage
An entry point `censo` will be installed by pip.
For information about command line options use the `-h` option.

The helper script to run CENSO without pip installation located in `CENSO` is also called via `censo`.

> Previous versions required the ``-i`` and ``--maxcores`` keywords. Since version 2.1.3 these are no longer required 
> and instead assume default values, ``crest_conformers.xyz`` and the total number of CPU cores on the machine 
> running CENSO, respectively.

CENSO can also be used as a package within Python to define custom pipelines. 
A basic setup for a CENSO run in a Python file could look like this:
```python
from censo.ensembledata import EnsembleData
from censo.configuration import configure
from censo.ensembleopt import prescreening, screening, optimization
from censo.properties import nmr
from censo.config import GeneralConfig
from censo.config.parallel_config import ParallelConfig
from censo.parallel import Config

# CENSO will put all files in the current working directory (os.getcwd())
# To output to different dirs you need to start subprocesses in different working directories using e.g. subprocess.Popen
input_path = "rel/path/to/your/inputfile" # path relative to the working directory
ensemble = EnsembleData(input_file=input_path) 
# the above can be used if you molecule is neutral and closed shell, otherwise
# it is necessary to proceed with e.g.
# ensemble = EnsembleData()
# ensemble.read_input(input_path, charge=-1, unpaired=1)

# If the user wants to use a specific rcfile:
config = configure(rcpath="/path/to/rcfile")

# Configure specific settings for parallelization
parallel_config = ParallelConfig(ncores=os.cpu_count(), omp=1)
# Note that os.cpu_count() can return None which is an invalid value for ncores

# After changing a setting you should revalidate
# Alternatively you could also enable assignment validation using:
# GeneralConfig.model_config["validate_assignment"] = True
# This can however temporarily lead to blocked combinations (like trying to use COSMORS with ORCA)
config.general.solvent = "dmso"
config = config.model_validate(config)

# Setup and run all the parts that the user wants to run
# Running the parts in order here, while it is also possible to use a custom order or run some parts multiple times
# Running a part will return an instance of the respective type
# References to the resulting part instances will be appended to a list in the EnsembleData object (ensemble.results)
# Note though, that currently this will lead to results being overwritten in your working directory
# (you could circumvent this by moving/renaming the folders)
timings = [part(ensemble, config, parallel_config=parallel_config) for part in [prescreening, screening, optimization, nmr]]
# parallel_config does not need to be provided

# You can also find all the results in the <part>.json output files
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
