# Installation
Can be installed using `pip` by running

    pip install .

If you want to install and run `CENSO` without `pip` you can add the `CENSO/src` directory to your `$PYTHONPATH` and add `CENSO/bin` to your `$PATH`.
> NOTE: you may want to change the `maxcores` setting in your `.censo2rc` to something more than 4.
> It indicates the maximum number of cores `CENSO` is allowed to use.

# Usage
Basic usage: 

    python3 -m censo -inp [path to ensemble input] ...

For information about command line options use the `-h` option.

If you want to run it via helper script after adding it to your `$PATH`:

    censo -inp [path to ensemble input] ...

CENSO can also be used as a package. A basic setup for a CENSO run in a Python file would look like this (WIP):
```python
from censo.ensembledata import EnsembleData
from censo.configuration import configure
from censo.ensembleopt import Prescreening, Screening, Optimization
from censo.properties import NMR

workdir = "/absolute/path/to/your/workdir" # CENSO will put all files in this directory
input_path = "rel/path/to/your/inputfile" # path relative to the working directory
ensemble = EnsembleData(workdir)
ensemble.read_input(input_path, charge=0, unpaired=0)

# If the user wants to use a specific rcfile:
configure("/abs/path/to/rcfile")

# Setup all the parts that the user wants to run
parts = [
    part(ensemble) for part in [Prescreening, Screening, Optimization, NMR]
]

# Run all the parts and collect their runtimes
part_timings = []
for part in parts:
    part_timings.append(part.run())

# If no Exceptions were raised, all the output can now be found in 'workdir'
# Data is given in a formatted plain text format (*.out) and and json format
# The files used in the computations for each conformer can be found in the folders 
# generated by each part, respectively (e.g. 'prescreening/CONF2/...')
```
