# Installation
Can be installed using `pip` by running

    pip install .

If you want to run CENSO without calling Python directly you can add `CENSO` to your `$PYTHONPATH` and add `CENSO/bin` to your `$PATH`.

# Usage
Basic usage: 

    python3 -m censo -inp [path to ensemble input] ...

For information about command line options use the `-h` option.

If you want to run it via helper script after adding it to your `$PATH`:

    censo -inp [path to ensemble input] ...

CENSO can also be used as a package. A basic setup for a CENSO run in a Python file would look like this (WIP):
```python
    test_dir = "/absolute/path/to/your/workdir"
    input_path = "rel/path/to/your/inputfile" # path relative to your workdir
    core = CensoCore(test_dir)
    core.read_input(input_path, charge=0, unpaired=0)
    
    # TODO TODO TODO
```
