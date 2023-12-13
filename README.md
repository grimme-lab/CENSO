# Installation
Can be installed using pip by running
    pip install .
If you want to run CENSO without calling Python directly you can add CENSO to your $PYTHONPATH and add CENSO/bin to your path.

# Usage
Basic usage: 
    python3 -m censo -inp [path to ensemble input] -chrg [charge] -u [number of unpaired electrons] -solvent [solvent name]

Or alternatively if you want to run it via helper script:
    censo -inp [path to ensemble input] -chrg [charge] -u [number of unpaired electrons] -solvent [solvent name]
