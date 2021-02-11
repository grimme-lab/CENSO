[![Latest Version](https://img.shields.io/github/v/release/grimme-lab/CENSO)](https://github.com/grimme-lab/CENSO/releases/latest)

CENSO - Commandline ENergetic SOrting of Conformer Rotamer Ensembles
====================================================================

This repository hosts the `CENSO` code for the refinement of Conformer Rotamer 
Ensembles (CRE) as obtained from `CREST`.



Installation
------------

There are several options possible. The easiest is to use the packaged censo programs
(by use of Pyinstaller) which can be found at the release section. The packaged 
censo is linked against GLIBC version 2.19 and will work for GLIBC version 2.19 and above.

Other options to use censo are shown below:


Download the git repository and run:

.. code::

    $ pip install --upgrade pip
    $ pip install --editable .
    $ censo arg1 arg2


Flexible Invocation
-------------------

1) Treating the censo directory as a package and as the main script::

    $ python3 -m censo arg1 arg2

2) Using the censo-runner.py wrapper::

    $ ./censo-runner.py arg1 arg2

3) After installation with pip::

    $ censo arg1 arg2



Getting started:
----------------

Create the remote configuration file .censorc where the user can adjust default
settings and provide paths to the external programs e.g. `xtb`, `crest`, `orca` ...

.. code::

    $ censo -newconfig
    $ cp censorc-new /home/$USER/.censorc
    # edit .censorc
    vi /home/$USER/.ensorc


Interactive Documentation can be accessed:

.. code::

    $ censo -tutorial


Explainations on the commandline arguments can be printed by:

.. code::

    $ censo --help

The molecule numbering from the input structure ensemble is kept throughout the 
entire program. There are several program parts which can be used to filter a structure 
ensemble:

0) Cheap prescreening (part0): Very fast DFT energies in order to improve upon the energy
    description of the SQM method used to generate the input structure ensemble.
    The (free) energies are evaluated on the input geometries (DFT unoptimized).

1) Prescreening (part1): Improved DFT energies and accurate solvation energies (if needed).
    The free energies are evaluated on the input geometries (DFT unoptimized).

2) Optimization (part2): efficient structure ensemble optimization and 
    free energy calculation on DFT optimized geometries.

3) Refinement (part3): Optional free energy refinement (on DFT optimized geometries).
    E.g. using hybrid DFA with large basis set.

4) NMR properties (part4): Optional calculation of shielding and coupling constants on 
    populated conformers.

5) Optical Rotation (part5): Optional calculation of optical rotatory dispersion 
    for the populated ensemble.


Usage:

.. code::

    # check if setting-combinations match:
    $ censo -inp structure_ensemble.xyz -part2 on -solvent h2o --checkinput
    # start the calculation:
    $ censo -inp structure_ensemble.xyz -part2 on -solvent h2o > censo.out 2> error.censo &

Requirements:
-------------

* newest xtb (currently: https://github.com/grimme-lab/xtb/releases/tag/v6.4.0 )
* newest cefine https://github.com/grimme-lab/cefine/releases
* ORCA > version 4.1


Further information (will be ordered later on):

* the file .censorc can be used in the current working directory and will be preferred to 
  the global configuration file in ~/.censorc
* a folder ~/.censo_assets/ will be created upon usage of censo
* ORCA has not been used extensively so please be careful, test calculations
  and report possible "bad" settings
* To be efficient COSMO-RS calculations are not performed with BP86 but whith the functionals
  for energy evaluation.




License
-------

LGPL3


Available solvation models:
---------------------------

Solvation models available for implicit effect on properties e.g. the 
geometry (SM). And "additive" solvation models which return a solvation contribution 
to free energy (Gibbs energy) of the choosen geometry (SMGSOLV).

.. csv-table:: 
    :header: "programs", "solvation models", "comment"
    
    "Turbomole","COSMO", "(SM)"
    "", "DCOSMO-RS","(SM)"
    "COSMO-RS","COSMO-RS","(SMGSOLV) (only solvent model for evaluation at different temperatures)"
    "ORCA", "CPCM", "(SM)"
    "","SMD","(SM)"
    "","SMD_GSOLV", "(SMGSOLV)"
    "xTB","GBSA_Gsolv","(SMGSOLV)"
    "","ALPB_Gsolv","(SMGSOLV)"



For Turbomole user:
-------------------

The amount of *ricore* for each calculation can be set in your `.cefinerc`. The same
holds for *maxcor* and/or *rpacor*.

.. code::

    $ echo "ricore  4000" > .cefinerc
    $ echo "maxcor  4000" >> .cefinerc
    $ echo "rpacor  4000" >> .cefinerc


Solvents:
---------

CENSO uses several QM-packages and not all solvents are available for all solvation
models throughout the QM-packages.
For this reason a user editable file is created in the folder:

    $  ~/.censo_assets/censo_solvents.json

which contains a dictionary of all available solvent models and solvents.
If a solvent is not available with a certain solvent model, the user can then choose
a replacement solvent. E.g. if CCl4 is not available choose CHCl3. 

.. figure:: docs/src/solvents.png
	:scale: 25%
	:align: center
	:alt: censo_solvents.json


The solvent file is directly used in `CENSO` and typos will cause calculations to crash!
Adding a new solvent is as easy as adding a new dictionary to the file.
