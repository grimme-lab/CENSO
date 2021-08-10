|GitHub release| |DOI|

.. |GitHub release| image:: https://img.shields.io/github/v/release/grimme-lab/CENSO
   :target: https://github.com/grimme-lab/CENSO/releases/latest

.. |DOI| image:: https://img.shields.io/badge/DOI-10.1021/acs.jpca.1c00971-blue
    :target: https://doi.org/10.1021/acs.jpca.1c00971


====================================================================
CENSO - Commandline ENergetic SOrting of Conformer Rotamer Ensembles
====================================================================

.. raw:: html

    <p align="center">
    <img src="docs/src/censo_logo_300dpi.png" alt="CENSO logo" width="500">
    </p>

This repository hosts the ``CENSO`` code for the refinement of Conformer Rotamer 
Ensembles (CRE) as obtained from ``CREST``.


Installation
============

There are several options possible. The easiest is to use the packaged censo programs
(by use of Pyinstaller) which can be found at the release section. The packaged 
censo is linked against GLIBC version 2.19 and will work for GLIBC version 2.19 and above.

For the packaged "binary" of CENSO, download it from the 
`release <https://github.com/grimme-lab/CENSO/releases/>`_ site, 
copy to your bin and make executable:

.. code::

    $ cp censo ~/bin/censo
    $ chmod u+x ~/bin/censo

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
================

Create the remote configuration file .censorc where the user can adjust default
settings and provide paths to the external programs e.g. `xtb`, `crest`, `orca` ...

.. code::

    $ censo -newconfig
    $ cp censorc-new /home/$USER/.censorc
    # edit .censorc
    vi /home/$USER/.ensorc


Documentation:
==============

Can be found following: https://xtb-docs.readthedocs.io/en/latest/CENSO_docs/censo.html

**Interactive Documentation can be accessed:**

.. code::

    $ censo -tutorial

Explainations on the commandline arguments can be printed by:

.. code::

    $ censo --help


Requirements:
-------------

* newest xtb >= v6.4.0
* newest cefine https://github.com/grimme-lab/cefine/releases
* ORCA > version 4.1 


Further information (will be ordered later on):

* the file .censorc can be used in the current working directory and will be preferred to 
  the global configuration file in ~/.censorc
* a folder ~/.censo_assets/ will be created upon usage of censo
* ORCA has not been tested extensively so please be careful, test calculations
  and report possible "bad" settings
* To save computational cost COSMO-RS calculations are not performed with BP86 but whith the 
  density functionals for energy evaluation.


Usage:
------

.. code::

    # check if setting-combinations match:
    $ censo -inp structure_ensemble.xyz -part2 on -solvent h2o --checkinput
    # start the calculation:
    $ censo -inp structure_ensemble.xyz -part2 on -solvent h2o > censo.out 2> error.censo &


Short overview:
---------------

The molecule numbering from the input structure ensemble is kept throughout the 
entire program. There are several program parts which can be used to filter a structure 
ensemble:

1) Cheap prescreening (part0): Very fast DFT energies in order to improve upon the energy
    description of the SQM method used to generate the input structure ensemble.
    The (free) energies are evaluated on the input geometries (DFT unoptimized).

2) Prescreening (part1): Improved DFT energies and accurate solvation energies (if needed).
    The free energies are evaluated on the input geometries (DFT unoptimized).

3) Optimization (part2): efficient structure ensemble optimization and 
    free energy calculation on DFT optimized geometries.

4) Refinement (part3): Optional free energy refinement (on DFT optimized geometries).
    E.g. using hybrid DFA with large basis set.

5) NMR properties (part4): Optional calculation of shielding and coupling constants on 
    populated conformers.

6) Optical Rotation (part5): Optional calculation of optical rotatory dispersion 
    for the populated ensemble.


For Turbomole user:
-------------------

The amount of *ricore* for each calculation can be set in your `.cefinerc`. The same
holds for *maxcor* and/or *rpacor*.

.. code::

    $ echo "ricore  4000" > .cefinerc
    $ echo "maxcor  4000" >> .cefinerc
    $ echo "rpacor  4000" >> .cefinerc

For ORCA user:
--------------

CENSO has been updated to work with the new ORCA 5.0.1 release. The user has to provide the matching
ORCA version number in the `.cefinerc` file in order for the correct ORCA input generation 
to work, e.g. 

.. code::

    ORCA: /tmp1/orca_5_0_1_linux_x86-64_openmpi411
    ORCA version: 5.0.1


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


Solvents:
---------

CENSO uses several QM-packages and not all solvents are available for all solvation
models throughout the QM-packages.
For this reason a user editable file is created in the folder:

    $  ~/.censo_assets/censo_solvents.json

which contains a dictionary of all available solvent models and solvents.
If a solvent is not available with a certain solvent model, the user can then choose
a replacement solvent. E.g. if CCl4 is not available choose CHCl3. 


.. raw:: html

    <p align="center">
    <img src="docs/src/solvents.png" alt="censo_solvents.json" width="700">
    </p>

The solvent file is directly used in `CENSO` and typos will cause calculations to crash!
Adding a new solvent is as easy as adding a new dictionary to the file.

Cite
----

General reference:

S. Grimme, F. Bohle, A. Hansen, P. Pracht, S. Spicher, and M. Stahn 
*J. Phys. Chem. A* **2021**, *125* (19), 4039–4054.

DOI: `10.1021/acs.jpca.1c00971 <https://doi.org/10.1021/acs.jpca.1c00971>`_. 

Reference is available in `bibtex format <./docs/reference.bib>`_.


License
-------

``CENSO`` is free software: you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

``CENSO`` is distributed in the hope that it will be useful,
but without any warranty; without even the implied warranty of
merchantability or fitness for a particular purpose. See the
GNU Lesser General Public License for more details.
