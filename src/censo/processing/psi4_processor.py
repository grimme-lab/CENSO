"""
Contains Psi4Proc class for calculating psi4 related properties of conformers.
"""

import typing

from ..utilities import Factory
from ..params import (
    Prog,
)
from .qm_processor import QmProc


@typing.final
class Psi4Proc(QmProc):
    pass


Factory.register_builder(Prog.PSI4, Psi4Proc)
