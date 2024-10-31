from .configuration import configure
from .params import DESCR
from .__version__ import __version__

print(DESCR)
configure()

from .cli import interface, cml_parser
from . import (
    configuration,
    ensembledata,
    datastructure,
    orca_processor,
    parallel,
    part,
    qm_processor,
    utilities,
    ensembleopt,
    properties,
)
