from typing import Dict

from censo.qm_processor import QmProc
from censo.orca_processor import OrcaProc
# from censo.tm_processor import TmProc

class ProcessorFactory:
    
    # for now these are the only available processor types
    __proctypes: Dict[str, type] = {
        "orca": OrcaProc,
    #    "tm": TmProc,
    }
    
    @classmethod
    def create_processor(cls, prog, *args, **kwargs) -> QmProc:
        """
        returns an instance of the requested processor type (mapping with the 'prog' setting)
        for now the QmProc class uses xtb as driver (this method should be changed if additional drivers are implemented)
        available processor types are mapped in ProcessorFactory.__proctypes

        example: create_processor(orca, external_paths, solvents_dict=...) 
        (lookup the constructors for the processor types for further documentation)
        """
        type_t: type = cls.__proctypes.get(prog, None)
        
        if not type_t is None:
            return type_t(*args, **kwargs)
        else:
            raise TypeError(f"No processor type was found for {prog} in {cls.__proctypes}.")