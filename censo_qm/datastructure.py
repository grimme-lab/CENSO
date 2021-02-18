"""
contains molecule_data class for storing all thermodyn. properties of the
conformer.
"""
from collections import OrderedDict
from .utilities import print


class MoleculeData:
    """
    molecule_data contains all thermodynamic properties of a conformer e.g.
    energy, gsolv, grrho
    """

    def __init__(
        self,
        rank,
        chrg=0,
        uhf=0,
        xtb_energy=None,
        xtb_energy_unbiased=None,
        xtb_free_energy=None,
        rel_xtb_energy=None,
        rel_xtb_free_energy=None,
        sym="c1",
        gi=1.0,
        removed=False,
        free_energy=0.0,
        temperature_info={"temperature": 298.15, "range": None},
        cheap_prescreening_sp_info={
            "energy": None,
            "info": "not_calculated",
            "method": None,
            "prev_methods": None,
        },
        cheap_prescreening_gsolv_info={
            "energy": None,
            "gas-energy": None,
            "solv-energy": None,
            "info": "not_calculated",
            "method": None,
            "prev_methods": None,
        },
        prescreening_sp_info={
            "energy": None,
            "info": "not_calculated",
            "method": None,
            "prev_methods": None,
        },
        lowlevel_sp_info={
            "energy": None,
            "info": "not_calculated",
            "method": None,
            "prev_methods": None,
        },
        highlevel_sp_info={
            "energy": None,
            "info": "not_calculated",
            "method": None,
            "prev_methods": None,
        },
        prescreening_grrho_info={
            "energy": None,
            "info": "not_calculated",
            "method": None,
            "fuzzythr": 0.0,
            "rmsd": None,
            "prev_methods": None,
        },
        lowlevel_grrho_info={
            "energy": None,
            "range": None,
            "info": "not_calculated",
            "method": None,
            "rmsd": None,
            "prev_methods": None,
        },
        lowlevel_hrrho_info={
            "energy": None,
            "range": None,
            "info": "not_calculated",
            "method": None,
            "rmsd": None,
            "prev_methods": None,
        },
        highlevel_grrho_info={
            "energy": None,
            "range": None,
            "info": "not_calculated",
            "method": None,
            "rmsd": None,
            "prev_methods": None,
        },
        highlevel_hrrho_info={
            "energy": None,
            "range": None,
            "info": "not_calculated",
            "method": None,
            "rmsd": None,
            "prev_methods": None,
        },
        prescreening_gsolv_info={
            "energy": None,
            "gas-energy": None,
            "info": "not_calculated",
            "method": None,
            "prev_methods": None,
        },
        lowlevel_gsolv_info={
            "energy": None,
            "gas-energy": None,
            "range": None,
            "info": "not_calculated",
            "method": None,
            "prev_methods": None,
        },
        lowlevel_gsolv_compare_info={
            "energy": None,
            "range": None,
            "info": "not_calculated",
            "method": None,
            "prev_methods": None,
            "std_dev": None,
        },
        highlevel_gsolv_info={
            "energy": None,
            "gas-energy": None,
            "range": None,
            "info": "not_calculated",
            "method": None,
            "prev_methods": None,
        },
        optimization_info={
            "energy": None,
            "convergence": "not_converged",
            "cregen_sort": "pass",  # pass and removed
            "info": "not_calculated",
            "cycles": 0,
            "ecyc": [],
            "decyc": [],
            "energy_rrho": 0.0,
            "method_rrho": None,
            "info_rrho": "not_calculated",
        },
        nmr_coupling_info={
            "info": "not_calculated",
            "method": None,
            "h_active": False,
            "c_active": False,
            "f_active": False,
            "si_active": False,
            "p_active": False,
        },
        nmr_shielding_info={
            "info": "not_calculated",
            "method": None,
            "h_active": False,
            "c_active": False,
            "f_active": False,
            "si_active": False,
            "p_active": False,
        },
        part_info={
            "part0": None,
            "part1": None,
            "part2": None,
            "part3": None,
            "part4": None,
            "part5": None,
        },
        comment=[],
        optical_rotation_info={
            "range": None,
            "info": "not_calculated",
            "method": None,
            "prev_methods": None,
        },
    ):
        """
        molecule_data: Creates a molecule instance where all thermodynamic data 
                       concerning the molecule is stored.
        Input:
         - rank [int] --> rank of the molecule in the input ensemble e.g. CONF(X)
         - temperature [float] --> evaluation at this temperature
         - trange [list(float)] --> list of temperatures for evaluation at 
                                multiple temperatures
         - chrg [int] --> charge of the molecule
         - uhf  [int] --> number of unpaired electrons
         - xtb_energy [float] a.u.-->  initial energy taken from the input ensemble
         - rel_xtb_energy [float] kcal/mol --> relative initial energy taken 
                                        from the input ensemble
         - sym [string] --> schoenflies notation of pointgroup
         - gi  [float]  --> degeneracy of conformer 
         - prescreening_sp [float] --> single point energy of the preescreening 
                                        procecure
         - lowlevel_sp  [float] --> single point energy of the optimization
         - highlevel_sp  [float] --> high level single point energy at hybrid 
                                    level with larger basis set
         - prescreening_grrho [float] --> thermostatistical contribution of 
                                            prescreening procedure
         - grrho [float] --> thermostatistical (G) contribution on optimized DFT
                             geometry
         - hrrho [float] --> thermostatistical (H) contribution on optimized DFT
                             geometry
         - part_info [string] --> partx : passed/refused/not_calculated
         

        *_info = {'info': calculated/not_calculated/failed/skipped/removed/prep-failed}
        """
        # check default arguments:
        for store in (
            prescreening_grrho_info,
            lowlevel_grrho_info,
            lowlevel_hrrho_info,
            highlevel_grrho_info,
            highlevel_hrrho_info,
            prescreening_gsolv_info,
            lowlevel_gsolv_info,
            lowlevel_gsolv_compare_info,
            highlevel_gsolv_info,
            optical_rotation_info,
            prescreening_sp_info,
            lowlevel_sp_info,
            highlevel_sp_info,
            cheap_prescreening_sp_info,
            cheap_prescreening_gsolv_info,
        ):
            if store.get("prev_methods", None) is None:
                store["prev_methods"] = {}
        if temperature_info.get("range") is None:
            temperature_info["range"] = []
        if xtb_energy is None:
            xtb_energy = 100.0
        if xtb_energy_unbiased is None:
            xtb_energy_unbiased = 100.0
        if xtb_free_energy is None:
            xtb_free_energy = 100.0
        if rel_xtb_energy is None:
            rel_xtb_energy = 100.0
        if rel_xtb_free_energy is None:
            rel_xtb_free_energy = 100.0
        if cheap_prescreening_sp_info.get("energy") is None:
            cheap_prescreening_sp_info["energy"] = 0.0
        if cheap_prescreening_gsolv_info.get("energy") is None:
            cheap_prescreening_gsolv_info["energy"] = 0.0
        if prescreening_sp_info.get("energy") is None:
            prescreening_sp_info["energy"] = 0.0
        if lowlevel_sp_info.get("energy") is None:
            lowlevel_sp_info["energy"] = 0.0
        if highlevel_sp_info.get("energy") is None:
            highlevel_sp_info["energy"] = 0.0
        if prescreening_grrho_info.get("energy") is None:
            prescreening_grrho_info["energy"] = 0.0
        if lowlevel_grrho_info.get("energy") is None:
            lowlevel_grrho_info["energy"] = 0.0
        self._initialize(lowlevel_grrho_info)
        if lowlevel_hrrho_info.get("energy") is None:
            lowlevel_hrrho_info["energy"] = 0.0
        self._initialize(lowlevel_hrrho_info)
        if prescreening_gsolv_info.get("energy") is None:
            prescreening_gsolv_info["energy"] = 0.0
        if lowlevel_gsolv_info.get("energy") is None:
            lowlevel_gsolv_info["energy"] = 0.0
        self._initialize(lowlevel_gsolv_info)
        if lowlevel_gsolv_compare_info.get("energy") is None:
            lowlevel_gsolv_compare_info["energy"] = 0.0
        # if lowlevel_gsolv_compare_info.get("std_dev") is None:
        #     lowlevel_gsolv_compare_info["std_dev"] = 0.0
        self._initialize(lowlevel_gsolv_compare_info)
        if highlevel_gsolv_info.get("energy") is None:
            highlevel_gsolv_info["energy"] = 0.0
        self._initialize(highlevel_gsolv_info)
        for key in part_info.keys():
            if part_info.get(key) is None:
                part_info[key] = "not_calculated"
        # highlevel_grrho
        if highlevel_grrho_info.get("energy") is None:
            highlevel_grrho_info["energy"] = 0.0
        self._initialize(highlevel_grrho_info)
        # highlevel_hrrho
        if highlevel_hrrho_info.get("energy") is None:
            highlevel_hrrho_info["energy"] = 0.0
        self._initialize(highlevel_hrrho_info)
        # optical_rotation_info
        self._initialize(optical_rotation_info)

        # exceptions:
        if not isinstance(rank, int):
            raise TypeError(
                "Please input an integer. The id is the rank of the "
                "molecule in the input ensemble!"
            )
        if type(temperature_info.get("temperature", None)) != float:
            raise TypeError(
                "Please input an float. Thermodynamic properties are"
                "evaluated at this temperature!"
            )
        if not isinstance(temperature_info.get("range", None), list):
            raise ValueError("Please provide a list with temperatures!")
        elif any([type(i) != float for i in temperature_info.get("range")]):
            raise TypeError("Please provide temperatures as float!")
        if not isinstance(chrg, int):
            raise TypeError("Please provide charge as integer!")
        if not isinstance(uhf, int):
            raise TypeError(
                "Please provide number of unpaired electrons as " "integer!"
            )
        if not isinstance(xtb_energy, float):
            raise TypeError("Please provide energy from input ensemble as float!")
        if not isinstance(rel_xtb_energy, float):
            raise TypeError(
                "Please provide rel. energy from input ensemble as " "float!"
            )
        if not isinstance(sym, str):
            raise TypeError("Please provide symmetry as string.")
        if not isinstance(gi, float):
            try:
                gi = float(gi)
            except (TypeError, ValueError):
                raise "Please provide gi as float!"
        if not isinstance(prescreening_sp_info.get("energy", None), float):
            raise TypeError("Please provide preescreening sinlge point as float!")
        if not isinstance(lowlevel_sp_info.get("energy", None), float):
            raise TypeError("Please provide low level sinlge point as float!")
        if not isinstance(highlevel_sp_info.get("energy", None), float):
            raise TypeError("Please provide high level sinlge point as float!")
        if not isinstance(prescreening_grrho_info.get("energy", None), float):
            raise TypeError("Please provide G_RRHO as float!")
        if type(lowlevel_grrho_info.get("energy", None)) != float:
            raise TypeError("Please provide G_RRHO as float!")
        if not isinstance(lowlevel_grrho_info["range"], dict):
            raise TypeError("Please input a dict with Grrho values!")
        if not isinstance(lowlevel_hrrho_info.get("energy", None), float):
            raise TypeError("Please provide H_RRHO as float!")
        if any([type(i) != float for i in lowlevel_grrho_info.get("range").values()]):
            raise TypeError("Please provide Grrho values as float!")
        if not isinstance(lowlevel_hrrho_info["range"], dict):
            raise TypeError("Please input a dict with Hrrho values!")
        if any([type(i) != float for i in lowlevel_hrrho_info.get("range").values()]):
            raise TypeError("Please provide Hrrho values as float!")
        if not isinstance(prescreening_gsolv_info.get("energy"), float):
            raise TypeError("Please provide Gsolv as float!")
        if not isinstance(lowlevel_gsolv_info.get("energy", None), float):
            raise TypeError("Please provide Gsolv as float!")
        if not isinstance(lowlevel_gsolv_info.get("range"), dict):
            raise TypeError("Please input a dict with Gsolv values!")
        if any([type(i) != float for i in lowlevel_gsolv_info.get("range").values()]):
            raise TypeError("Please provide Gsolv values as float!")
        if type(highlevel_gsolv_info.get("energy", None)) != float:
            raise TypeError("Please provide Gsolv as float!")
        if not isinstance(highlevel_gsolv_info.get("range", None), dict):
            raise TypeError("Please input a dict with Gsolv values!")
        if any([type(i) != float for i in highlevel_gsolv_info.get("range")]):
            raise TypeError("Please provide Gsolv values as float!")
        if not isinstance(removed, bool):
            raise TypeError("Please provide removed with boolean true/false.")
        if any([type(i) != str for i in part_info.values()]):
            raise TypeError("Please provide part_info settings as str!")
        # assignment:
        self.id = rank  # this is the rank from the input ensemble
        self.temperature_info = temperature_info  # temperature for general evaluation
        self.chrg = chrg
        self.uhf = uhf
        self.xtb_energy = xtb_energy
        self.xtb_energy_unbiased = xtb_energy_unbiased
        self.xtb_free_energy = xtb_free_energy
        self.rel_xtb_energy = rel_xtb_energy
        self.rel_xtb_free_energy = rel_xtb_free_energy
        self.sym = sym
        self.gi = gi
        self.cheap_prescreening_gsolv_info = cheap_prescreening_gsolv_info
        self.cheap_prescreening_sp_info = cheap_prescreening_sp_info
        self.prescreening_sp_info = prescreening_sp_info
        self.lowlevel_sp_info = lowlevel_sp_info
        self.highlevel_sp_info = highlevel_sp_info
        self.prescreening_grrho_info = prescreening_grrho_info
        self.lowlevel_grrho_info = lowlevel_grrho_info
        self.lowlevel_hrrho_info = lowlevel_hrrho_info
        self.highlevel_grrho_info = highlevel_grrho_info
        self.highlevel_hrrho_info = highlevel_hrrho_info
        self.prescreening_gsolv_info = prescreening_gsolv_info
        self.lowlevel_gsolv_info = lowlevel_gsolv_info
        self.lowlevel_gsolv_compare_info = lowlevel_gsolv_compare_info
        self.highlevel_gsolv_info = highlevel_gsolv_info
        self.optimization_info = optimization_info
        self.nmr_coupling_info = nmr_coupling_info
        self.nmr_shielding_info = nmr_shielding_info
        self.removed = removed
        self.free_energy = free_energy
        self.part_info = part_info
        self.comment = comment
        self.optical_rotation_info = optical_rotation_info

    def _initialize(self, attr=None):
        """
        json saves keys as string. Convert some keys to float.
        """
        if attr is not None:
            if attr.get("range") is None:
                attr["range"] = {}
            else:  # check if keys are float
                if isinstance(attr["range"], dict):
                    new = {}
                    for key, value in attr["range"].items():
                        new[float(key)] = value
                    attr["range"] = new
            for method in attr["prev_methods"]:
                if isinstance(attr["prev_methods"][method].get("range"), dict):
                    new = {}
                    for key, value in attr["prev_methods"][method]["range"].items():
                        new[float(key)] = value
                    attr["prev_methods"][method]["range"] = new
                else:
                    attr["prev_methods"][method]["range"] = {}

    def reset_range_info(self, trange=None):
        """
        Reset all dictionaries concerned with a temperature range and
        set info to not calculated. (This is needed if the temperature range was
        not calculated in a previous run and is requested in a current run).
        trange -> list with temperatures
        """
        attributes = [
            "lowlevel_grrho_info",
            "lowlevel_hrrho_info",
            "lowlevel_gsolv_info",
            "highlevel_gsolv_info",
        ]
        # reset only if not all temperatures are found which are needed in trange
        for data in attributes:
            reset = False
            if trange is not None:
                for temp in trange:
                    if getattr(self, data)["range"].get(temp, None) is None:
                        reset = True
            if reset:
                getattr(self, data)["info"] = "not_calculated"
                # keep only value at "normal" temperature
                getattr(self, data)["range"] = {
                    self.temperature_info["temperature"]: getattr(self, data)["energy"]
                }
        # END---

    def save_prev(self, attr, method):
        """
        save dictionary with all information of
        previously calculated data
        """
        # store data under 'prev_methods'[method]
        tmp = {method: {}}
        attributes = vars(MoleculeData(0)).get(attr)
        if getattr(self, attr)["info"] != "not_calculated":
            for key in attributes.keys():
                if key != "prev_methods":
                    tmp[method][key] = getattr(self, attr).get(key)
            getattr(self, attr)["prev_methods"].update(tmp)

    def load_prev(self, attr, method, saveto=None):
        """
        load dictionary with all information from
        previously calculated data,
        if not previously calculated load presets
        self --> conf object
        attr --> e.g. lowlevel_sp_info
        method --> method identifier e.g. func/basis[sm]
        saveto --> optional if desired to save data somewhere else
                   e.g. highlevel_sp_info
        """
        attributes = vars(MoleculeData(0)).get(attr)
        # check if calculated previously
        if getattr(self, attr)["prev_methods"].get(method, None) is not None:
            tmp = {}
            for key in attributes.keys():
                if key != "prev_methods":
                    tmp[key] = getattr(self, attr)["prev_methods"][method].get(key)
            if saveto is not None:
                attr = saveto
            getattr(self, attr).update(tmp)
        else:
            # if not calculated previously reset
            for key, value in attributes.items():
                if key != "prev_methods":
                    getattr(self, attr)[key] = value

    def provide_runinfo(self):
        """
        Write dictionary with molecule data information:
        """
        runinfo = []
        for key in vars(MoleculeData(0)).keys():
            runinfo.append((key, getattr(self, key)))
        return OrderedDict(runinfo)

    def calc_free_energy(self, e=None, solv=None, rrho=None, t=None, out=False):
        """
        Calculate free energy for molecule either at normal temperature,
        or if the temperature is not None from the range of temperatures.
        if out=False free energy is written to self.free_energy
        if out=True free energy is simply returned 
        """
        if t is None:
            try:
                f = 0.0
                if e is not None:
                    if e in ("xtb_energy", "xtb_energy_unbiased"):
                        f += getattr(self, e, 0.0)
                    else:
                        f += getattr(self, e, {"energy": 0.0})["energy"]
                if solv is not None:
                    f += getattr(self, solv, {"energy": 0.0})["energy"]
                if rrho is not None:
                    f += getattr(self, rrho, {"energy": 0.0})["energy"]
                if not out:
                    self.free_energy = f
                else:
                    return f
            except Exception as error:
                print("ERROR in _calc_free_energy: ", error)
                if not out:
                    self.free_energy = None
                else:
                    return f
        else:
            try:
                f = 0.0
                if e is not None:
                    if e in ("xtb_energy", "xtb_energy_unbiased"):
                        f += getattr(self, e, 0.0)
                    else:
                        f += getattr(self, e)["energy"]
                if solv is not None:
                    f += getattr(self, solv)["range"].get(t, 0.0)
                if rrho is not None:
                    f += getattr(self, rrho)["range"].get(t, 0.0)
                if not out:
                    self.free_energy = f
                else:
                    return f
            except (Exception, KeyError) as error:
                print("ERROR in _calc_free_energy: ", error)
                if not out:
                    self.free_energy = None
                else:
                    return f
