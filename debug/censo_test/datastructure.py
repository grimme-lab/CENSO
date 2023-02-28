"""
contains molecule_data class for storing all thermodyn. properties of the
conformer.
"""
from collections import OrderedDict
from .utilities import print
from .cfg import R, AU2KCAL, rot_sym_num
from math import log


# TODO - give information regarding status in the run in attributes e.g. above or below threshold
# TODO - devide into "mutable" and "immutable" attributes
# TODO - wtf is going on with this code
class MoleculeData:
    """
    molecule_data contains all thermodynamic properties of a conformer e.g.
    energy, gsolv, grrho
    """

    def __init__(self, id: int, xyz):
        """
        Creates a molecule instance where all thermodynamic data 
        concerning the molecule is stored.
        """
        


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
            "highlevel_grrho_info",
            "highlevel_hrrho_info",
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

    # def calc_free_energy(self, e=None, solv=None, rrho=None, t=None, out=False):
    #     """
    #     Calculate free energy for molecule either at normal temperature,
    #     or if the temperature is not None from the range of temperatures.
    #     if out=False free energy is written to self.free_energy
    #     if out=True free energy is simply returned
    #     """
    #     if t is None:
    #         try:
    #             f = 0.0
    #             if e is not None:
    #                 if e in ("xtb_energy", "xtb_energy_unbiased"):
    #                     f += getattr(self, e, 0.0)
    #                 else:
    #                     f += getattr(self, e, {"energy": 0.0})["energy"]
    #             if solv is not None:
    #                 f += getattr(self, solv, {"energy": 0.0})["energy"]
    #             if rrho is not None:
    #                 f += getattr(self, rrho, {"energy": 0.0})["energy"]
    #             if not out:
    #                 self.free_energy = f
    #             else:
    #                 return f
    #         except Exception as error:
    #             print("ERROR in _calc_free_energy: ", error)
    #             if not out:
    #                 self.free_energy = None
    #             else:
    #                 return f
    #     else:
    #         try:
    #             f = 0.0
    #             if e is not None:
    #                 if e in ("xtb_energy", "xtb_energy_unbiased"):
    #                     f += getattr(self, e, 0.0)
    #                 else:
    #                     f += getattr(self, e)["energy"]
    #             if solv is not None:
    #                 f += getattr(self, solv)["range"].get(t, 0.0)
    #             if rrho is not None:
    #                 f += getattr(self, rrho)["range"].get(t, 0.0)
    #             if not out:
    #                 self.free_energy = f
    #             else:
    #                 return f
    #         except (Exception, KeyError) as error:
    #             print("ERROR in _calc_free_energy: ", error)
    #             if not out:
    #                 self.free_energy = None
    #             else:
    #                 return f

    def calc_free_energy(
        self, e=None, solv=None, rrho=None, t=None, out=False, consider_sym=None
    ):
        """
        Calculate free energy for molecule either at normal temperature,
        or if the temperature is not None from the range of temperatures.
        if out=False free energy is written to self.free_energy
        if out=True free energy is simply returned 
        """
        if t is None:
            print(
                "xxxxxxxxxxxxx No temperature provided in calc_free_energy xxxxxxxxxxxxxxx"
            )
        try:
            f = 0.0
            if e is not None:
                if e in ("xtb_energy", "xtb_energy_unbiased"):
                    if getattr(self, e, 0.0) is not None:
                        f += getattr(self, e, 0.0)
                    else:
                        f += 0.0
                else:
                    f += getattr(self, e)["energy"]
            if solv is not None:
                if solv in ("cheap_prescreening_gsolv_info", "prescreening_gsolv_info"):
                    f += (
                        getattr(self, solv)
                        .get("range", {})
                        .get(t, getattr(self, solv, {"energy": 0.0})["energy"])
                    )
                else:
                    f += getattr(self, solv)["range"].get(t, 0.0)
            if rrho is not None:
                f += self.get_mrrho(
                    t, rrho=rrho, consider_sym=consider_sym, symnum=self.symnum
                )
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

    def _get_sym_num(self, sym=None, linear=None):
        """Get rotational symmetry number from Schoenfließ symbol"""
        if sym is None:
            sym = "c1"
        if linear is None:
            linear = False
        symnum = 1
        if linear and "c" in sym.lower()[0]:
            symnum = 1
            return symnum
        elif linear and "d" in sym.lower()[0]:
            symnum = 2
            return symnum
        for key in rot_sym_num.keys():
            if key in sym.lower():
                symnum = rot_sym_num.get(key, 1)
                break
        return symnum

    def calc_entropy_sym(self, temperature, symnum=None):
        """ RTln(sigma) rotational entropy"""
        if symnum is None:
            symnum = self.symnum
        return R / AU2KCAL * temperature * log(symnum)

    def get_mrrho(
        self, temperature, rrho=None, consider_sym=None, symnum=None, direct_input=0.0
    ):
        """ return mRRHO with or without RTln(sigma) (rot entropy)"""
        f = 0.0
        if rrho is not None:
            if rrho == "direct_input":
                f = direct_input
            elif rrho in ("prescreening_grrho_info",):
                f += (
                    getattr(self, rrho)
                    .get("range", {})
                    .get(temperature, getattr(self, rrho, {"energy": 0.0})["energy"])
                )
            elif rrho in ("rrho_optimization",):
                f += getattr(self, "optimization_info").get("energy_rrho", 0.0)
            else:
                f += getattr(self, rrho)["range"].get(temperature, 0.0)
            if consider_sym is None:
                consider_sym = True
            if not consider_sym:
                # sym is considered but consider_sym off
                if symnum is None:
                    symnum = self.symnum
                f += -self.calc_entropy_sym(temperature, symnum=symnum)
            return f
        else:
            return f
