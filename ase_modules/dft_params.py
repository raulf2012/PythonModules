"""DFT calculator's default parameters.

Development Notes:
  * TODO Automatically scale bands with system size
"""

#| - Import Modules
import os
import json
#__|

class DFT_Params:
    """Base class for DFT parameters encapsulation."""

    #| - DFT_Params ************************************************************

    def __init__(self):
        """Initialize base class."""
        #| - __init__
        self.file_name = "dft-params"
        self.compute_env = os.environ.get("COMPENV")
        #__|

    def load_params(self, dir=".", update_params=True):
        """Import JSON file containing parameters in dictionary.

        Args:
            dir:
            update_params:
        """
        #| - load_params
        try:
            data = open(dir + "/dft-params.json").read()
            data = json.loads(data)

            if update_params is True:
                self.update_params(data)
            else:
                self.params = data

        except:
            # data = None
            pass



        #| - OLD
        # try:
        #     # data=open(dir + "/" + name).read()
        #     data=open(dir + "/ase_modules/" + name).read()
        #     data = json.loads(data)
        # except:
        #     print("Could not load file")
        #     data = None
        #
        # return data
        #__|

        #__|

    def PythonPath(self):
        """Return PYTHONPATH system variable_lst.

        Checks dir is checked for default_espresso_params file
        """
        #| - PythonPath
        pyth_paths = os.environ["PYTHONPATH"].split(os.pathsep)


        pyth_path = [x for x in pyth_paths if "00_PythonModules" in x]

        #FIX | Assertion breaks whenver .bashrc is source multiple times, b.c.
        # redundant paths are duplicated in the PYTHONPATH variable

        # assert len(pyth_path) == 1, "Unique PYTHONPATH dir not found"

        return pyth_path[0]

        #__|

    def update_params(self, new_params, user_update=True):
        """Update parameter dict with new keys or new values for old keys.

        Args:
            new_params:
            user_update:
        """
        #| - update_params
        self.params.update(new_params)

        if user_update:
            mod_d_new = {x: True for x in new_params}
            self.mod_dict.update(mod_d_new)
        #__|

    def write_params(self, path=".", overwrite=False):
        """Write parameters to file in getcwd.

        # TODO change to path_i

        Args:
            path:
            overwrite:s
        """
        #| - write_params

        def json_dump_command(params, file):
            #| - json_dump_command
            json.dump(params, file, indent=2, skipkeys=True)
            #__|

        # num_files = [fle for fle in os.listdir(".") if "dft-params" in fle]
        num_files = [fle for fle in os.listdir(path) if "dft-params" in fle]
        num_files = len(num_files)

        if overwrite is False:
            if not os.path.exists(self.file_name + ".json"):
                with open(path + "/" + self.file_name + ".json", "w") as fle:
                    # json.dump(self.params, fle, indent=2)
                    json_dump_command(self.params, fle)

            else:
                fle_name = path + "/" + self.file_name + "_" + \
                    str(num_files) + ".json"
                with open(fle_name, "w") as fle:
                    json_dump_command(self.params, fle)
                    # json.dump(self.params, fle, indent=2)

        else:
            with open(path + "/" + self.file_name + ".json", "w+") as fle:
                json_dump_command(self.params, fle)
                # json.dump(self.params, fle, indent=2)
        #__|

    #__| ***********************************************************************

class VASP_Params(DFT_Params):
    """Useful method to define VASP parameters for DFT job."""

    #| - VASP_Params ***********************************************************
    def __init__(self, load_defaults=True):
        """TMP_docstring.

        TEMP TEMP
        """
        #| - __init__
        # self.pymoddir = self.PythonPath()
        DFT_Params.__init__(self)

        if load_defaults:
            self.params = self.default_params()
        else:
            self.params = {}

        self.mod_dict = self.create_mod_dict()

        # self.params = self.load_params(self.pymoddir,
        # "default_espresso_params.json")
        #__|

    def default_params(self):
        """User-defined default DFT parameters."""
        #| - default_params

        # Do calculation
        my_encut = 500
        my_enaug = 750
        my_kpts = (1, 1, 1)
        my_ediff = 0.00001
        my_nsw = 200  # 0 for ASE, > 0 for VASP internal
        my_prec = "Normal"  # "Normal" or "Accurate"
        my_istart = 0  # 0 = new job, 1:
        my_npar = 4
        my_ibrion = 1  # RMM-DIIS
        my_isif = 2
        my_fmax = -0.001  # ev/ang, tighter converergence
        my_ispin = 2

        params = {}

        params["ivdw"] = 0
        # params["ivdw"]      = 12

        # Dipole Correction
        params["ldipol"] = False
        # params["idipol"]    = 4
        # params[dipol]       = [0.5, 0.5, 0,5]  # Center of cell

        params["kpts"] = my_kpts
        params["encut"] = my_encut
        params["enaug"] = my_enaug

        params["ispin"] = my_ispin
        params["nsw"] = my_nsw
        params["prec"] = my_prec
        params["istart"] = my_istart
        params["isif"] = my_isif

        # IBRION
        params["ibrion"] = my_ibrion
        #

        params["ismear"] = 0
        params["sigma"] = 0.05

        params["icharg"] = 2
        params["lasph"] = True
        params["voskown"] = 1
        params["algo"] = "Normal"
        params["lreal"] = False

        # convergence Parameters
        params["ediff"] = my_ediff
        params["ediffg"] = my_fmax

        # Bader Charge Analysis
        params["laechg"] = False

        # Output Options
        params["nwrite"] = 1
        params["lcharg"] = True
        params["lwave"] = False

        # Symmetry
        params["isym"] = 0

        # Functional Type

        # Different VASP versions use 'pbe' vs 'PBE'
        if self.compute_env == "aws":
            pbe_val = "pbe"
        elif self.compute_env == "slac":
            pbe_val = "PBE"
        else:
            pbe_val = "PBE"
        # print(pbe_val)  # TEMP_PRINT
        params["xc"] = pbe_val  # TEMP
        params["gga"] = " PE"  # sets the PBE functiona

        # Compuational Parameters
        params["npar"] = my_npar
        params["nsim"] = 1
        params["nelmin"] = 4
        params["nelm"] = 100
        params["lplane"] = True

        # Dielectric Matrix density functional perturbation theory
        params["lepsilon"] = False

        return(params)
        #__|

    def create_mod_dict(self):
        """
        Create modification book-keepig dictionary.

        Dictionary which keeps track of whether parameter has been updated
        externally in script or if it kept at its default value.
        This is useful to know for the dependent-parameters.
        """
        #| - create_mod_dict
        param_keys = list(self.params.keys())
        mod_dict = dict((key, False) for key in param_keys)

        return(mod_dict)
        #__|


    #__| ***********************************************************************


class Espresso_Params(DFT_Params):
    """Useful method to define quantum espresso parameters for DFT job."""

    #| - Espresso_Params *******************************************************
    def __init__(self, load_defaults=True):
        """TMP_docstring.

        TEMP TEMP
        """
        #| - __init__
        # self.pymoddir = self.PythonPath()
        DFT_Params.__init__(self)

        if load_defaults:
            self.params = self.default_params()
        else:
            self.params = {}

        self.mod_dict = self.create_mod_dict()

        # self.params = self.default_params()
        # self.params = self.load_params(self.pymoddir,
        # "default_espresso_params.json")
        #__|

    def create_mod_dict(self):
        """
        Create modification book-keepig dictionary.

        Dictionary which keeps track of whether parameter has been updated
        externally in script or if it kept at its default value.
        This is useful to know for the dependent-parameters.
        """
        #| - create_mod_dict
        param_keys = list(self.params.keys())
        mod_dict = dict((key, False) for key in param_keys)

        return(mod_dict)
        #__|

    def default_params(self):  # ***********************************************
        """User-defined default DFT parameters."""
        #| - default_params
        params = {}
        params["pw"] = 500            # plane-wave cutoff
        params["dw"] = 5000           # density cutoff
        params["dipole"] = {"status": True}    # Turn on only for slabs not bulk
        params["xc"] = "BEEF-vdW"     # exchange-correlation functional
        params["kpts"] = (3, 3, 1)    # k-points for hexagonal symm in 2-D mater
        # 20 extra bands besides bands needed for val elec
        params["nbands"] = -20
        params["spinpol"] = False     # Spin-polarized calculation

        params["sigma"] = 0.1  # Should be low for spin calculations
        # params["psppath"] = "/home/vossj/suncat/psp/gbrv1.5pbe/"
        # pseudopotential path

        params["beefensemble"] = False

        # Parallelization <----------------------------------------------------
        # params["parflags"] = "-npool "
        params["parflags"] = None

        # Convergence Criteria <----------------------------------------------
        params["convergence"] = {
            "energy": 1e-5,  # convergence parameters

            # TF (Normal) or local-TF (Better for inhomogeneous systems)
            "mixing_mode": "local-TF",
            "mixing": 0.2,
            "nmix": 20,
            "maxsteps": 500,
            "diag": "david",
            # "mix": 4,  # This parameter doesn't do anything (Pretty sure)
            }


        # File Output <--------------------------------------------------------
        params["output"] = {"removesave": True}  # Aayush, saves almost nothing
        # params["output"] = {
        #     "avoidio": False,
        #     "removewf": True,
        #     "wf_collect": False
        #     }

        params["outdir"] = "calcdir"

        return(params)
        #__|

    def test_check(self):
        """
        Automatically tries to set params based on other dependent parameters.

        Ex.) If spinpol == True, then sigma should be smaller, 0.01ish
        (Charlotte stold me)

        Only set this dependent param if the user has not explicitly done so,
        is what the condition for setting a parameter automatically is.
        """
        #| - test_check

        #| - Setting dw to 10 * pw
        if self.mod_dict["dw"] is False:
            dw_i = 10. * self.params["pw"]
            self.update_params({"dw": dw_i}, user_update=False)
        #__|

        #| - Decreasing Sigma for Spin Polarization Calculations
        if self.mod_dict["sigma"] is False:
            if self.params["spinpol"] is True:
                sigma_i = 0.02
                self.update_params({"sigma": sigma_i}, user_update=False)
        #__|

        #| - Removing Beef-Ensemble if XC-Functional Not BEEF
        xc_list = self.params["xc"]
        if "beefensemble" in self.params:
            if self.params["beefensemble"] is True and "BEEF" not in xc_list:
                print("Functional not compatible with BEEF-ensemble method")
                self.update_params({"beefensemble": False}, user_update=False)
        else:
            pass
        #__|

        #__|

    #__| ***********************************************************************
