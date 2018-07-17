#!/usr/bin/env python

"""Methods for ASE scripts, mostly DFT scripts.

Author: Raul A. Flores

Development Notes:
    TODO Master print statements should have "****..." to be easily readable
    TODO Delete all set_mag_mom_to_0 references (Already commented them out)
    TODO Add work function analysis (Ask Colin about Karen's method)
"""

#| - IMPORT MODULES
import sys
import os
import glob
import copy
import json
import math
import pickle as pickle
import numpy as np
import shutil

from scipy.stats import norm

# NOTE Changed for python3.6
# import matplotlib.pyplot as plt
import matplotlib as plt

from ase.io import read, write, Trajectory
from ase import io
from ase.dft.kpoints import ibz_points, get_bandpath
from ase.vibrations import Vibrations
from ase.thermochemistry import HarmonicThermo

# pymatgen
from pymatgen.core.periodic_table import _pt_data as periodic_table_dict

# My Modules
from misc_modules.numpy_methods import angle_between
from ase_modules.dft_params import Espresso_Params

from quantum_espresso.qe_methods import estimate_magmom

from bader_charge.bader import bader
#__|

#| - METHODS

def update_FINISHED(text, filename=".FINISHED.new"):
    """Update finished job/processes log file with message.

    TODO Change filename --> .FINISHED when the migration to new system is
    complete

    Args:
        text:
        filename:
    """
    #| - update_FINISHED
    if os.path.exists("./" + filename):
        append_write = "a"  # append if already exists
    else:
        append_write = "w"  # make a new file if not

    with open(filename, append_write) as fle:
        fle.write(text)
        fle.write("\n")
    #__|

#__|

#| - Parse DFT Job Parameters *************************************************

def set_QE_calc_params(
    atoms=None,
    params={},
    load_defaults=True,
    init_inst=True,
    ):
    """Set Quantum Espresso calculation parameters to atoms object.

    Handles reading, and setting of dft calculator parameters to atoms object.

    Args:
        atoms:
        params:
        load_defaults:
        init_inst:
            Whether or not to institiate an espresso instance
    """
    #| - set_QE_calc_params
    from espresso import espresso

    mess = "Loading QE Parameters From File "
    mess += "**********************************************"
    print(mess); sys.stdout.flush()

    espresso_params_inst = Espresso_Params(load_defaults=load_defaults)

    if os.path.isfile("dft-params.json"):
        params_file = json.load(open("dft-params.json"))
        espresso_params_inst.update_params(params_file)

    elif os.path.isfile("dir_dft_params/dft-params.json"):
        params_file = json.load(open("dir_dft_params/dft-params.json"))
        espresso_params_inst.update_params(params_file)

    espresso_params_inst.update_params(params)  # Create params dict if wanted

    espresso_params_inst.test_check()
    espresso_params_inst.write_params()
    espresso_params = espresso_params_inst.params

    calc = None
    if init_inst:
        calc = espresso(**espresso_params)

    # NOTE Calculator set for the 1st time here
    # atoms.set_calculator(calc=calc)

    return(calc, espresso_params)
    #__|

#__| **************************************************************************

#| - Ionic Optimization *******************************************************

def ionic_opt(
    atoms,
    calc=None,
    traj=None,
    espresso_params=None,
    mode="opt",
    fmax=0.05,
    maxsteps=100000000,
    run_beef_an=True,
    run_bader_an=True,
    ):
    """Run ionic dft relaxation on atoms object.

    Development Notes:
        * What to do when the job has already been previously completed?
        Should I run a single point calculation just to make sure the
        calculator is fully operational/loaded?

        TODO .FINSISHED file should include information about what kind of
        optimization was performed

        TODO Implemetn qn.replay_trajectory('qn.traj') functionality

    Args:
        atoms: ASE atoms object
        calc: Calculator object
        espresso_params: Quantum Espresso job parameters dictionary
        mode:
        fmax: Force convergence criteria
        maxsteps: Maximum number of ionic steps
        run_beef_an:
            Attempts to run Beef-vdW ensemble of energies
            Must have ran calculation with appropriate paramters to begin with
    """
    #| - ionic_opt
    from espresso import espresso
    from ase.optimize import QuasiNewton

    mess = "Running DFT Calculation "
    mess += "******************************************************"
    print(mess)
    sys.stdout.flush()

    # TODO Skip calculation if previously converged
    # COMBAK Uncomment this section
    #| - Checking if Previous Calculation Has Been Completed
    # filename = ".FINISHED.new"
    # if os.path.exists("./" + filename):
    #     with open(filename, "r") as fle:
    #         lines = [line.strip() for line in fle.readlines()]
    #         if "ionic_opt" in lines:
    #             print("ionic_opt | Optimization previusly completed "
    #                 "(Running a single-point calculation)"
    #                 )
    #             mode = "sp"
    #__|

    #| - Setting Optimization Specific Espresso Parameters

    # espresso_params_copy = copy.deepcopy(espresso_params)

    params_opt = {
        "output": {
            "avoidio": True,
            "removesave": True,
            "removewf": True,
            "wf_collect": False,
            },

        "outdir": "calcdir_opt",
        }

    # espresso_params_copy.update(params_opt)
    # calc_opt = espresso(**espresso_params_copy)
    # atoms.set_calculator(calc_opt)

    calc_opt, espresso_params_opt = set_QE_calc_params(
        params=params_opt,
        )

    calc_opt = espresso(**espresso_params_opt)
    atoms.set_calculator(calc_opt)
    #__|

    reduce_magmoms(atoms)

    if mode == "opt":
        #| - Regular Optimization
        mess = "Running regular optimization "
        print(mess); sys.stdout.flush()

        qn = QuasiNewton(
            atoms,
            # trajectory="out_opt.traj",
            logfile="qn.log",
            )

        if traj is not None:
            qn.attach(traj)  # COMBAK Test feature (restarting traj files)

        if os.path.exists("prev.traj"):
            qn.replay_trajectory("prev.traj")

        qn.run(
            fmax=fmax,
            steps=maxsteps,
            )
        #__|

    elif mode == "easy_opt":
        #| - Easy Optimization -> Full Optimization
        mess = "Running easy optimization scheme"
        print(mess); sys.stdout.flush()

        #| - Easy Optimization Settings
        espresso_params_copy = copy.deepcopy(espresso_params)

        easy_params = {
            "pw": 400,
            "dw": 4000,
            "spinpol": False,
            "kpts": (3, 3, 1),

            "convergence": {
                "energy": 5e-5,  # convergence parameters
                "mixing": 0.05,
                "nmix": 20,
                "mix": 4,
                "maxsteps": 400,
                "diag": "david",
                },
            "outdir": "calcdir_easy",
            }

        espresso_params_copy.update(easy_params)
        easy_calc = espresso(**espresso_params_copy)
        #__|

        #TODO: Find a way to freeze all atoms but adsorbates
        # SIMPLE RELAXATION #################
        print("ionic_opt | Running Easy Relaxation"); sys.stdout.flush()

        magmoms = atoms.get_initial_magnetic_moments()

        # set_mag_mom_to_0(atoms)

        # FIXME I should be able to just use the "set_initial_magnetic_moments"
        # method here. That should automatically set to magmoms to 0 from the
        # QE parameters dict
        atoms.set_initial_magnetic_moments(np.zeros(len(atoms)))

        atoms.set_calculator(easy_calc)

        qn = QuasiNewton(
            atoms,
            trajectory="out_opt_easy.traj",
            logfile="qn.log",
            )

        if os.path.exists("prev.traj"):
            qn.replay_trajectory("prev.traj")

        qn.run(fmax=fmax)

        # FULL RELAXATION #################
        print("ionic_opt | Running Full Relaxation"); sys.stdout.flush()
        # atoms.set_calculator(calc)
        atoms.set_calculator(calc_opt)

        print(magmoms)
        atoms.set_initial_magnetic_moments(magmoms)

        qn = QuasiNewton(
            atoms,
            trajectory="out_opt.traj",
            logfile="qn.log",
            )

        if os.path.exists("prev.traj"):
            qn.replay_trajectory("prev.traj")

        qn.run(fmax=fmax)
        #__|

    elif mode == "sp":
        #| - Single Point Calculation
        mess = "Running Single-Point Calculation"
        print(mess); sys.stdout.flush()

        atoms.get_potential_energy()
        write("out_opt.traj", atoms)
        #__|

    estimate_magmom(
        path_i=".",
        atoms=atoms,
        log="calcdir_opt/log",
        )

    elec_e = atoms.get_potential_energy()

    outdir = "dir_opt"
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    e_out = outdir + "/elec_e.out"
    with open(e_out, "w") as fle:
        fle.write(str(elec_e) + "\n")

    # if mode != "sp":
    #     #| - Always Run Single-Point Calculation with Full IO
    #     mess = "Running Post-run Single-Point Calculation"
    #     print(mess); sys.stdout.flush()
    #     atoms.get_potential_energy()
    #     #__|

    update_FINISHED("ionic_opt")

    if run_beef_an:
        an_beef_ensemble(atoms)

    if run_bader_an:

        #| - Running initial single-point calculation
        params_bader = {
            "output": {
                "avoidio": False,
                "removesave": True,
                "removewf": True,
                "wf_collect": False,
                },

            "outdir": "calcdir_bader",
            }

        calc_bands, espresso_params_bands = set_QE_calc_params(
            params=params_bader,
            )

        calc_bands = espresso(**espresso_params_bands)
        atoms.set_calculator(calc_bands)

        mess = "Running single-point calc with high io "
        print(mess); sys.stdout.flush()
        atoms.get_potential_energy()
        print("finished single-point"); sys.stdout.flush()
        #__|


        bader(atoms, spinpol=espresso_params_opt["spinpol"], run_exec=True)


#__| **************************************************************************

#__| **************************************************************************

#| - Magnetic Moments *********************************************************

def set_init_mag_moms(
    atoms,
    preference="bader",
    espresso_params=None,
    magmoms=None,
    read_from_file=False,
    ):
    """Set initial magnetic moments to atoms object using several methods.

    Set inital magnetic moments to atoms object. If the atoms object has
    previously had a bader or pdos analysis performed those magnetic moments
    will be available under the atoms.info dict ("pdos_magmoms" and
    "bader_magmoms" dict keys).

    # TODO Implement the "average" setting for the preference variable

    Args:
        atoms:
            Atoms object to be used
        preference:
            "bader"
            "pdos"
            "average"
        magmoms: <list> <dict>
            If specified as a list, set the atom's initial magnetic
            moments to list. If a element: magmom_i dict is given the magmoms
            will be initilized from this dict.
    """
    #| - set_init_mag_moms
    mess = "Setting Inital Magnetic Moments "
    mess += "**********************************************"
    print(mess); sys.stdout.flush()

    preference_map_dict = {
        "bader": "bader_magmoms",
        "pdos": "pdos_magmoms"
        }

    magmom_keys = ["pdos_magmoms", "bader_magmoms"]
    magmoms_master_dict = {}
    for magmom_key in magmom_keys:
        if magmom_key in atoms.info.keys():
            magmoms_master_dict[magmom_key] = atoms.info[magmom_key]

    magmom_keys = magmoms_master_dict.keys()

    preferred_an = preference_map_dict[preference]

    if magmoms is not None:

        if type(magmoms) == list:
            print("set_init_mag_moms | Using given magmoms"); sys.stdout.flush()
            magmoms_i = magmoms

        # TODO Implement this featrure in simple_mag_moms, providing dict will
        # update the init magmom dict
        elif type(magmoms) == dict:
            text = ("set_init_mag_moms | "
                    "Using simple method for initial magnetic moments "
                    "with updated values for dict - NOT IMPLEMENTED - 180425"
                    )
            print(text); sys.stdout.flush()
            magmoms_i = simple_mag_moms(atoms)

    else:
        spinpol_calc = calc_spinpol(atoms)

        if espresso_params is not None:
            if "spinpol" in list(espresso_params):
                spinpol_calc = espresso_params["spinpol"]

        #| - Spin-polarization turned off, set magmoms to 0
        if not spinpol_calc:
            print("set_init_mag_moms | Spin-polarization turned off")
            sys.stdout.flush()
            mag_mom_list = atoms.get_initial_magnetic_moments()
            magmoms_i = np.zeros(len(mag_mom_list))
        #__|

        # COMBAK This is a poor way of enforcing the analysis method preference
        # Assumes that there are only 2 analysis methods
        elif preferred_an in magmoms_master_dict.keys():
            text = ("set_init_mag_moms | "
                    "Using preferred method for initial magnetic moments "
                    "(" + preferred_an + ")"
                    )
            print(text); sys.stdout.flush()
            magmoms_i = magmoms_master_dict[preferred_an]

        elif len(magmoms_master_dict.keys()) == 1:
            an_method = magmoms_master_dict.keys()[0]
            text = ("set_init_mag_moms | "
                    "Using " + an_method + " for initial magnetic moments")
            print(text); sys.stdout.flush()

            magmoms_i = magmoms_master_dict[an_method]

        #| - Use simple method
        else:
            text = ("set_init_mag_moms | "
                    "Using simple method for initial magnetic moments")
            print(text); sys.stdout.flush()
            magmoms_i = simple_mag_moms(atoms)
        #__|

        if read_from_file:
            text = ("set_init_mag_moms | "
                    "Using magmom_init.in file for initial magnetic moments")
            print(text); sys.stdout.flush()

            magmoms_tmp = read_magmoms_from_file()
            if magmoms_tmp is not None:
                magmoms_i = magmoms_tmp


    #| - Check That Length of Magmoms_list == len(atoms)
    if not len(magmoms_i) == len(atoms):
        text = ("Length of magmoms doesn't match the number of atoms!!!!!!!!!!!"
                "\n Will use simple method to assign initial magmoms")

        print(text); sys.stdout.flush()
        magmoms_i = simple_mag_moms(atoms)
    #__|


    magmoms_i_new = increase_abs_val_magmoms(atoms, magmoms_i)

    atoms.set_initial_magnetic_moments(magmoms_i_new)

    #| - Printing Magnetic Moments
    print("set_init_mag_moms | Initial Magnetic Moments:")
    for atom in atoms:
        elem = atom.symbol
        magmom = atom.magmom
        print(elem + ": " + str(magmom))

    #__|

    reduce_magmoms(atoms)
    #__|

def increase_abs_val_magmoms(atoms, magmoms_list, increase_amount=0.8):
    """Increase absolute value of magmoms for atoms object.

    # COMBAK Shouldn't raise initial guess for light atoms (Or don't raise as
    much)

    Will not modify magmoms == 0.

    Select light atoms will have their |magmom| raised by a set small amount to
    not oversaturate the atomic magnetic moment (init_magmom > # valence e)

    Args:
        magmoms_list:
        increase_amount:
    """
    #| - increase_abs_val_magmoms
    inc = increase_amount

    light_atoms_list = {
        "H": 0.1,
        "He": 0.3,
        "Li": 0.4,
        "Be": 0.5,
        "B": 0.6,
        "C": 0.7,
        "N": 0.8,
        "O": 0.9,
        "F": 1.0,
        }

    new_magmom_list = []
    # for magmom in magmoms_list:
    for atom, magmom in zip(atoms, magmoms_list):

        elem_i = atom.symbol

        if elem_i in light_atoms_list.keys():
            inc = light_atoms_list[atom.symbol]

        if magmom < 0.:
            new_magmom = abs(magmom) + inc
            new_magmom = -new_magmom

        elif magmom > 0.:
            new_magmom = abs(magmom) + inc

        else:
            new_magmom = 0.

        # Making sure H initial magmom isn't > 1
        if elem_i == "H":
            if new_magmom > 0.99:
                print("Setting H magmom to 1 - RF - 180423 - TEMP")
                new_magmom = 1.

        new_magmom_list.append(new_magmom)


    return(new_magmom_list)
    #__|

def calc_spinpol(atoms):
    """Return whether spin polarization should be turned on or off.

    The atoms object must have the appropriate calculator object attachedd with
    the relevent parameters declared

    Args:
        atoms:
    """
    #| - calc_spinpol
    spinpol = True
    if hasattr(atoms, "calc"):
        if atoms.calc is not None:
            spin_key = "spinpol"
            params_dict = atoms.calc.__dict__
            if spin_key in params_dict:
                spin_pol = params_dict[spin_key]
                if spin_pol is False:
                    spinpol = False
                    # set_mag_mom_to_0(atoms)
                    # return(spinpol)

    return(spinpol)
    #__|

def simple_mag_moms(atoms):
    """Implement simple procedure to set initial guess for magnetic moments.

    Args:
        atoms
    """
    #| - simple_mag_moms

    #| - High Spin Dictionary
    master_dict_high_spin = {

        "Cr": 5, "Mn": 5,
        "Cu": 1, "Ag": 0, "Au": 0,  # 3d10
        "Ni": 3, "Pd": 3, "Pt": 3,  # 3d8, 4s1
        "Co": 4, "Rh": 4, "Ir": 4,  # 3d7, 4s1
        "Fe": 5, "Ru": 5, "Os": 5,  # 3d6, 4s1
        "Mo": 5, "Tc": 5, "Hf": 2,
        "Ta": 3, "W": 4, "Re": 5, "Ti": 2,
        }

    spin_dict = {
        "O": 2.5,
        "H": 0.2,
        "C": 0.2,
        "N": 0.2,
        }
    #__|

    #| - Find cations
    cations = []
    for atom in atoms:
        if atom.symbol in master_dict_high_spin.keys():
            # if master_dict_high_spin.has_key(atom.symbol):
            cations.append(atom.symbol)
    #__|

    #| - Cation Magnetic Moments Dictionary
    cation_magmom_dict = {}
    for cation in cations:
        cation_magmom_dict[cation] = master_dict_high_spin[cation] * 1.2 + 0.5
    #__|

    #| - Setting Magnetic Moments of Atoms
    magmoms = atoms.get_initial_magnetic_moments()

    for atom in atoms:
        if atom.symbol in cations:
            magmoms[atom.index] = cation_magmom_dict[atom.symbol]
            continue

        for elem_i in spin_dict.keys():
            if atom.symbol == elem_i:
                magmoms[atom.index] = spin_dict[elem_i]
                break

            else:
                magmoms[atom.index] = 0.0
    magmoms = np.array(magmoms)

    return(magmoms)
    #__|

    #__|

def reduce_magmoms(atoms, ntypx=10):
    """Reduce number of unique magnetic moments of atoms object.

    Reduce the number of unique magnetic moments by combining those that are
    most similar among atoms with the same atomic symbol. This is necessary for
    atoms objects with more than 10 types of magmom/symbol pairs because QE only
    accepts a maximum of 10 types of atoms.
    """
    #| - reduce_magmoms
    syms = set(atoms.get_chemical_symbols())

    master_dict = {}
    for sym in syms:
        master_dict[sym] = {}

    for atom in atoms:
        if atom.magmom in master_dict[atom.symbol]:
            master_dict[atom.symbol][atom.magmom].append(atom.index)
        else:
            master_dict[atom.symbol][atom.magmom] = [atom.index]

    ntyp = 0
    for sym in syms:
        ntyp += len(master_dict[sym].keys())

    while ntyp > ntypx:
        magmom_pairs = {}
        for sym in syms:
            magmoms = master_dict[sym].keys()
            if not len(magmoms) > 1: continue
            min_delta = 1e6
            min_pair = ()
            for i, magmom1 in enumerate(magmoms):
                for j, magmom2 in enumerate(magmoms):
                    if not i < j: continue
                    delta = np.abs(magmom1 - magmom2)
                    if delta < min_delta:
                        min_delta = delta
                        min_pair = (magmom1, magmom2)

            assert min_delta != 1e6
            assert min_pair != ()
            magmom_pairs[sym] = min_pair

        min_delta = 1e6
        min_sym = ""
        for sym in magmom_pairs:
            delta = np.abs(magmom_pairs[sym][0] - magmom_pairs[sym][1])
            if delta < min_delta:
                min_delta = delta
                min_sym = sym

        assert min_delta != 1e6
        assert min_sym != ""
        if min_delta > 0.5:
            warn = "WARNING, reducing pair of magmoms whose difference is "
            print(warn + "%.2f" % min_delta)

        if np.abs(magmom_pairs[min_sym][0]) > np.abs(magmom_pairs[min_sym][1]):
            master_dict[min_sym][magmom_pairs[min_sym][0]].extend(
                master_dict[min_sym][magmom_pairs[min_sym][1]])
            del master_dict[min_sym][magmom_pairs[min_sym][1]]
        else:
            master_dict[min_sym][magmom_pairs[min_sym][1]].extend(
                master_dict[min_sym][magmom_pairs[min_sym][0]])
            del master_dict[min_sym][magmom_pairs[min_sym][0]]

        ntyp = 0
        for sym in syms:
            ntyp += len(master_dict[sym].keys())

    #reassign magmoms
    for sym in syms:
        for magmom in master_dict[sym]:
            for index in master_dict[sym][magmom]:
                atoms[index].magmom = magmom
    #__|

def read_magmoms_from_file(file_name="magmom_init.in"):
    """Read inital magmoms from a file.

    Args:
        file_name:
    """
    #| - read_magmoms_from_file
    print(os.path.isfile(file_name))
    magmoms = None
    if os.path.isfile(file_name):
        try:
            with open(file_name, "r") as fle:
                magmoms = [float(i.strip()) for i in fle.readlines()]
        except:
            print("Couldn't read init magmom file")

    return(magmoms)
    #__|

# def compare_magmoms(self):
def compare_magmoms():
    """Compare spin states of two atoms objects.

    (I think I got this from Colin on 180413)

    Part of bigger script
    /home/colinfd/usr/bin/get_G.py

    TODO Port this code to my workflow

    Author: Colin Dickens
    """
    #| - compare_magmoms
    # def nearest_atom(atoms, position):
    #     """Returns atom nearest to position."""
    #     #| - nearest_atom
    #     position = np.array(position)
    #     dist_list = []
    #     for atom in atoms:
    #         dist = np.linalg.norm(position - atom.position)
    #         dist_list.append(dist)
    #
    #     return atoms[np.argmin(dist_list)]
    #     #__|
    #
    # if len(self.ads_atoms) >= len(self.slab_atoms):
    #     ads = self.ads_atoms
    #     slab = self.slab_atoms
    #     indexed_by = self.slab
    #     not_indexed_by = self.ads
    # else:
    #     slab = self.ads_atoms
    #     ads = self.slab_atoms
    #     indexed_by = self.ads
    #     not_indexed_by = self.slab
    #
    # delta_magmoms = []
    # ads_indices_used = []
    # for atom in slab:
    #     ads_atom = nearest_atom(ads, atom.position)
    #     if not self.quiet:
    #         if ads_atom.symbol != atom.symbol:
    #             print("WARNING! MAGMOM COMPARISON FAILURE")
    #     ads_indices_used.append(ads_atom.index)
    #     delta_magmoms.append(atom.magmom - ads_atom.magmom)
    #
    # ads_indices_not_used = []
    # for i in range(len(ads)):
    #     if i not in ads_indices_used:
    #         ads_indices_not_used.append(i)
    #
    # self.delta_magmoms = zip(range(len(slab)), delta_magmoms)
    # self.delta_magmoms.sort(key=lambda x: abs(x[1]), reverse=True)
    #
    # common = ""
    # uncommon = ""
    # for i in range(8):
    #     atom = slab[self.delta_magmoms[i][0]]
    #     common += "%s%d: %.2f\t" % (atom.symbol,
    #         atom.index,
    #         self.delta_magmoms[i][1],
    #         )
    #
    # for i in ads_indices_not_used:
    #     uncommon += "%s%d: %.2f\t"%(ads[i].symbol,ads[i].index,ads[i].magmom)
    #
    # if self.quiet:
    #     return
    # print("~" * 6 + "MAGNETIC MOMENT COMPARISON" + "~" * 6)
    # print("Largest magmom discrepancies (indexed by %s)" % indexed_by)
    # print(common)
    # print("Magnetic moments only present in %s" % not_indexed_by)
    # print(uncommon + "\n")
    #__|

#__| **************************************************************************

#| - Density of States ********************************************************
def an_pdos(
    atoms,
    # calc,
    dos_kpts,
    espresso_params,
    ):
    """Perform projected density of states (PDOS) analysis.

    # TODO Clean up, makes too much data

    Args:
        atoms:
        dos_kpts:
        espresso_params:
    """
    #| - an_pdos
    from espresso import espresso

    mess = "Running PDOS Analysis "
    mess += "********************************************************"
    print(mess); sys.stdout.flush()

    outdir = "dir_pdos"

    params_dos = {
        "output": {
            "avoidio": False,
            "removesave": True,
            "removewf": False,
            "wf_collect": True,
            },

        "outdir": "calcdir_dos",
        }

    calc_opt, espresso_params_dos = set_QE_calc_params(
        params=params_dos,
        )

    calc_dos = espresso(**espresso_params_dos)
    atoms.set_calculator(calc_dos)

    reduce_magmoms(atoms)

    mess = "Running single-point calc with high io "
    print(mess); sys.stdout.flush()
    atoms.get_potential_energy()
    print("finished single-point"); sys.stdout.flush()

    cwd = os.getcwd()
    dos = atoms.calc.calc_pdos(
        nscf=True,
        kpts=dos_kpts,
        Emin=-20.0,
        Emax=20.0,
        ngauss=0,
        sigma=0.2,
        DeltaE=0.01,
        tetrahedra=False,
        slab=True,
        )
    os.chdir(cwd)

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    pdos_out = outdir + "/dos.pickle"
    with open(pdos_out, "w") as fle:
        pickle.dump(dos, fle)

    # Set Magnetic Moments To Atoms Object From PDOS Intergration
    spin_pdos(
        atoms,
        pdos_pkl=pdos_out,
        outdir=outdir,
        spinpol=espresso_params["spinpol"],
        )

    atoms.write(outdir + "/out_pdos.traj")
    update_FINISHED("an_pdos")
    #__|

def spin_pdos(
    atoms,
    pdos_pkl=None,
    valence_dict=None,
    nscf=False,
    kpts=None,
    outdir=None,
    save_pkl=False,
    spinpol=False,
    write_charges=True,
    **kwargs
    ):
    """Calculate spin moments on each atom from PDOS analysis.

    Calculate more accurate charges/magnetic moments from pdos and assign to
    atom.charge/atom.magmom. If pdos_pkl not defined, it will be calculated
    (atoms object must have a real calculator attached). If pdos_pkl is
    defined, will attempt to load from pickle file, but valence_dict must be
    specified! Specify calculation directory as outdir for easy cleanup.

    Args:
        atoms:
        pdos_pkl:
        valence_dict:
        nscf:
        kpts:
        outdir:
        save_pkl:
        spinpol:
        write_charges:
        kwargs:
    """
    #| - spin_pdos
    valence_dict = {
        "Cu": 11, "C": 4, "O": 6, "H": 1, "Li": 1,
        "Rh": 17, "Co": 9, "Pd": 10, "Pt": 10,
        "Ni": 1, "Fe": 16, "N": 5, "Ru": 16,
        }

    #| - Reading PDOS File, Otherwise Creates It
    # FIXME I don't like that it can create PDOS file, this is handled by my
    # an_pdos method.

    if pdos_pkl:
        print("spin_pdos | pdos_pkl is not None")  # TEMP For testing purposes
        assert valence_dict is not None, "MUST SPECIFY valence_dict"
        single_point_calc = True
        pdos = pickle.load(open(pdos_pkl))
        nvalence_dict = valence_dict

    else:
        print("spin_pdos | pdos_pkl is None")  # TEMP For testing purposes
        single_point_calc = False
        # dict with (chemical symbol) : (num valence)
        nvalence_dict = atoms.calc.get_nvalence()[1]
        # double k-points for higher quality pdos --> O(single point calc)
        if nscf:
            print("spin_pdos | TEMP1")
            pdos = atoms.calc.calc_pdos(nscf=True, kpts=kpts, **kwargs)
        else:  # no single point calc, should take 1 or 2 minutes
            print("spin_pdos | TEMP2")
            pdos = atoms.calc.calc_pdos(**kwargs)
    #__|

    #| - Finding Index of Fermi Level
    for i_ind, e_i in enumerate(pdos[0]):
        if e_i > 0:
            fi = i_ind  # index of fermi level
            break
    #__|

    #| - Analysing PDOS For Magnetic Moments and Charge of All Atoms
    if spinpol:
        #| - Spin Polarlized Calculation
        magmom_list = []
        charge_list = []
        for i, atom in enumerate(atoms):

            #| - Integrating Up and Down Spin PDOS
            spin_up = 0; spin_down = 0
            for sym in pdos[2][i]:
                spin_up += np.trapz(pdos[2][i][sym][0][:fi], x=pdos[0][:fi])
                spin_down += np.trapz(pdos[2][i][sym][1][:fi], x=pdos[0][:fi])

            #__|

            #| - Update Atoms Magmom and Charge Properties
            ##Update magmom
            if np.abs(spin_up - spin_down) > 1e-4:
                magmom_i = spin_up - spin_down
            else:
                magmom_i = 0.

            magmom_list.append(magmom_i)
            atom.magmom = magmom_i

            ##Update charge
            charge_i = nvalence_dict[atom.symbol] - (spin_up + spin_down)
            if write_charges:
                # atom.charge = nvalence_dict[atom.symbol]-(spin_up+spin_down)
                atom.charge = charge_i

            charge_list.append(charge_i)
            #__|

        print("PDOS MAGMOMS: " + str(atoms.get_initial_magnetic_moments()))
        reduce_magmoms(atoms)

        atoms.info.update({"magmom_set": True})
        atoms.info.update({"pdos_magmoms": magmom_list})
        atoms.info.update({"pdos_charges": charge_list})

        pickle.dump(magmom_list, open("%s/magmom_list.pickle" % outdir, "w"))
        pickle.dump(charge_list, open("%s/charge_list.pickle" % outdir, "w"))
        #__|

    else:
        #| - Non-Spin Polarized Calculation
        charge_list = []
        for i, atom in enumerate(atoms):

            #| - Integrating PDOS For Charges
            charge = 0
            for sym in pdos[2][i]:
                charge += np.trapz(pdos[2][i][sym][0][:fi], x=pdos[0][:fi])
            #__|

            #| - Update Atoms Charges
            # Update charge
            charge_i = nvalence_dict[atom.symbol] - (charge)
            if write_charges:
                # atom.charge = nvalence_dict[atom.symbol] - (charge)
                atom.charge = charge_i
            charge_list.append(charge_i)
            #__|

            atoms.info.update({"pdos_charges": charge_list})

            pickle.dump(
                charge_list,
                open(
                    "%s/charge_list.pickle" % outdir,
                    "w",
                    )
                )
        #__|

    print("PDOS CHARGES: " + str(atoms.get_initial_charges()))
    #__|

    #| - Writing Output To File
    if outdir and not single_point_calc:
        #Save only the Lowdin charges in the pdos log file
        light_lines = []

        f = open("%s/pdos.log" % outdir)
        lines = f.readlines()
        f.close()

        skip = False
        for line in lines:
            try:
                if line.split()[0] == "k":
                    skip = True
                elif line.split()[0] == "Lowdin":
                    skip = False
            except: continue

            if not skip:
                light_lines.append(line)

        f = open("%s/pdos_Lowdin.log" % outdir, "w")
        f.writelines(light_lines)
        f.close()
        os.system("rm %s/pdos.log" % outdir)

    if save_pkl:
        pickle.dump(pdos, open("pdos.pkl", "w"))
    #__|

    #__|

#__| **************************************************************************

#| - Band Structure ***********************************************************

def an_bands(atoms, bands_kpts, espresso_params):
    """Perform band analysis on atoms object.

    Development Notes:
        * TODO Creates too much data, clean up
        * TODO Check if bands have been calculated earlier

    Args:
        atoms:
        bands_kpts:
        espresso_params:
    """
    #| - an_bands
    from espresso import espresso

    mess = "Executing Band Structure Analysis "
    mess += "********************************************"
    print(mess); sys.stdout.flush()

    # FIXME This should be handled by envoking the
    # set_initial_magnetic_moments method
    spinpol_calc = calc_spinpol(atoms)
    if spinpol_calc is False:
        print("an_bands | Spin-polarization turned off, setting magmoms to 0")
        atoms.set_initial_magnetic_moments(np.zeros(len(atoms)))
        # set_mag_mom_to_0(atoms)  # COMBAK Remove this if working

    #| - Running initial single-point calculation
    params_bands = {
        "output": {
            "avoidio": False,
            "removesave": False,
            "removewf": True,
            "wf_collect": False,
            },

        "kpts": bands_kpts,
        # "outdir": "dir_bands"

        "outdir": "calcdir_bands",
        }

    calc_bands, espresso_params_bands = set_QE_calc_params(
        params=params_bands,
        )

    calc_bands = espresso(**espresso_params_bands)
    atoms.set_calculator(calc_bands)

    mess = "Running single-point calc with high io "
    print(mess); sys.stdout.flush()
    atoms.get_potential_energy()
    print("finished single-point"); sys.stdout.flush()
    #__|

    # # calc_spinpol(atoms)  # COMBAK I don't think this is doing anything
    # espresso_params.update(
    #     {
    #         "kpts": bands_kpts,
    #         "outdir": "dir_bands"
    #         },
    #     )
    # atoms.calc = espresso(**espresso_params)
    # mess = "Running single-point calculation"
    # print(mess); sys.stdout.flush()
    # atoms.get_potential_energy()

    atoms.calc.save_flev_chg("charge_den.tgz")
    atoms.calc.load_flev_chg("charge_den.tgz")

    # COMBAK How to handle this more generally?
    ip = ibz_points["fcc"]
    points = ["Gamma", "X", "W", "K", "L", "Gamma"]
    bzpath = [ip[p] for p in points]
    kpts, x, X = get_bandpath(bzpath, atoms.cell, npoints=300)

    mess = "Calculating band structure"
    print(mess); sys.stdout.flush()

    # Note Need current dir because calc_bandstructure moves to calc_dir
    cwd = os.getcwd()
    energies = atoms.calc.calc_bandstructure(kpts, atomic_projections=True)
    os.chdir(cwd)

    if not os.path.exists("dir_bands"):
        os.makedirs("dir_bands")

    with open("dir_bands/band_disp.pickle", "w") as fle:
        pickle.dump((points, kpts, x, X, energies), fle)

    # COMBAK This broke when file already existed, tried a fix - 180405 - RF
    shutil.move("charge_den.tgz", "dir_bands/charge_den.tgz")
    atoms.write("dir_bands/out_bands.traj")

    update_FINISHED("an_bands")
    #__|

#__| **************************************************************************

#| - Beef Ensemble of Energies ************************************************
# def an_beef_ensemble(atoms, xc):  # COMBAK
def an_beef_ensemble(atoms):
    """Perform BEEF ensemble of enery analysis.

    Atoms must have an initialized calculator object attached.

    FIXME The xc parameter is no longer needed it can be gleamed from the
    atoms' calculator object

    Args:
        atoms:
        xc:
    """
    #| - an_beef_ensemble
    mess = "Executing BEEF Ensemble Analysis "
    mess += "*********************************************"
    print(mess); sys.stdout.flush()

    beefensemble_on = atoms.calc.beefensemble
    xc = atoms.calc.xc

    if beefensemble_on is False:
        print("The ase-espresso calculator has the beefensemble"
            " parameter turned off!, analysis will not run!!")
        return(None)

    allowed_xc = ["BEEF", "BEEF-vdW"]
    # if xc != "BEEF" or xc != "BEEF-vdW":
    if xc not in allowed_xc:
        print("The exchange-correlation functional has to be either "
            "BEEF or BEEF-vdW to do the BEEF ensemble analysis")
        return(None)

    # if xc == "BEEF" or xc == "BEEF-vdW":

    # calc = atoms.calc
    from ase.dft.bee import BEEFEnsemble

    energy = atoms.get_potential_energy()
    # ens = BEEFEnsemble(atoms=atoms, e=energy, xc="beefvdw")
    # NOTE The BEEFEnsemble class should parse xc from atoms object
    # COMBAK
    ens = BEEFEnsemble(atoms=atoms, e=energy)
    ens_e = ens.get_ensemble_energies()

    if not os.path.exists("dir_beef_ensemble"):
        os.makedirs("dir_beef_ensemble")

    with open("dir_beef_ensemble/ensemble.pickle", "w") as fle:
        pickle.dump(ens_e, fle)

    update_FINISHED("an_beef_ensemble")
    #__|

def plot_beef_ensemble(
    folder_dir="dir_beef_ensemble",
    file_name="ensemble.pickle",
    file_out="beef_ens_hist.png",
    ):
    """Create plot of distribution of energies from BEEF ensemble.

    Args:
        folder_dir:
        file_name:
        file_out:
    """
    #| - plot_beef_ensemble
    file_loc = folder_dir + "/" + file_name

    data = pickle.load(open(file_loc, "r"))

    mu, std = norm.fit(data)
    plt.hist(data, bins=25, normed=False, alpha=0.6, color="g")

    xmin, xmax = plt.xlim()
    x = np.linspace(xmin, xmax, 100)
    p = norm.pdf(x, mu, std)
    plt.plot(x, p, 'k', linewidth=2)
    title = "Fit results: mu = %.2f,  std = %.2f" % (mu, std)
    plt.title(title)

    out_file = folder_dir + "/" + file_out
    plt.savefig(out_file)
    # plt.show()
    #__|

#__| **************************************************************************

#| - Vibrational Analysis *****************************************************
def an_ads_vib(
    atoms,
    espresso_params=None,
    ads_index_list=None,
    # thermochem_corrections=True,
    thermochem_corrections="harmonic",  # "harmonic" or "IG"
    remove_imag_modes=True,

    Temperature=300.0,
    Pressure=100000.0,
    symmetrynumber=2,
    spin=0,
    linear=True,
    ):
    """Adsorbate vibrational analysis.

    Returns the ASE Vibrations instance in case further analysis is needed.

    Args:
        atoms:
        ads_index_list:
        thermochem_corrections:
        remove_imag_modes: Removes imaginary modes
    """
    #| - an_ads_vib
    # from espresso import espresso
    from espresso.vibespresso import vibespresso

    mess = "Starting vibrational analysis "
    mess += "************************************************"
    print(mess); sys.stdout.flush()

    #| - Setting Adsorbate Index List
    if ads_index_list is not None:
        pass
    elif "adsorbates" in atoms.info.keys():
        print("Adsorbate info present! Good good.")
        ads_index_list = atoms.info["adsorbates"]
    else:
        print("an_ads_vib | Adsorbate index info couldn't be parsed from atoms")
        print("an_ads_vib | Will vibrate all atoms!!!!!!!! (Probably not good)")
        pass
    #__|

    #| - Removing Empty Pickle Files
    pckl_file_list = glob.glob("dir_vib/*.pckl*") + glob.glob("*.pckl*")
    for pckl_file in pckl_file_list:
        if os.stat(pckl_file).st_size == 0:
            os.remove(pckl_file)
            print("an_ads_vib | " + pckl_file + " empty, so removed")
    #__|

    #| - Copy vib.pckl files back to root dir (For restarting)
    for fle in glob.glob("dir_vib/*.pckl*"):
        fle_name = fle.split("/")[-1]
        shutil.move(fle, fle_name)
    #__|

    #| - Running initial single-point calculation
    params_vib = {
        "output": {
            "avoidio": False,
            "removesave": False,
            "removewf": False,
            "wf_collect": True,
            },

        # "kpts": bands_kpts,
        # "outdir": "dir_bands"

        "outdir": "calcdir_vib",
        }

    calc_vib, espresso_params_vib = set_QE_calc_params(
        params=params_vib,
        )

    calc_vib = vibespresso(
        outdirprefix="out_vib",
        **espresso_params_vib)

    # atoms.set_calculator(calc)
    # calc_vib = espresso(**espresso_params_vib)

    atoms.set_calculator(calc_vib)

    # mess = "Running single-point calc with high io "
    # print(mess); sys.stdout.flush()
    # atoms.get_potential_energy()
    # print("finished single-point"); sys.stdout.flush()
    #__|

    set_init_mag_moms(
        atoms,
        preference="bader",
        espresso_params=espresso_params_vib,
        )

    # set_init_mag_moms(atoms)

    vib = Vibrations(atoms, indices=ads_index_list)
    vib.run()

    # print("an_ads_vib | Getting vibrational energies")
    vib_e_list = vib.get_energies()

    vib.summary(log="vib_summ.out")

    #| - Copy Files to dir_vib Folder
    if not os.path.exists("dir_vib"):
        os.makedirs("dir_vib")
    shutil.move("vib_summ.out", "dir_vib/vib_summ.out")

    dest_dir = "dir_vib"
    for fle in glob.glob(r'*.pckl*'):
        shutil.move(fle, dest_dir + "/" + fle)


    for fle in glob.glob(r'*out_vib*'):
        shutil.move(fle, dest_dir + "/" + fle)
    #__|

    #| - Remove Imaginary Modes
    # Removes them, not just making them 0
    if remove_imag_modes:
        vib_e_list = vib_e_list.real
        vib_e_list = [vib_i for vib_i in vib_e_list if vib_i != 0.]
    #__|

    if thermochem_corrections == "harmonic":
        thermochem_harm_corr(vib_e_list)

    elif thermochem_corrections == "IG":
        thermochem_IG_corr(
            vib_e_list,
            Temperature=Temperature,
            Pressure=Pressure,
            potentialenergy=0.,
            symmetrynumber=symmetrynumber,
            spin=spin,
            atoms=atoms,
            linear=linear,
            )

    # Saving Vibrational Mode List
    with open("dir_vib/vib_modes.pickle", "w") as fle:
        pickle.dump(vib_e_list, fle)

    #| - Saving Vibrations Class Instance
    file_name = "dir_vib/vib_inst.pickle"
    if not os.path.exists(file_name):
        pass
    else:
        os.remove(file_name)

    # Pickling error occured when pickling the atoms object
    if "atoms" in vars(vib): del vars(vib)["atoms"]

    with open(file_name, "w") as fle:
        pickle.dump(vib, fle)
    #__|

    update_FINISHED("an_ads_vib")

    # return(vib_e_list)
    return(vib)
    #__|

def thermochem_harm_corr(
    vib_e_list,
    Temperature=300.0,
    ):
    """Thermochemical free energy corrections from vibrational analysis.

    Args:
        vib_e_list:
            List of vibrational modes in eV
        Temperature:
    """
    #| - thermochem_harm_corr
    mess = "Starting thermochemical harmonic energy contributions "
    mess += "***********************"
    print(mess); sys.stdout.flush()

    print("Calculating thermochemical corrections @ " + str(Temperature) + "K")

    vib_e_list = np.array(vib_e_list)

    # Remove imaginary frequencies
    vib_e_list = vib_e_list.real
    vib_e_list = [vib_i for vib_i in vib_e_list if vib_i != 0.]

    ht = HarmonicThermo(vib_e_list)
    F_energy = ht.get_helmholtz_energy(Temperature)

    if not os.path.exists("dir_vib"):
        os.makedirs("dir_vib")

    with open("dir_vib/gibbs_corr.out", "w") as fle:
        fle.write(str(F_energy))
        fle.write("\n")
    #__|

def thermochem_IG_corr(
    vib_energies,
    Temperature=300.0,
    Pressure=100000.0,
    potentialenergy=0.,
    symmetrynumber=None,
    spin=None,
    atoms=None,
    linear=True,
    ):
    """Calculate free energy to gas phase molecule from vibrational modes.

    Args:
        vib_e_list:
        Temperature=300.0:
        Pressure=100000.0:
        potentialenergy=0.:
        symmetrynumber=None:
        spin=None:
        atoms=None:
    """
    #| - thermochem_IG_corr
    mess = "Starting thermochemical ideal gas energy contributions "
    mess += "***********************"
    print(mess); sys.stdout.flush()

    from ase.thermochemistry import IdealGasThermo

    # try:
    #     with open("dir_vib/vib_inst.pickle", "r") as fle:
    #         vib = pickle.load(fle)
    # except:
    #     pass
    #
    # vib_energies = vib.get_energies()

    if atoms is not None:
        potentialenergy = atoms.get_potential_energy()

    if linear is True:
        lin = "linear"
    else:
        lin = "nonlinear"

    thermo = IdealGasThermo(
        vib_energies,
        lin,  # linear/nonlinear
        potentialenergy=potentialenergy,
        atoms=atoms,

        # H2O, H2, O2: 2,
        symmetrynumber=symmetrynumber,  # symmetry numbers from point group
        spin=spin,  # 1 for O2, 0 for H2O and H2
        )

    G = thermo.get_gibbs_energy(temperature=Temperature, pressure=Pressure)

    if not os.path.exists("dir_vib"):
        os.makedirs("dir_vib")

    with open("dir_vib/g_energy.out", "w") as fle:
        fle.write(str(G))
    #__|


#__| **************************************************************************

#| - Job Cleanup
def clean_up_dft():
    """Clean up files after DFT job.

    I'm using a lot of try statement because I dont want this to break anywhere
    """
    #| - clean_up_dft
    with open(".FINISHED", "w") as fle:
        fle.write("\n")

    # def update_FINISHED(text, filename=".FINISHED.new"):

    update_FINISHED("job_completed")

    if not os.path.exists("__misc__"):
        os.makedirs("__misc__")

    #| - Moving Sherlock Node Files to __misc__ folder
    # Sherlock didn't like this:
    # Open RTE was unable to open the hostfile:
    #     /scratch/users/flores12/03_graph_N_Fe/01_opt_struct/N_doped_graph_Fe/1-att/__test__/1-att/_5/uniqnodefile.21286552
    # Check to make sure the path and filename are correct.


    # try:
    #     if "COMPENV" in os.environ:
    #         compenv = os.environ["COMPENV"]
    #         if compenv == "sherlock":
    #             os.system("mv *nodefile.* __misc__")
    #             os.system("mv *uniqnodefile..* __misc__")
    #
    # except:
    #     pass
    #__|

    #| - Moving DFT Parameter Files to dir_dft_params
    try:
        if not os.path.exists("dir_dft_params"):
            os.makedirs("dir_dft_params")

        # os.system("ls")
        # os.system("pwd")
        os.system("mv *dft-params* dir_dft_params")

    except:
        pass
    #__|

    #__|

#__|

#| - Atoms File Operations ****************************************************
def read_atoms_from_file(filename=None, try_restart=True):
    """Read atoms object from file.

    Checks several file names

    TODO Options and methods to read and check output trajectory files

    Args:
        filename: optional atoms file, will attempt to read first.
        try_restart: Restart job by reading output atoms/trajectory object
    """
    #| - read_atoms_from_file
    mess = "Reading Atoms Object From File "
    mess += "***********************************************"
    print(mess)

    atoms = None
    traj = None

    #| - Restart From Output Atoms Object
    restart_file_name_list = [
        "out.traj",
        "out_opt.traj",
        "out_final.traj",
        ]

    if try_restart:
        for file_name in restart_file_name_list:
            if os.path.isfile(file_name):

                try:
                    atoms = read(file_name)
                    traj = Trajectory(file_name, "a", atoms)

                    mess = "read_atoms_from_file | Restarting from "
                    mess += file_name
                    print(mess); sys.stdout.flush()

                    break

                except:
                    pass
    #__|

    #| - Starting From Fresh Atoms Object
    file_name_list = [
        "init.traj",
        "init.POSCAR",
        "manual_start.traj",
        "POSCAR",
        ]

    if atoms is None:

        if filename is not None:
            file_name_list.insert(0, filename)

        for file_name in file_name_list:
            if os.path.isfile(file_name):
                atoms = read(file_name)
                # COMBAK Test feature
                traj = Trajectory("out_opt.traj", "w", atoms)

                mess = "read_atoms_from_file | Starting from "
                mess += file_name
                print(mess); sys.stdout.flush()

                break
    #__|

    if atoms is None:
        raise IOError("No atoms file found")

    return(atoms, traj)
    #__|

def convert_atoms_object(atoms_filename, out_file):
    """Convert atoms objects to new file format.

    Args:
        atoms_filename:
        out_file:
    """
    #| - convert_atoms_object
    atoms = read(atoms_filename)

    #| - Convert to New Trajectory File Format
    if out_file.split(".")[-1] == "traj":
        write(out_file, atoms)

        #| - Using Trajectory Class Directly (Not Right?)
        # traj = Trajectory(out_file, mode="w", atoms=atoms)
        # traj.set_description({"author": "Created by Raul Flores"})
        # traj.write()
        #__|

    #__|

    #__|

#__| **************************************************************************

#| - Atoms Geometry Methods ***************************************************
def angle_between_lattice_vectors(atoms, vector_0=0, vector_1=1):
    """Calculate angle between cell lattice vectors.

    Calculates the angle between the first 2 lattice vectors of a computational
    cell in degrees.

    Args:
        atoms:
        vector_0:
        vector_1:
    """
    #| - angle_between_lattice_vectors
    v1 = atoms.cell[vector_0]
    v2 = atoms.cell[vector_1]

    angle = angle_between(v1, v2)
    angle = math.degrees(angle)

    return(angle)
    #__|

def magnitude_of_lattice_vectors(atoms):
    """Return magnitude of three lattice vectors.

    Args:
        atoms:
    """
    #| - magnitude_of_lattice_vectors
    v1 = atoms.cell[0]
    v2 = atoms.cell[1]
    v3 = atoms.cell[2]

    mag1 = np.linalg.norm(v1)
    mag2 = np.linalg.norm(v2)
    mag3 = np.linalg.norm(v3)

    out_tup = (mag1, mag2, mag3)

    return(out_tup)
    #__|

def find_diff_between_atoms_objects(atoms_A, atoms_B):
    """Find indices of atoms that are unique to atoms_A and atoms_B.

    Given two atoms objects (atoms_A and atoms_B), finds the atoms that are in
    atoms_B but not in atoms_A and conversly the atoms that are in atoms_A but
    not in atoms_B

    The criteria for equivalence between atoms is that their positions and
    element type are the same

    This method should only be used for cases where atoms are added or removed,
    if a relaxation is performed then the methods will fail to find the union
    of atoms between both atoms objects since there positions will no longer be
    exact

    Args:
        atoms_A:
        atoms_B:
    """
    #| - find_diff_between_atoms_objects

    #| - __old__

    # #| - Import Modules
    # from ase import io
    # #__|
    #
    # #| - Script Inputs
    # atoms_A_filename = "A.traj"
    # atoms_B_filename = "B.traj"
    # #__|
    #
    # #| - Reading Atoms Objects
    # atoms_A = io.read(atoms_A_filename)
    # atoms_B = io.read(atoms_B_filename)
    # #__|

    #__|

    #| - Building the Identical Atom Index List for Both Atoms Objects
    atoms_A_ind_list = []
    atoms_B_ind_list = []
    for atom_A in atoms_A:
        for atom_B in atoms_B:

            # Comparing Positions
            pos_A = atom_A.position
            pos_B = atom_B.position
            pos_AB_comp = pos_A == pos_B

            # Comparing Element Type
            elem_A = atom_A.symbol
            elem_B = atom_B.symbol
            elem_AB_comp = elem_A == elem_B

            if all(pos_AB_comp) is True and elem_AB_comp is True:
                atoms_A_ind_list.append(atom_A.index)
                atoms_B_ind_list.append(atom_B.index)
    #__|

    atoms_A_unique_ind_list = []
    for atom_A in atoms_A:
        if atom_A.index not in atoms_A_ind_list:
            atoms_A_unique_ind_list.append(atom_A.index)

    atoms_B_unique_ind_list = []
    for atom_B in atoms_B:
        if atom_B.index not in atoms_B_ind_list:
            atoms_B_unique_ind_list.append(atom_B.index)

    return(atoms_A_unique_ind_list, atoms_B_unique_ind_list)

    #__|

#__| **************************************************************************

#| - Modify Atoms Object ******************************************************

def move_atoms_of_element_i(atoms, element, new_position, dim="z"):
    """Modify positions of all atoms of certain element.

    Args:
        atoms:
        element:
        new_position:
        dim:
    """
    #| - move_atoms_of_element_i
    if type(atoms) == str:
        atoms = read(atoms)
    else:
        pass

    #| - Converting Input Dimension to Integer Index
    if dim == "x":
        dim = 0
    elif dim == "y":
        dim = 1
    elif dim == "z":
        dim = 2
    #__|

    elem_i_list = []
    for atom in atoms:
        if atom.symbol == element:

            atom.position[dim] = atom.position[dim] + new_position
            elem_i_list.append(atom)

    # return(atoms)
    #__|

def displace_overlayer(
    atoms,
    x_ind,
    y_ind,
    mesh_size_x,
    mesh_size_y,
    element="C",
    save_file=False,
    ):
    """Displace atoms in an overlayer.

    Args:
        x_ind:
        y_ind:
        mesh_size_x:
        mesh_size_y:
        element:
        save_file:
    """
    #| - displace_overlayer
    atoms = copy.deepcopy(atoms)

    x_frac = 1. * x_ind / mesh_size_x
    y_frac = 1. * y_ind / mesh_size_y

    # atoms = io.read("dir_atoms/init.traj")
    # mag_lv0 = np.linalg.norm(atoms.cell[0])
    # mag_lv1 = np.linalg.norm(atoms.cell[1])

    lv0 = atoms.cell[0]
    lv1 = atoms.cell[1]


    disp = lv0 * x_frac + lv1 * y_frac

    move_atoms_of_element_i(atoms, element, disp[0], dim="x")
    move_atoms_of_element_i(atoms, element, disp[1], dim="y")

    fle_name = str(x_ind).zfill(2) + "_" + str(y_ind).zfill(2) + "_new.traj"

    if save_file:
        atoms.write(fle_name)

    return(atoms)
    #__|

def change_vacuum(atoms, vacuum):
    """Change the amount of vacuum in a slab.

    Assumes that the vacuum is in the z-direction orthogonal to
    x and y axis

    Args:
        atoms:
        vacuum:
    """
    #| - change_vacuum
    if type(atoms) == str:
        atoms = read(atoms)
    else:
        pass

    cell_copy = copy.deepcopy(atoms.cell)

    cell_copy[2][2] = vacuum

    atoms.set_cell(cell_copy)
    atoms.center()

    return(atoms)
    #__|

#__| **************************************************************************

#| - Atoms Information Methods ************************************************


def number_of_constrained_atoms(atoms):
    """Count number of constrained atoms in atoms object.

    Args:
        atoms:
    """
    #| - number_of_constrained_atoms
    if type(atoms) == str:
        atoms = read(atoms)
    else:
        pass

    N_constraints = len(atoms.constraints)

    return(N_constraints)
    #__|

def highest_position_of_element(atoms, element_symbol):
    """Return highest z-value for given element type.

    Args:
        atoms:
        element_symbol
    """
    #| - highest_position_of_element

    #| - SCRIPT INPUTS
    # element_symbol
    element_name = element_symbol
    # atoms_file_name = atoms
    #__|

    if type(atoms) == str:
        atoms = read(atoms)
    else:
        pass

    elem_atom_list = []
    for atom_i in atoms:
        if atom_i.symbol == element_name:
            elem_atom_list.append(atom_i.position)

    elem_atom_list = np.array(elem_atom_list)
    elem_atom_list = np.sort(elem_atom_list[:, 2])
    highest_z_pos = elem_atom_list[-1]

    return(highest_z_pos)
    #__|

def number_of_atoms(atoms):
    """Return atom count dictionary.

    DEPRECATED DEPRECATED DEPRECATED DEPRECATED DEPRECATED DEPRECATED
    create_species_element_dict

    Args:
        atoms
    """
    #| - number_of_atoms
    atoms_sym_list = atoms.get_chemical_symbols()
    unique_atom_symbols = list(set(atoms_sym_list))

    atom_dict = {}
    for atom_sym_i in unique_atom_symbols:
        atom_i_cnt = atoms_sym_list.count(atom_sym_i)

        atom_dict[atom_sym_i] = atom_i_cnt

    print("THIS HAS BEEN DEPRECATED to create_species_element_dict")
    print("THIS HAS BEEN DEPRECATED to create_species_element_dict")
    print("THIS HAS BEEN DEPRECATED to create_species_element_dict")
    print("THIS HAS BEEN DEPRECATED to create_species_element_dict")

    # return(atom_dict)
    #__|

def create_species_element_dict(
    atoms,
    include_all_elems=False,
    elems_to_always_include=None,
    ):
    """Create dict from an atoms object with element: element number entries.

    If 'include_all_elems' is True then 'elems_to_always_include' must be None

    Args:
        atoms:
        include_all_elems: <boolean> or <list>
            False: Does not add additional elements other than the ones
            present in the atoms object

            True: Includes all elements in the periodic table, with 0 values
            for the elements not present in the atoms object

        elems_to_always_include:
            List: List of elements to include in the final dict, not including
            the elements present in the atoms object.

    """
    #| - create_species_element_dict
    from misc_modules.misc_methods import merge_two_dicts

    all_elements = list(periodic_table_dict)

    chem_syms = atoms.get_chemical_symbols()
    chem_syms_unique = set(chem_syms)

    species_elem_dict = {}
    for elem_i in chem_syms_unique:
        num_elem_i = chem_syms.count(elem_i)

        species_elem_dict[elem_i] = num_elem_i

    #| - Include All Elements in the periodic table
    if include_all_elems or elems_to_always_include is not None:
        all_non_occuring_elements = list(
            filter(
                lambda x: x not in set(list(species_elem_dict)), all_elements
                )
            )
#         print(all_non_occuring_elements)

#         if elems_to_always_include is not None:
        if elems_to_always_include is not None and type(elems_to_always_include) == list:
            non_occuring_elements = [i for i in all_non_occuring_elements if i in elems_to_always_include]
        else:
            non_occuring_elements = all_non_occuring_elements

        non_occuring_species_elem_dict = dict(
            zip(
                non_occuring_elements,
                [0 for i in non_occuring_elements],
                )
            )

#         non_occuring_species_elem_dict = dict(
#             zip(
#                 all_non_occuring_elements,
#                 [0 for i in all_non_occuring_elements],
#                 )
#             )


        species_elem_dict = merge_two_dicts(
            non_occuring_species_elem_dict,
            species_elem_dict,
            )
    #__|

    return(species_elem_dict)
    #__|

#__| **************************************************************************

#| - Visualization ************************************************************
def create_gif_from_traj(
    traj_name="qn.traj",
    path_i=".",
    image_range="0:5",  # "all"
    delay=10,
    rotation="0x, 0y, 0z",
    ):
    """Create gif animation from trajectory file.

    Args:
        traj_name:
        path_i:
        image_range:
        delay:

    """
    #| - create_gif_from_traj

    #| - Method  Parameters
    atoms_file_name = "out_movie"
    fold_name = "images"
    #__|

    #| - Creating png *********************************************************
    if image_range == "all":
        ind_range = ":"
    else:
        ind_range = image_range

    atoms_traj = io.read(path_i + "/" + traj_name, index=ind_range)

    folder_i = fold_name
    if not os.path.isdir(path_i + folder_i):
        os.makedirs(path_i + "/" + folder_i)

    for index, atoms_i in enumerate(atoms_traj):
        textures = ["simple" for i in range(len(atoms_i))]

        name = "TMP"
        name_i = fold_name + "/" + str(index).zfill(3) + "_" + name + "_top"
        io.write(
            path_i + "/" + name_i + ".pov",
            atoms_i,
            rotation=rotation,
            run_povray=True,
            textures=textures,
            canvas_height=1000,
            display=False
            )

        print("Creating IMAGE: " + path_i + "/" + name_i + ".png")

        os.remove(path_i + "/" + name_i + ".pov")
        os.remove(path_i + "/" + name_i + ".ini")
    #__| **********************************************************************

    #| - Converting Images to GIF *********************************************
    root_dir = os.getcwd()
    os.chdir(path_i + "/" + fold_name)

    bash_command = "/usr/bin/convert "
    bash_command += "-delay " + str(delay) + " "
    bash_command += "-loop 0 "
    # bash_command += "*png animation.gif"
    bash_command += "*png " + atoms_file_name + ".gif"

    print("###### CREATING GIF ####### | " + atoms_file_name)
    os.system(bash_command)
    os.system("mv *.gif ..")
    os.chdir(root_dir)
    #__| **********************************************************************

    #__|

def create_gif_from_atoms_movies(
    atoms_file="Default",
    path_i=".",
    delay=10,
    ):
    """Create png images from an multi-atoms atoms object.

    TODO This method is way to specific, should be chopped up
    TODO Use create_gif_from_traj method

    Args:
        atoms_file:
        path_i:
        delay:
    """
    #| - create_images_from_atoms_movies

    #| - SCRIPT PARAMETERS
    fold_name = "images"
    #__|

    #| - Read Atoms File with *.traj File Name

    if atoms_file == "Default":
        filenames = next(os.walk(path_i))[2]

        traj_files = [fle for fle in filenames if ".traj" in fle]
        alist = {}
        if len(traj_files) == 1:
            atoms_file = traj_files[0]

            atoms_file_path = path_i + "/" + atoms_file

            # alist[atoms_file] = read(atoms_file, index=":")

    else:
        pass

    atoms_file_name = copy.deepcopy(atoms_file)
    atoms_file_name = atoms_file_name.split(".")
    atoms_file_name.remove("traj")
    atoms_file_name = ".".join(atoms_file_name)

    alist[atoms_file] = read(atoms_file_path, index=":")
    #__|

    folder_i = fold_name
    if not os.path.isdir(path_i + folder_i):
        os.makedirs(path_i + folder_i)

    for i, file in enumerate(alist):
        for index, atoms_i in enumerate(alist[file]):
            textures = ["simple" for i in range(len(atoms_i))]

            if "traj" in file:
                name = file[0:-5]
            elif "xyz" in file:
                name = file[0:-4]
            elif "POSCAR" in file:
                name = file[0:-7]

            name_i = fold_name + "/" + str(index).zfill(2) + "_" + name + "_top"
            write(
                path_i + "/" + name_i + ".pov",
                atoms_i,
                # alist[file],
                rotation="0x, 0y, 0z",
                run_povray=True,
                textures=textures,
                canvas_height=1000,
                display=False
                )

            print("Creating IMAGE: " + path_i + name_i + ".png")

            os.remove(path_i + "/" + name_i + ".pov")
            os.remove(path_i + "/" + name_i + ".ini")


    #| - Converting Images to GIF
    root_dir = os.getcwd()

    os.chdir(path_i + "/" + fold_name)

    bash_command = "/usr/bin/convert "
    bash_command += "-delay " + str(delay) + " "
    bash_command += "-loop 0 "
    # bash_command += "*png animation.gif"
    bash_command += "*png " + atoms_file_name + ".gif"

    print("###### CREATING GIF ####### | " + atoms_file_name)
    os.system(bash_command)

    os.system("mv *.gif ..")

    os.chdir(root_dir)
    #__|

    # TODO - Remove png files after creating gif !!!!!
    #__|

#__| **************************************************************************

#| - MISC

def max_force(atoms):
    """Return largest force on any atom.

    Args:
        atoms
    """
    #| - max_force
    forces = atoms.get_forces()

    sum = 0.0
    largest = 0.0
    for a in range(len(atoms)):
        force = np.sqrt(
            forces[a][0] ** 2 +
            forces[a][1] ** 2 +
            forces[a][2] ** 2
            )

        sum += force
        if(force > largest):
            largest = force

    return(largest, sum)
    #__|



#__|
