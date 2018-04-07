#!/usr/bin/env python

"""Methods for ASE scripts, mostly DFT scripts.

Author: Raul A. Flores

Development Notes:
    TODO Master print statements should have "****..." to be easily readable
    TODO Delete all set_mag_mom_to_0 references (Already commented them out)
"""

#| - Table of Contents
"""
set_QE_calc_params:
ionic_opt:
set_mag_mom_to_0:
increase_abs_val_magmoms:
calc_spinpol:
simple_mag_moms:
set_init_mag_moms:
reduce_magmoms:
an_pdos:
spin_pdos:
an_bands:
an_beef_ensemble:
plot_beef_ensemble:
an_ads_vib:
thermochemical_corrections:
read_atoms_from_file:
convert_atoms_object:
angle_between_lattice_vectors:
magnitude_of_lattice_vectors:
move_atoms_of_element_i:
displace_overlayer:
change_vacuum:
number_of_atoms:
number_of_constrained_atoms:
highest_position_of_element:
create_gif_from_atoms_movies:
"""
#__|

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
import matplotlib.pyplot as plt

from ase.io import read, write
from ase.dft.kpoints import ibz_points, get_bandpath

from ase.vibrations import Vibrations
from ase.thermochemistry import HarmonicThermo

# My Modules
from misc_modules.numpy_methods import angle_between
from ase_modules.dft_params import Espresso_Params
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
    atoms,
    params={},
    load_defaults=True,
    ):
    """Set Quantum Espresso calculation parameters to atoms object.

    Handles reading, and setting of dft calculator parameters to atoms object.

    Args:
        atoms:
        params:
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

    calc = espresso(**espresso_params)
    atoms.set_calculator(calc=calc) # NOTE Calculator set for the 1st time here

    return(calc, espresso_params)
    #__|

#__| **************************************************************************

#| - Ionic Optimization *******************************************************

def ionic_opt(
    atoms,
    calc,
    espresso_params=None,
    mode="opt",
    fmax=0.05,
    ):
    """
    Run ionic dft relaxation on atoms object.


    Development Notes:
        * What to do when the job has already been previously completed?
        Should I run a single point calculation just to make sure the calculator is fully
        operational/loaded?

        TODO .FINSISHED file should include information about what kind of optimization was
        performed

    Args:
        atoms:
        calc:
        espresso_params:
        mode:
        fmax:
    """
    #| - ionic_opt
    from espresso import espresso
    from ase.optimize import QuasiNewton

    mess = "Running DFT Calculation "
    mess += "******************************************************"
    print(mess)
    sys.stdout.flush()

    #| - Checking if Previous Calculation Has Been Completed
    filename = ".FINISHED.new"
    if os.path.exists("./" + filename):
        with open(filename, "r") as fle:
            lines = [line.strip() for line in fle.readlines()]
            if "ionic_opt" in lines:
                print("ionic_opt | Optimization previusly completed "
                    "(Running a single-point calculation)"
                    )
                mode = "sp"
    #__|

    atoms.set_calculator(calc)
    # TODO Skip calculation if previously converged

    if mode == "opt":
        #| - Regular Optimization
        qn = QuasiNewton(
            atoms,
            trajectory="out_opt.traj",
            logfile="qn.log",
            )

        if os.path.exists("prev.traj"):
            qn.replay_trajectory("prev.traj")

        qn.run(fmax=fmax)
        #__|

    elif mode == "easy_opt":
        #| - Easy Optimization -> Full Optimization

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
        # FIXME I should be able to just use the "set_initial_magnetic_moments" method here
        # That should automatically set to magmoms to 0 from the QE parameters dict
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
        atoms.set_calculator(calc)

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
        print("ionic_opt | Running Single-Point Calculation"); sys.stdout.flush()

        atoms.get_potential_energy()
        write("out.traj", atoms)
        #__|

    #| - Update .FINISHED File

    update_FINISHED("ionic_opt")

    # filename = ".FINISHED.new"
    #
    # if os.path.exists("./" + filename):
    #     append_write = "a"  # append if already exists
    # else:
    #     append_write = "w"  # make a new file if not
    #
    # with open(filename, append_write) as fle:
    #     fle.write("ionic_opt")
    #     fle.write("\n")
    #__|

#__| **************************************************************************

#__| **************************************************************************

#| - Magnetic Moments *********************************************************

def set_init_mag_moms(atoms, preference="bader", magmoms=None):
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
        magmoms:
            If specified, set the atom's initial magnetic moments to  this
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
        print("set_init_mag_moms | Using given magmoms"); sys.stdout.flush()
        magmoms_i = magmoms
    else:
        spinpol_calc = calc_spinpol(atoms)

        #| - Spin-polarization turned off, set magmoms to 0
        if not spinpol_calc:
            print("set_init_mag_moms | Spin-polarization turned off"); sys.stdout.flush()
            mag_mom_list = atoms.get_initial_magnetic_moments()
            magmoms_i = np.zeros(len(mag_mom_list))
        #__|


        # COMBAK This is a poor way of enforcing the analysis method preference
        # Assumes that there are only 2 analysis methods
        elif preferred_an in magmoms_master_dict.keys():
            text = (
                        "set_init_mag_moms | "
                        "Using preferred method for initial magnetic moments "
                        "(" + preferred_an + ")"
                        )
            print(text); sys.stdout.flush()
            magmoms_i = magmoms_master_dict[preferred_an]

        elif len(magmoms_master_dict.keys()) == 1:
            an_method = magmoms_master_dict.keys()[0]
            text = ("set_init_mag_moms | "
                    "Using " + an_method + " for initial magnetic moments")
            print(text);

            magmoms_i = magmoms_master_dict[an_method]

        #| - __old__
        # elif preference == "bader":
        #     if "bader_magmoms" in magmoms_master_dict.keys():
        #         text = ("set_init_mag_moms | "
        #                 "Initial magnetic moments set from bader analysis")
        #         print(text)
        #         magmoms_i = magmoms_master_dict["bader_magmoms"]
        #
        # elif preference == "pdos" and "pdos_magmoms" in magmoms_master_dict.keys():
        #     text = ("set_init_mag_moms | "
        #             "Initial magnetic moments set from PDOS analysis")
        #     print(text)
        #     magmoms_i = magmoms_master_dict["pdos_magmoms"]
        #__|


        #| - Use simple method
        else:
            text = ("set_init_mag_moms | "
                    "Using simple method for initial magnetic moments")
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

# COMBAK | Remove this method, it is too simple
# def set_mag_mom_to_0(atoms):
#     """Set magnetic moments to 0 for all atoms in atoms object.
#
#     Args:
#         atoms:
#     """
#     #| - set_mag_mom_to_0
#     mag_mom_list = atoms.get_initial_magnetic_moments()
#     new_mag_mom_list = np.zeros(len(mag_mom_list))
#     atoms.set_initial_magnetic_moments(new_mag_mom_list)
#     #__|

# def increase_abs_val_magmoms(magmoms_list, increase_amount=0.2):
def increase_abs_val_magmoms(atoms, magmoms_list, increase_amount=0.2):
    """Increase absolute value of magmoms for atoms object.

    # COMBAK Shouldn't raise initial guess for light atoms (Or don't raise as
    much)

    Select light atoms will have their |magmom| raised by a set small amount to
    not oversaturate the atomic magnetic moment (init_magmom > # valence e)

    Args:
        magmoms_list:
        increase_amount:
    """
    #| - increase_abs_val_magmoms
    inc = increase_amount

    light_atoms_list = {
        "H": 0.2,
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
        if atom.symbol in light_atoms_list.keys():
            inc = light_atoms_list[atom.symbol]

        if magmom < 0.:
            new_magmom = abs(magmom) + inc
            new_magmom_list.append(-new_magmom)

        elif magmom > 0.:
            new_magmom = abs(magmom) + inc
            new_magmom_list.append(new_magmom)

        else:
            new_magmom_list.append(0.)

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
        elif atom.symbol == "O":
            magmoms[atom.index] = 2.5
        else:
            magmoms[atom.index] = 0.0

    return(np.array(magmoms))
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
    mess = "Running PDOS Analysis "
    mess += "********************************************************"
    print(mess); sys.stdout.flush()


    outdir = "dir_pdos"
    # atoms.set_calculator(calc=calc)

    dos = atoms.calc.calc_pdos(
        nscf=True,
        kpts=dos_kpts,
        Emin=-15.0,
        Emax=15.0,
        ngauss=0,
        sigma=0.2,
        DeltaE=0.01,
        tetrahedra=False,
        slab=True,
        )

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
        "Cu": 11, "C": 4, "O": 6, "H": 1,
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
            ##Update charge
            charge_i = nvalence_dict[atom.symbol] - (charge)
            if write_charges:
                # atom.charge = nvalence_dict[atom.symbol] - (charge)
                atom.charge = charge_i
            charge_list.append(charge_i)
            #__|

            atoms.info.update({"pdos_charges": charge_list})

            pickle.dump(charge_list, open("%s/charge_list.pickle" % outdir, "w"))
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

    # FIXME This should be handled by envoking the set_initial_magnetic_moments method
    spinpol_calc = calc_spinpol(atoms)
    if spinpol_calc is False:
        print("set_init_mag_moms | Spin-polarization turned off")
        atoms.set_initial_magnetic_moments(np.zeros(len(atoms)))
        # set_mag_mom_to_0(atoms)  # COMBAK Remove this if working

    # calc_spinpol(atoms)  # COMBAK I don't think this is doing anything

    espresso_params.update(
        {
            "kpts": bands_kpts,
            "outdir": "dir_bands"
            }
        )

    atoms.calc = espresso(**espresso_params)
    atoms.get_potential_energy()

    atoms.calc.save_flev_chg("charge_den.tgz")
    atoms.calc.load_flev_chg("charge_den.tgz")

    # COMBAK How to handle this more generally?
    ip = ibz_points["fcc"]
    points = ["Gamma", "X", "W", "K", "L", "Gamma"]

    bzpath = [ip[p] for p in points]

    kpts, x, X = get_bandpath(bzpath, atoms.cell, npoints=300)
    energies = atoms.calc.calc_bandstructure(kpts)

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
    mess += "************************************************"
    print(mess); sys.stdout.flush()

    beefensemble_on = atoms.calc.beefensemble
    xc = atoms.calc.xc

    if beefensemble_on is False:
        print("The ase-espresso calculator has the beefensemble"
                            " parameter turned off!, analysis will not run!!")
        return(None)

    if xc != "BEEF" or xc != "BEEF-vdW":
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
def an_ads_vib(atoms, ads_index_list=None):
    """Adsorbate vibrational analysis.

    Args:
        atoms:
        ads_index_list:
    """
    #| - an_ads_vib

    mess = "Starting vibrational analysis "
    mess += "************************************************"
    print(mess)

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

    # if len(ads_index_list) == 0:
    #     print("an_ads_vib | len(ads_index_list) == 0, something wrong")
    #     return(None)
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

    set_init_mag_moms(atoms)
    vib = Vibrations(atoms, indices=ads_index_list)
    vib.run()

    print("an_ads_vib | Getting vibrational energies")
    vib_e_list = vib.get_energies()

    # COMBAK I think that the summary is succesfully being outputed to
    # vib_summ.out, in which case all of these prints statements can be deleted
    print("an_ads_vib | printing summary??"); sys.stdout.flush()
    tmp = vib.summary(log="vib_summ.out")
    print(tmp)  # TEMP Just trying to get this info to print
    print("an_ads_vib | printing summary??"); sys.stdout.flush()

    #| - Copy Files to dir_vib Folder
    if not os.path.exists("dir_vib"):
        os.makedirs("dir_vib")
    shutil.move("vib_summ.out", "dir_vib/vib_summ.out")

    dest_dir = "dir_vib"
    for fle in glob.glob(r'*.pckl*'):
        shutil.move(fle, dest_dir + "/" + fle)
    #__|

    thermochemical_corrections(vib_e_list)

    update_FINISHED("an_ads_vib")

    return(vib)
    #__|

def thermochemical_corrections(vib_e_list, Temperature=300.0):
    """Thermochemical free energy corrections from vibrational analysis.

    Args:
        vib_e_list:
            List of vibrational modes in eV
        Temperature:
    """
    #| - thermochemical_corrections
    print("Calculating thermochemical corrections @ " + str(Temperature) + "K")

    # Remove imaginary frequencies
    vib_e_list = vib_e_list.real
    vib_e_list = [vib_i for vib_i in vib_e_list if vib_i != 0.]

    print("Vibrations list:")
    print(vib_e_list)

    ht = HarmonicThermo(vib_e_list)
    F_energy = ht.get_helmholtz_energy(Temperature)

    if not os.path.exists("dir_vib"):
        os.makedirs("dir_vib")

    with open("dir_vib/gibbs_corr.out", "w") as fle:
        fle.write(str(F_energy))
        fle.write("\n")
    #__|

#__| **************************************************************************


#| - Atoms File Operations ****************************************************
def read_atoms_from_file(filename=None):
    """Read atoms object from file.

    Checks several file names

    Args:
        filename: optional atoms file, will attempt to read first.
    """
    #| - read_atoms_from_file
    mess = "Reading Atoms Object From File "
    mess += "***********************************************"
    print(mess)

    atoms = None
    file_name_list = [
        "init.traj",
        "init.POSCAR",
        "manual_start.traj",
        "out.traj",
        "POSCAR",
        ]

    if filename is not None:
        file_name_list.insert(0, filename)

    for file_name in file_name_list:
        if os.path.isfile(file_name):
            atoms = read(file_name)
            break

    if atoms is None:
        raise IOError("No atoms file found")

    #| - OLD
    #
    # file_name = "init.traj"
    # if os.path.isfile(file_name):
    #     atoms = read(file_name)
    #
    # file_name = "init.POSCAR"
    # elif os.path.isfile(file_name):
    #     atoms = read(file_name)
    #
    # file_name = "manual_start.traj"
    # elif os.path.isfile(file_name):
    #     atoms = read(file_name)
    #
    # file_name = "out.traj"
    # elif os.path.isfile(file_name):
    #     atoms = read(file_name)
    #
    # file_name = "POSCAR"
    # elif os.path.isfile(file_name):
    #     atoms = read(file_name)
    #
    # if atoms == None:
    #     raise IOError("No atoms file found")



    # elif os.path.isfile("init.POSCAR"):
    #     atoms = read("init.POSCAR")
    #
    # elif os.path.isfile("init.traj"):
    #     atoms = read("init.traj")
    #
    # elif os.path.isfile("manual_start.traj"):
    #     atoms = read("manual_start.traj")
    #
    # elif os.path.isfile("out.traj"):
    #     atoms = read("out.traj")
    #
    # elif os.path.isfile("POSCAR"):
    #     atoms = read("POSCAR")
    #
    # if atoms == None:
    #     raise IOError("No atoms file found")
    #__|

    return(atoms)
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

def number_of_atoms(atoms):
    """Return atom count dictionary.

    Args:
        atoms
    """
    #| - number_of_atoms
    # atoms = io.read(path + "/out.traj")

    atoms_sym_list = atoms.get_chemical_symbols()

    unique_atom_symbols = list(set(atoms_sym_list))

    atom_dict = {}
    for atom_sym_i in unique_atom_symbols:
        atom_i_cnt = atoms_sym_list.count(atom_sym_i)

        atom_dict[atom_sym_i] = atom_i_cnt

    return(atom_dict)
    #__|

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

#__| **************************************************************************

#| - Visualization ************************************************************
def create_gif_from_atoms_movies(
    atoms_file="Default",
    path_i=".",
    delay=10,
    ):
    """Create png images from an multi-atoms atoms object.

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
            textures = ['simple' for i in range(len(atoms_i))]

            if 'traj' in file:
                name = file[0:-5]
            elif 'xyz' in file:
                name = file[0:-4]
            elif 'POSCAR' in file:
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
