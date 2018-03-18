#!/usr/bin/env python
"""Methods for ASE scripts, mostly DFT scripts.
I'm adding this here and pushing to repo - 180317"""

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
from ase.io.trajectory import Trajectory
from ase.dft.kpoints import ibz_points, get_bandpath

from ase.vibrations import Vibrations
from ase.thermochemistry import HarmonicThermo

# My Modules
from misc_modules.numpy_methods import angle_between
from ase_modules.dft_params import Espresso_Params
#__|


#| - Parse DFT Job Parameters


def set_QE_calc_params(atoms, params={}):
    """
    """
    #| - set_QE_calc_params
    from espresso import espresso

    print("Loading QE parameters from file")
    espresso_params_inst = Espresso_Params(load_defaults=True)

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
    atoms.set_calculator(calc=calc)

    return(calc, espresso_params)
    #__|

#__|

#| - Ionic Optimization


def ionic_opt(
    atoms,
    calc,
    espresso_params,
    mode="opt",
    fmax=0.05,
    ):

    """
    """
    #| - ionic_opt
    from espresso import espresso
    from ase.optimize import QuasiNewton

    print("Running DFT calculation"); sys.stdout.flush()


    atoms.set_calculator(calc)

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
        print("Running Easy Relaxation"); sys.stdout.flush()

        magmoms = atoms.get_initial_magnetic_moments()

        set_mag_mom_to_0(atoms)
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
        print("Running Full Relaxation"); sys.stdout.flush()
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
        atoms.get_potential_energy()
        write("out.traj", atoms)
        #__|

    #__|

#__|

#| - Magnetic Moments *********************************************************


def set_mag_mom_to_0(atoms):
    """
    """
    #| - set_mag_mom_to_0
    mag_mom_list = atoms.get_initial_magnetic_moments()
    new_mag_mom_list = np.zeros(len(mag_mom_list))
    atoms.set_initial_magnetic_moments(new_mag_mom_list)
    #__|


# def increase_abs_val_magmoms(magmoms_list, increase_amount=1.1):
def increase_abs_val_magmoms(magmoms_list, increase_amount=0.2):
    """
    """
    #| - increase_abs_val_magmoms
    inc = increase_amount

    new_magmom_list = []
    for magmom in magmoms_list:
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
    """
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
    """
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
        "Ta": 3, "W":  4, "Re": 5, "Ti": 2,
        }
    #__|

    #| - Find cations
    cations = []
    for atom in atoms:
        if master_dict_high_spin.has_key(atom.symbol):
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


def set_init_mag_moms(atoms, preference="bader", magmoms=None):
    """
    Set inital magnetic moments to atoms object. If the atoms object has
    previously had a bader or pdos analysis performed those magnetic moments
    will be available under the atoms.info dict ("pdos_magmoms" and
    "bader_magmoms" dict keys).

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
    print("Setting inital magnetic moments")

    magmom_keys = ["pdos_magmoms", "bader_magmoms"]
    magmoms_master_dict = {}
    for magmom_key in magmom_keys:
        if magmom_key in atoms.info.keys():
            magmoms_master_dict[magmom_key] = atoms.info[magmom_key]

    magmom_keys = magmoms_master_dict.keys()

    spinpol_calc = calc_spinpol(atoms)
    if not spinpol_calc:
        print("set_init_mag_moms | Spin-polarization turned off")
        set_mag_mom_to_0(atoms)

    elif magmoms is not None:
        magmoms_i = magmoms

    elif preference == "bader" and "bader_magmoms" in magmom_keys:
        text = ("set_init_mag_moms | "
                "Initial magnetic moments set from bader analysis")
        print(text)
        magmoms_i = magmoms_master_dict["bader_magmoms"]
        # atoms.set_initial_magnetic_moments(magmoms_master_dict["bader_magmoms"])

    elif preference == "pdos" and "pdos_magmoms" in magmom_keys:
        text = ("set_init_mag_moms | "
                "Initial magnetic moments set from PDOS analysis")
        print(text)
        magmoms_i = magmoms_master_dict["pdos_magmoms"]
        # atoms.set_initial_magnetic_moments(magmoms_master_dict["pdos_magmoms"])

    else:
        text = ("set_init_mag_moms | "
                "Using simple method for initial magnetic moments")
        print(text)
        magmoms_i = simple_mag_moms(atoms)
        # atoms.set_initial_magnetic_moments(magmoms_simple)

    magmoms_i_new = increase_abs_val_magmoms(magmoms_i)
    atoms.set_initial_magnetic_moments(magmoms_i_new)
    #__|


def reduce_magmoms(atoms, ntypx=10):
    """
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
            print "WARNING, reducing pair of magmoms whose difference is %.2f"%min_delta

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
    """
    """
    #| - an_pdos
    print("Running PDOS analysis"); sys.stdout.flush()
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

    if not os.path.exists("dir_pdos"):
        os.makedirs("dir_pdos")

    pdos_out = "dir_pdos/dos.pickle"
    with open(pdos_out, "w") as fle:
        pickle.dump(dos, fle)

    # Set Magnetic Moments To Atoms Object From PDOS Intergration
    spin_pdos(atoms, pdos_out, spinpol=espresso_params["spinpol"])
    # write("dir_pdos/pdos.traj", atoms)

    atoms.write("out.traj")
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
    """
    Calculate more accurate charges/magnetic moments from pdos and assign to
    atom.charge/atom.magmom. If pdos_pkl not defined, it will be calculated
    (atoms object must have a real calculator attached). If pdos_pkl is
    defined, will attempt to load from pickle file, but valence_dict must be
    specified! Specify calculation directory as outdir for easy cleanup.
    """
    #| - spin_pdos

    valence_dict = {
        "Cu": 11, "C": 4,  "O": 6,  "H": 1,
        "Rh":17,  "Co":9,  "Pd":10, "Pt":10,
        "Ni":1,   "Fe":16, "N":5,   "Ru":16,
        }

    #| - Reading PDOS File, Otherwise Creates It
    if pdos_pkl:
        assert valence_dict is not None, "MUST SPECIFY valence_dict"
        single_point_calc = True
        pdos = pickle.load(open(pdos_pkl))
        nvalence_dict=valence_dict
    else:
        single_point_calc = False
        nvalence_dict = atoms.calc.get_nvalence()[1] #dict with (chemical symbol) : (num valence)
        if nscf: #double k-points for higher quality pdos --> O(single point calc)
            pdos = atoms.calc.calc_pdos(nscf=True,kpts=kpts,**kwargs)
        else: #no single point calc, should take 1 or 2 minutes
            pdos = atoms.calc.calc_pdos(**kwargs)
    #__|

    #| - Finding Index of Fermi Level
    for i_ind, e_i in enumerate(pdos[0]):
        if e_i > 0:
            fi = i_ind #index of fermi level
            break
    #__|

    #| - Analysing PDOS For Magnetic Moments and Charge of All Atoms
    if spinpol:
        #| - Spin Polarlized Calculation
        magmom_list = []
        charge_list = []
        for i, atom in enumerate(atoms):

            #| - Integrating Up and Down Spin PDOS
            spin_up = 0; spin_down = 0;
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
                # atom.charge = nvalence_dict[atom.symbol] - (spin_up + spin_down)
                atom.charge = charge_i

            charge_list.append(charge_i)
            #__|

        print("PDOS MAGMOMS: " + str(atoms.get_initial_magnetic_moments()))
        reduce_magmoms(atoms)

        atoms.info.update({"magmom_set": True})
        atoms.info.update({"pdos_magmoms": magmom_list})
        atoms.info.update({"pdos_charges": charge_list})
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

        f = open("%s/pdos_Lowdin.log"%outdir,"w")
        f.writelines(light_lines)
        f.close()
        os.system("rm %s/pdos.log"%outdir)

    if save_pkl:
        pickle.dump(pdos, open("pdos.pkl", "w"))
    #__|

    #__|

#__| **************************************************************************

#| - Band Structure ***********************************************************
def an_bands( atoms, bands_kpts, espresso_params, ):
    """ """
    #| - an_bands
    from espresso import espresso

    print("Executing Band Structure Analysis"); sys.stdout.flush()

    spinpol_calc = calc_spinpol(atoms)
    if spinpol_calc == False:
        print("set_init_mag_moms | Spin-polarization turned off")
        set_mag_mom_to_0(atoms)

    calc_spinpol(atoms)

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

    ip = ibz_points["fcc"]
    points = ["Gamma", "X", "W", "K", "L", "Gamma"]

    bzpath = [ip[p] for p in points]

    kpts, x, X = get_bandpath(bzpath, atoms.cell, npoints=300)
    energies = atoms.calc.calc_bandstructure(kpts)

    if not os.path.exists("dir_bands"):
        os.makedirs("dir_bands")

    with open("dir_bands/band_disp.pickle", "w") as fle:
        pickle.dump((points, kpts, x, X, energies), fle)

    shutil.move("charge_den.tgz", "dir_bands")

    atoms.write("out.traj")
    #__|

#__| **************************************************************************

#| - Beef Ensemble of Energies ************************************************
def an_beef_ensemble(atoms, xc):
    """
    """
    #| - an_beef_ensemble
    if xc == "BEEF" or xc == "BEEF-vdW":
        calc = atoms.calc
        from ase.dft.bee import BEEFEnsemble

        energy = atoms.get_potential_energy()
        ens = BEEFEnsemble(atoms=atoms, e=energy, xc="beefvdw")
        ens_e = ens.get_ensemble_energies()

        if not os.path.exists("dir_beef_ensemble"):
            os.makedirs("dir_beef_ensemble")
        else: pass

        with open("dir_beef_ensemble/ensemble.pickle", "w") as fle:
            pickle.dump(ens_e, fle)
    #__|

def plot_beef_ensemble(
    folder_dir="dir_beef_ensemble",
    file_name="ensemble.pickle",
    file_out="beef_ens_hist.png",
    ):
    """
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

def an_ads_vib(atoms, ads_index_list=[]):
    """

    """
    #| - an_ads_vib
    if len(ads_index_list) == 0:
        print("Must define adsorbate(s) index(s)")
        return(None)

    # ads_index_list = [1]
    # FIXME - Would be nice to remove pickle files that are empty

    #| - Copy vib.pckl files back to root dir (For restarting)
    for fle in glob.glob(r'dir_vib/*.pckl*'):
        fle_name = fle.split("/")[-1]
        shutil.move(fle, fle_name)
    #__|

    print("Starting vibratinal analysis")
    set_init_mag_moms(atoms)

    # list indices of atoms to be vibrated, this one is particular to OH*,
    # i. e. [80,81]
    vib = Vibrations(atoms, indices=ads_index_list)
    vib.run()
    vib_e_list = vib.get_energies()

    vib.summary(log="vib_summ.out")

    #| - Copy Files to dir_vib Folder
    if not os.path.exists("dir_vib"):
        os.makedirs("dir_vib")
    shutil.move("vib_summ.out", "dir_vib/vib_summ.out")

    dest_dir = "dir_vib"
    for fle in glob.glob(r'*.pckl*'):
        # print(fle)
        shutil.move(fle, dest_dir + "/" + fle)
    #__|

    thermochemical_corrections(vib_e_list)

    return(vib)
    #__|


def thermochemical_corrections(vib_e_list, Temperature=300.0):
    """

    Args:
        vib_e_list:
            List of vibrational modes in eV
    """
    #| - thermochemical_corrections
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






#| - Atoms File Operations


def read_atoms_from_file():
    """
    Read atoms object from file. Checks several file names.
    """
    #| - read_atoms_from_file
    print("Reading atoms object from file")

    atoms = None

    file_name_list = [
        "init.traj",
        "init.POSCAR",
        "manual_start.traj",
        "out.traj",
        "POSCAR",
        ]

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
    """
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

#__|

#| - Atoms Geometry Methods ***************************************************


def angle_between_lattice_vectors(atoms, vector_0=0, vector_1=1):
    """
    Calculates the angle between the first 2 lattice vectors of a computational
    cell in degrees.
    """
    #| - angle_between_lattice_vectors
    v1 = atoms.cell[vector_0]
    v2 = atoms.cell[vector_1]

    angle = angle_between(v1, v2)
    angle = math.degrees(angle)

    return(angle)
    #__|


def magnitude_of_lattice_vectors(atoms):
    """
    Returns magnitude of three lattice vectors
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

#| - Modify Atoms Object


def move_atoms_of_element_i(atoms, element, new_position, dim="z"):
    """
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

    """
    """
    #| - displace_overlayer
    atoms = copy.deepcopy(atoms)

    x_frac = 1. * x_ind / mesh_size_x
    y_frac = 1. * y_ind / mesh_size_y

    # atoms = io.read("dir_atoms/init.traj")

    mag_lv0 = np.linalg.norm(atoms.cell[0])
    mag_lv1 = np.linalg.norm(atoms.cell[1])

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
    """
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

#__|

#| - Atoms Information Methods ************************************************


def number_of_atoms(atoms):
    """
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
    """
    """
    #| - number_of_constrained_atoms
    if type(atoms) == str:
        atoms = read(atoms)
    else:
        pass

    N_constraints = len(atoms.constraints)

    # print(N_constraints)
    return(N_constraints)
    #__|


def highest_position_of_element(atoms, element_symbol):
    """
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
    """Creates png images from an multi-atoms atoms object in the
    """
    #| - create_images_from_atoms_movies

    #| - SCRIPT PARAMETERS
    fold_name = "images"
    #__|


    #| - Reading Atoms Objects - OLD
    # default = 'qn.traj'
    # alist = {}
    #
    # if sys.argv[1:] == []:
    #     try:
    #         alist[default] = read(default)
    #     except:
    #         print 'Usage: python screenshot.py traj1 traj2 ...'
    # else:
    #     for arg in sys.argv[1:]:
    #         try:
    #             alist[arg] = read(arg, index=":")
    #         except:
    #             print 'Invalid traj file: ' + arg
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

            # print(path_i)
            # print(name_i)
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









#| - __old__


# def pdos_analysis():
#     """
#     OLD VERSION DON'T USE
#     """

#| - pdos_analysis
# from espresso import espresso
#
# a = read("qn.traj")
#
# calc = espresso(pw=500,             #plane-wave cutoff
#     dw=5000,            #density cutoff
#     xc="BEEF-vdW",      #exchange-correlation functional
#     kpts=(3,3,1),       # ark - k-points for hexagonal symmetry in 2-D materials
#     nbands=-20,         #20 extra bands besides the bands needed for valence electrons
#     spinpol = True,     # ark - added spinpolarizatoin
#     sigma=0.1,
#     psppath="/home/vossj/suncat/psp/gbrv1.5pbe",    #pseudopotential path
#     convergence= {
#         "energy": 1.e-5, #convergence parameters
#         "mixing": 0.1,
#         "nmix": 20,
#         "mix": 4,
#         "maxsteps": 500,
#         "diag": "david"
#         },
#     output = {"removesave": False},
#     outdir="calcdir"
#     )    #output directory for Quantum Espresso files
#
# # attach the espresso calculator to the surface
# a.set_calculator(calc)
# a.get_potential_energy()
#
# pdos=calc.calc_pdos(nscf=True,
#                     kpts=(6,6,1),
#                     Emin=-10.0,
#                     Emax=10.0,
#                     tetrahedra=True,
#                     sigma=0.2,
#                     DeltaE=0.01)
#
# pickle.dump(pdos,open("pdos.pkl", "w"))
#__|

#__|
