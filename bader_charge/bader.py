#!/usr/bin/env python

"""Bader charge analysis methods.

Author(s): Colin Dickins wrote most of this; Raul A. Flores
"""

#| - IMPORT MODULES
import sys
import os

import numpy as np
from ase.io import write
#__|

def cd2cube(atoms, spin=""):
    """
    Takes charge density pickle file from ase-qe and writes a cube file ready
    for bader analysis to the current directory
    """
    #| - cd2cube
    # cd2cube(atoms.calc.extract_charge_density(spin="up")[2], atoms)

    if spin == "":
        cd = atoms.calc.extract_charge_density()[2]
    else:
        cd = atoms.calc.extract_charge_density(spin=spin)[2]

    file_name = "density" + spin + ".cube"

    nx, ny, nz = np.shape(cd)
    #cut away periodic image planes to correct QE output
    u = nx - 1
    v = ny - 1
    w = nz - 1
    cd2 = np.empty((u, v, w), np.float)
    for i in range(u):
        for j in range(v):
            cd2[i][j][:] = cd[i][j][:w]

    write(file_name, atoms, data=cd2)

    # edit density.cube grid size if odd number of grid points to
    # correct for old versions of ASE
    bohr = 0.52917721092
    cell = atoms.get_cell()
    da = cell[0] / (u * bohr)
    db = cell[1] / (v * bohr)
    dc = cell[2] / (w * bohr)

    f = open(file_name, "r")
    lines = f.readlines()
    f.close()

    line3 = "%5.0f    %.6f    %.6f    %.6f\n" % (u, da[0], da[1], da[2])
    line4 = "%5.0f    %.6f    %.6f    %.6f\n" % (v, db[0], db[1], db[2])
    line5 = "%5.0f    %.6f    %.6f    %.6f\n" % (w, dc[0], dc[1], dc[2])

    lines[3] = line3
    lines[4] = line4
    lines[5] = line5

    f = open(file_name, "w")
    f.writelines(lines)
    f.close()
    #__|

def cleanup(suffix="", save_cube=True):
    """
    """
    #| - cleanup
    if not os.path.exists("dir_bader"):
        os.makedirs("dir_bader")

    if not save_cube:
        os.system("rm density" + suffix + ".cube")
    else:
        os.system("mv density" + suffix + ".cube dir_bader")

    os.system("mv ACF.dat dir_bader/ACF%s.dat" % suffix)
    # os.system("rm AVF.dat")
    # os.system("rm BCF.dat")
    os.system("mv bader.out dir_bader/bader%s.out" % suffix)

    os.system("mv AVF.dat dir_bader/AVF.dat")
    os.system("mv BCF.dat dir_bader/BCF.dat")
    #__|

def bader_exec(atoms, spin=""):
    """
    Run bader executable on cube density file
    """
    #| - bader_exec
    bash_comm = "bader density" + spin + ".cube >> bader.out"
    os.system(bash_comm)

    #| - Spin Polarlized Calculation
    if spin == "up":

        f = open("ACF.dat"); lines = f.readlines(); f.close()
        for i, line in enumerate(lines[2:-4]):
            line = line.split()
            atoms[i].magmom = float(line[4])
            atoms[i].charge = -float(line[4])
        cleanup(suffix=spin)

    elif spin == "down":
        f = open("ACF.dat"); lines = f.readlines(); f.close()

        magmom_list = []
        charge_list = []
        for i, line in enumerate(lines[2:-4]):
            line = line.split()
            atoms[i].magmom -= float(line[4])
            val_i = atoms.calc.get_nvalence()[1][atoms[i].symbol]
            atoms[i].charge -= float(line[4]) - val_i

            magmom_list.append(atoms[i].magmom)
            charge_list.append(atoms[i].charge)

        atoms.info.update({"magmom_set": True})
        atoms.info.update({"bader_magmoms": magmom_list})
        atoms.info.update({"bader_charges": charge_list})

        cleanup(suffix=spin)
    #__|

    #| - Non-Spin Polarized Calculation
    elif spin == "":
        charge_list = []
        f = open("ACF.dat"); lines = f.readlines(); f.close()
        for i, line in enumerate(lines[2:-4]):
            line = line.split()

            charge_i = atoms.calc.get_nvalence()[1][atoms[i].symbol]
            atoms[i].charge = charge_i - float(line[4])
            charge_list.append(atoms[i].charge)

        atoms.info.update({"bader_charges": charge_list})
        cleanup()
    #__|

    #__|

def bader(atoms, spinpol=False, outdir=None, run_exec=True):
    """Perform bader charge analysis on atoms.

    Calculate charge density using atoms.calc, calculate bader charges on each
    atom in atoms, and assign to atom.data["bader_charge"].

    If spinpol: also assign atom.data["bader_magmom"].

    Args:
        atoms: ASE atoms object
        spinpol: Spin polarized calculation
        outdir: Output directory location
        run_exec: Whether to run bader executable or just create preliminary
        file (some clusters don't/can't have the bader fortran code)
    """
    #| - bader
    mess = "Executing Bader Analysis "
    mess += "*****************************************************"
    print(mess); sys.stdout.flush()

    #| - Don't Run Bader Executable on AWS
    if "COMPENV" not in os.environ:
        print("COMPENV env. var. doesn't exits, probably in AWS?")
        print("Bader executable turned off")

        run_exec = False
    else:
        pass
    #__|

    calc = atoms.calc

    #| - Using Spin Polarization
    if spinpol:

        #| - Spin up
        cd2cube(atoms, spin="up")
        if run_exec:
            bader_exec(atoms, spin="up")
        #__|

        #| - Spin down
        cd2cube(atoms, spin="down")
        if run_exec:
            bader_exec(atoms, spin="down")
        #__|

        print("BADER MAGMOMS: " + str(atoms.get_initial_magnetic_moments()))

    #__|

    #| - Not Spin Polarized
    else:
        cd2cube(atoms)

        if run_exec:
            bader_exec(atoms)
    #__|

    print("BADER CHARGES: " + str(atoms.get_initial_charges()))

    # if run_exec:
    #     write("dir_bader/bader.traj", atoms)

    if outdir:
        os.system("rm %s/charge.log" % outdir)

    atoms.set_calculator(calc=calc)
    atoms.write("out.traj")
    #__|
