#| - IMPORT MODULES
import numpy as np
import sys
import os
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

    nx,ny,nz = np.shape(cd)
    #cut away periodic image planes to correct QE output
    u=nx-1
    v=ny-1
    w=nz-1
    cd2 = np.empty((u,v,w), np.float)
    for i in range(u):
        for j in range(v):
            cd2[i][j][:] = cd[i][j][:w]

    write(file_name, atoms, data=cd2)

    #edit density.cube grid size if odd number of grid points to correct for old versions of ASE
    bohr = 0.52917721092
    cell = atoms.get_cell()
    da = cell[0]/(u*bohr)
    db = cell[1]/(v*bohr)
    dc = cell[2]/(w*bohr)

    f = open(file_name,"r")
    lines = f.readlines()
    f.close()

    line3 = "%5.0f    %.6f    %.6f    %.6f\n"%(u,da[0],da[1],da[2])
    line4 = "%5.0f    %.6f    %.6f    %.6f\n"%(v,db[0],db[1],db[2])
    line5 = "%5.0f    %.6f    %.6f    %.6f\n"%(w,dc[0],dc[1],dc[2])

    lines[3] = line3
    lines[4] = line4
    lines[5] = line5

    f = open(file_name,"w")
    f.writelines(lines)
    f.close()
    #__|

def cleanup(suffix="", save_cube=False):
    """
    """
    #| - cleanup
    if not os.path.exists("dir_bader"):
        os.makedirs("dir_bader")

    if not save_cube:
        os.system("rm density" + suffix + ".cube")

    os.system("mv ACF.dat dir_bader/.ACF%s.dat"%suffix)
    os.system("rm AVF.dat")
    os.system("rm BCF.dat")
    os.system("mv bader.out dir_bader/.bader%s.out"%suffix)
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
            atoms[i].charge -= float(line[4]) - atoms.calc.get_nvalence()[1][atoms[i].symbol]

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
        for i,line in enumerate(lines[2:-4]):
            line = line.split()
            atoms[i].charge = atoms.calc.get_nvalence()[1][atoms[i].symbol] - float(line[4])
            charge_list.append(atoms[i].charge)

        atoms.info.update({"bader_charges": charge_list})
        cleanup()
    #__|

    #__|

def bader(atoms, spinpol=False, outdir=None, run_exec=True):
    """
    Calculate charge density using atoms.calc, calculate bader charges on each
    atom in atoms, and assign to atom.data["bader_charge"].
    If spinpol: also assign atom.data["bader_magmom"].
    """
    #| - bader
    print("Executing Badder analysis"); sys.stdout.flush()
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
        os.system("rm %s/charge.log"%outdir)

    atoms.set_calculator(calc=calc)
    atoms.write("out.traj")
    #__|
