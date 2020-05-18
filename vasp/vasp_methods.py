# | - IMPORT MODULES
from raman_dft.vasp_raman_job_methods import get_modes_from_OUTCAR
# from raman_dft.vasp_raman_job_methods import parse_poscar

from ase_modules.ase_methods import create_gif_from_atoms_movies

from ase import io
import copy
import numpy as np
import os
import itertools
# __|

def num_of_atoms_OUTCAR_tmp(outcar_fh):
    """Parses OUTCAR for number of atoms in atoms object
    """
    # | - num_of_atoms_OUTCAR
    outcar_fh.seek(0)
    num_atoms = 0
    while True:
        line = outcar_fh.readline()
        if not line:
            break

        if "ion  position" in line:
            while True:
                next_line = outcar_fh.readline()

                if next_line == "\n":
                    break

                # ion_entry = next_line.split()[0].isdigit()
                num_atoms += 1
    return(num_atoms)
    # __|

def create_vib_modes_atoms(
    path_i=".",
    file_name="OUTCAR",
    modes_to_run="all",
    create_gifs=True,
    step_size=0.8,
    number_images=21,
    ):
    """Creates atoms movie for each vibrational mode.

    Args:
        path_i:
        file_name:
        modes_to_run:
        create_gifs:
        step_size:
        number_images:
    """
    # | - create_vib_modes_atoms

    # | - SCRIPT INPUTS
    # step_size = 1.2
    # number_images = 31
    vis_dir = "/mode_movies"
    # path_i = "/mnt/c/Users/raul_desktop/Dropbox/01_acad_folder
    # /01_grad_school/01_norskov/04_comp_clusters/00_scripts/
    # 09_raman_dft/view_modes"
    # __|

    # | - Reading in OUTCAR
    file_name = "OUTCAR"
    with open(path_i + "/" + file_name, "r") as fle:
        atoms = io.read(path_i + "/" + file_name)
        N_atoms = atoms.get_number_of_atoms()
        pos = atoms.positions

        eigvals, eigvecs, norms = get_modes_from_OUTCAR(fle, path_i=path_i)
    # __|

    # | - Creating Visualization Directory
    if not os.path.isdir(path_i + vis_dir):
        os.makedirs(path_i + vis_dir)
    # __|

    iterator = enumerate(itertools.izip(eigvals, eigvecs, norms))
    for index, (eigval_i, eigvec_i, norm_i) in iterator:

        # if not modes_to_run == "all":
        #     eigval_in_list = False
        #     for mode_to_run in modes_to_run:
        #         if eigval_i
        #         # if not eigval_i not in

        disps = np.linspace(-1, 1, num=number_images / 2)
        disps = np.append(disps, np.flip(disps, 0))

        master_pos_lst = []
        for disp_step in disps:
            pos_lst = []
            for k in range(N_atoms):

                p_k = pos[k]
                e_ik = eigvec_i[k]
                ss = step_size
                ds = disp_step
                l_lst = range(3)

                pos_disp = [p_k[l] + e_ik[l] * ss * ds / norm_i for l in l_lst]
                # pos_disp = [ pos[k][l] + eigvec_i[k][l] * step_size
                # * disp_step / norm_i for l in range(3)]
                pos_lst.append(pos_disp)

            atoms_i = copy.deepcopy(atoms)
            atoms_i.set_positions(pos_lst)

            master_pos_lst.append(atoms_i)

        eig_str = str(round(eigval_i, 3)).zfill(3)

        eig_str_lst = str(round(eigval_i, 3)).split(".")
        eig_str = eig_str_lst[0].zfill(4) + "." + eig_str_lst[1].ljust(3, "0")

        lead_num = str(index).zfill(2)
        folder_i = path_i + vis_dir + "/" + lead_num + "_" + eig_str + "_cm-1/"

        if not os.path.isdir(folder_i):
            os.makedirs(folder_i)


        out_traj_file = folder_i + "out-" + eig_str + ".traj"
        if os.path.isfile(out_traj_file):
            continue

        print("#####################################################")
        print("Creating: " + folder_i + "out-" + eig_str + ".traj")
        print("#####################################################")


        io.write(out_traj_file, master_pos_lst)

        if create_gifs:
            create_gif_from_atoms_movies(
                atoms_file="Default",
                path_i=folder_i,
                delay=10,
                )

    with open(".FINISHED", "w") as fle:
        fle.write("")
    # __|

def create_vdw_kernel_symlink():
    """
    If on the SLAC cluster, symlinks the vdw vasp kernel into the job directory,
    otherwise does nothing
    """
    # | - create_vdw_kernel_symlink
    if os.getenv("COMPENV") == "slac":
        print("TEMP - 180313 !@#")
        if not (os.path.exists("vdw_kernel.bindat")):
            os.symlink(
                "/nfs/slac/g/suncatfs/sw/vasp/vdw_kernel.bindat",
                "vdw_kernel.bindat",
                )


            # target = "/nfs/slac/g/suncatfs/sw/vasp/vdw_kernel.bindat"
            # tmpLink = "temp_vdw_kernel"
            # linkName = "vdw_kernel.bindat"
            # os.symlink(target, tmpLink)
            # os.rename(tmpLink, linkName)

    elif os.getenv("COMPENV") == "sherlock":
        pass

    if os.getenv("AWS_BATCH_JOB_ID") is None:
        pass

    # __|


def parse_incar(incar_list):
    """Manipulate INCAR data into python dictionary with correct data types.

    Args:
        incar_list:
            INCAR file in python list where each line represents a line from
            the file.
    """
    # | - parse_incar
    incar_1 = [line for line in incar_list if " = " in line]

    incar_dict = {}
    for line in incar_1:
        line_i = line.split("=")
        mess = "Each incar row should have 1 equals sign which "
        mess += "will be parsed into the LHS (key) and RHS (value)"
        assert len(line_i) == 2, mess
        incar_dict.update({line_i[0].strip(): line_i[1].strip()})

    # | - Incar keys list
    # incar_keys = [
    #     "ENCUT",
    #     "AMIX_MAG",
    #     "BMIX_MAG",
    #     "BMIX",
    #     "SIGMA",
    #     "AMIX",
    #     "EDIFF",
    #     "EDIFFG",
    #     "PREC",
    #     "GGA",
    #     "ALGO",
    #     "ISMEAR",
    #     "NPAR",
    #     "LDAUPRINT",
    #     "NELM",
    #     "IBRION",
    #     "IDIPOL",
    #     "ISIF",
    #     "ISPIN",
    #     "INIMIX",
    #     "NSW",
    #     "LORBIT",
    #     "LMAXMIX",
    #     "KPAR",
    #     "LDAUTYPE",
    #     "DIPOL",
    #     "LDAU",
    #     "LVTOT",
    #     "LDIPOL",
    #     "LASPH",
    #     "LREAL",
    #     "LDAUL",
    #     "LDAUU",
    #     "LDAUJ",
    #     ]
    # __|

    # | - Incar Types Dict

    incar_types_dict = {
        "ENCUT": "float",
        "AMIX_MAG": "float",
        "BMIX_MAG": "float",
        "BMIX": "float",
        "SIGMA": "float",
        "AMIX": "float",
        "EDIFF": "float",
        "EDIFFG": "float",

        # string
        "PREC": "string",
        "GGA": "string",
        "ALGO": "string",

        # float
        "ISMEAR": "integer",
        "NPAR": "integer",
        "LDAUPRINT": "integer",
        "NELM": "integer",
        "IBRION": "integer",
        "IDIPOL": "integer",
        "ISIF": "integer",
        "ISPIN": "integer",
        "INIMIX": "integer",
        "NSW": "integer",
        "LORBIT": "integer",
        "LMAXMIX": "integer",
        "KPAR": "integer",
        "LDAUTYPE": "integer",

        # [a, b, c]
        "DIPOL": "list",

        # True/False
        "LDAU": "boolean",
        "LVTOT": "boolean",
        "LDIPOL": "boolean",
        "LASPH": "boolean",

        # string
        "LREAL": "string",

        # [a, b]
        "LDAUL": "list",
        "LDAUU": "list",
        "LDAUJ": "list",
        }
    # __|

    # | - Formatting Dict to Proper Data Types
    formatted_incar_dict = {}
    for key, value in incar_dict.items():

        # if key == "LDIPOL":
        #     print(key)
        #     print(value)
        #     print("___;_-___--__")

        if key in incar_types_dict:

            data_type = incar_types_dict[key]

            if data_type == "float":
                value_new = float(value)
            elif data_type == "string":
                # Basically leave it alone
                value_new = value
            elif data_type == "integer":
                value_new = int(value)
            elif data_type == "list":
                # Figure this out later
                value_new = value
            elif data_type == "boolean":
                if value == ".FALSE.":
                    value_new = False
                elif value == ".TRUE.":
                    value_new = True
                else:
                    value_new = "ERROR 0847589347"

            else:
                value_new = value
        else:
            value_new = value

        formatted_incar_dict.update({key: value_new})
    # __|

    return(formatted_incar_dict)

    # return(incar_dict)
    # __|
