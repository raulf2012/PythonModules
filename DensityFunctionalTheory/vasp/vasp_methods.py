"""Vasp methods."""

# | - Import Modules
import os
import sys

import copy
import itertools
from pathlib import Path

import numpy as np
from ase import io

# #########################################################
#  from raman_dft.vasp_raman_job_methods import get_modes_from_OUTCAR
# from raman_dft.vasp_raman_job_methods import parse_poscar
from ase_modules.ase_methods import create_gif_from_atoms_movies
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


def read_incar(path_i, verbose=False):
    """
    """
    #| - read_incar
    incar_path = os.path.join(
        path_i,
        "INCAR")
    my_file = Path(incar_path)
    if my_file.is_file():
        with open(incar_path, "r") as f:
            incar_lines = f.read().splitlines()

        incar_dict = parse_incar(incar_lines)
    else:
        incar_dict = None
        if verbose:
            print("Couldn't get INCAR dict, maybe is not at path specified")
            print("path_i", "\n", path_i)

        # nsw_i = incar_dict["NSW"]
        # nelm_i = incar_dict["NELM"]
        # incar_parsed = True

    return(incar_dict)
    #__|


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


# from read_bad_outcar import read_vasp_out
from vasp.read_bad_outcar import read_vasp_out

def read_bad_OUTCAR(filename):
    """
    """
    #| - read_bad_OUTCAR
    traj = read_vasp_out(
        # filename="OUTCAR",
        # "OUTCAR",
        filename,
        index=":",
        )

    return(traj)
    #__|


def get_irr_kpts_from_outcar(path_i):
    """Get irreducible number of k-points from OUTCAR file.

    Printed towards beginning of file before SCF cycles begin.
    """
    #| - get_irr_kpts_from_outcar
    outcar_path = os.path.join(path_i, "OUTCAR")

    N = 1000

    my_file = Path(outcar_path)
    if my_file.is_file():
        num_lines = sum(1 for line in open(outcar_path))
        if num_lines < N:
            N = num_lines - 1

        with open(outcar_path) as myfile:
            outcar_lines = [next(myfile) for x in range(N)]
    else:
        return(None)


    # Parse lines looking for irr. kpts
    line_tmp = None
    for line_i in outcar_lines:
        frag_i = "irreducible k-points"
        if frag_i in line_i:
            # print(line_i)
            line_tmp = line_i
            break

    irr_kpts = None
    if line_tmp is not None:
        # Extract number of irr. kpts from line
        line_split = [i for i in line_tmp.split(" ") if i != ""]

        irr_kpts = None
        found_found = False
        for i in line_split:

            # I know that the number of irreducible k-points comes after the "Found"
            if i == "Found":
                found_found = True

            if i.isnumeric() and found_found:
                irr_kpts = int(i)

    return(irr_kpts)
    #__|


def read_ase_sort_dat(path_i=None):
    """
    """
    #| - read_ase_sort_dat
    # path_i = os.path.join(
    #     gdrive_path_i,
    #     "ase-sort.dat")

    with open(path_i, "r") as f:
        lines = f.read().splitlines()


    # REMOVE
    # atom_index_mapping = dict()

    sort_list = []
    resort_list = []
    for line_i in lines:
        line_i = [i for i in line_i.split(" ") if i != ""]
        line_i = [int(i) for i in line_i]

        sort_list.append(line_i[0])
        resort_list.append(line_i[1])

        # REMOVE
        # atom_index_mapping[line_i[0]] = line_i[1]

    atom_index_mapping = dict(zip(
        resort_list,
        list(range(len(resort_list))),
        ))

    return(
        atom_index_mapping,
        sort_list,
        resort_list,
        )
    #__|


def parse_dipole_OUTCAR(file_path):
    """Parse dipole data from OUTCAR generated from LCALCPOL=True run.
    """
    # | - parse_dipole_OUTCAR
    file_i = file_path

    lookup = 'Ionic dipole moment'
    with open(file_i) as f:
        for num, line in enumerate(f, 1):
            if lookup in line:
                start_line = num

    file_lines = []
    with open(file_i) as f:
        for i, line in enumerate(f):
            if i > start_line - 2:
                file_lines.append(line)
            if i > start_line + 6:
                break

    file_lines_2 = []
    for line in file_lines:

        # | - Parse for ionic dipole moment
        phrase = 'Ionic dipole moment'
        if phrase in line:
            ion_line = line

            match_ind = ion_line.index(phrase)
            tmp0 = ion_line[match_ind + len(phrase):]

            phrase1 = 'p[ion]=('
            match_ind = tmp0.index(phrase1)
            tmp1 = tmp0[match_ind + len(phrase1):]

            phrase2 = ') electrons'
            match_ind = tmp1.index(phrase2)
            tmp2 = tmp1[:match_ind]

            ion_dipole_mom = [float(i) for i in tmp2.split(' ') if i != '']
        # __|

        # | - Parse for electronic dipole moment
        phrase = 'Total electronic dipole moment'
        if phrase in line:
            elec_line = line

            match_ind = elec_line.index(phrase)
            tmp0 = elec_line[match_ind + len(phrase):]

            phrase1 = 'p[elc]=('
            match_ind = tmp0.index(phrase1)
            tmp1 = tmp0[match_ind + len(phrase1):]

            phrase2 = ') electrons'
            match_ind = tmp1.index(phrase2)
            tmp2 = tmp1[:match_ind]

            elc_dipole_mom = [float(i) for i in tmp2.split(' ') if i != '']
        # __|

    out_dict = dict()
    out_dict['elc_dipole_mom'] = elc_dipole_mom
    out_dict['ion_dipole_mom'] = ion_dipole_mom
    return(out_dict)
    # __|


def get_BEC_data(outcar_path):
    """Parse Born Effective Charge data from OUTCAR (LEPSILON calculation)

    Note: There is also a pymatgen instance for this

        from pymatgen.io.vasp.outputs import Outcar
        OUTCAR_INSTANCE = Outcar(outcar_path)
        OUTCAR_INSTANCE.read_lepsilon()
        OUTCAR_data = OUTCAR_INSTANCE.as_dict()
        BEC_data = OUTCAR_data['born']

    """
    # | - get_BEC_data

    # Read in OUTCAR  lines

    with open(outcar_path, 'r') as f:
        outcar_lines = f.read().splitlines()

    # Read OUTCAR as ASE atoms object
    # atoms = io.read(outcar_path)

    try:
        atoms = io.read(outcar_path)
    except:
        # Construct POSCAR path, use to read ase atoms object instead of OUTCAR
        poscar_path = os.path.join(
            '/'.join(outcar_path.split('/')[0:-1]),
            'POSCAR')
        atoms = io.read(poscar_path)

    num_atoms = atoms.arrays['numbers'].shape[0]

    start_index = None
    BEC_outcar_str = 'BORN EFFECTIVE CHARGES (including local field effects) (in e, cummulative output)'
    for index, line in enumerate(outcar_lines):
        if BEC_outcar_str in line:
            start_index = index

    # Relevent lines from OUTCAR containing BEC data
    outcar_lines_2 = outcar_lines[
        start_index + 2:
        start_index + 4 * num_atoms + 2
        ]

    # #####################################################
    def chunks(lst, n):
        for i in range(0, len(lst), n):
            yield lst[i:i + n]

    correct_atom_indices = False

    BEC_data = dict()
    for index, i in enumerate(chunks(outcar_lines_2, 4)):

        atom_index = int(i[0].split(' ')[-1])
        if index == 0 and atom_index == 1:
            correct_atom_indices = True
        if correct_atom_indices:
            atom_index = atom_index - 1

        BEC_tensor = []
        for j in i[1:]:
            BEC_j = [float(i) for i in [i for i in j.split(' ') if i != ''][1:]]
            BEC_tensor.append(BEC_j)
        BEC_data[atom_index] = BEC_tensor

    # for key, val in BEC_data.items():
    #     print(44 * '-')
    #     print(np.array(val))

    return(BEC_data)
    # __|
