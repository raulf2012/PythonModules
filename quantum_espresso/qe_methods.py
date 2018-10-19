#!/usr/bin/env python

"""Quantum Espresso methods and code.

Author: Raul A. Flores
"""

#| - Import Modules
import os
import pandas as pd

import numpy as np
#__|

#| - Log File Methods

def number_of_atoms(path_i=".", log="log"):
    """Return number of atoms from QE log file.

    Args:
        path_i:
        log:
    """
    #| - number_of_atoms
    file_name = path_i + "/" + log
    with open(file_name, "r") as fle:
        fle.seek(0)  # just in case

        while True:
            line = fle.readline()

            if "Cartesian axes" in line:
                fle.readline().strip()  # Blank line
                fle.readline()  # Column headers

                atom_list = []
                while True:
                    data_line_i = fle.readline().strip()

                    if data_line_i == "":
                        break

                    atom_list.append(data_line_i)

            if not line:
                break

    num_atoms = len(atom_list)

    return(num_atoms)
    #__|

def tot_abs_magnetization(path_i=".", log="log"):
    """Return total and absolute magnetization vs SCF iteration.

    Abs:
        path_i
        log
    """
    #| - tot_abs_magnetization
    fle = open(log, "r")
    fle.seek(0)  # just in case

    tot_mag_list = []
    abs_mag_list = []
    while True:
        line = fle.readline()

        # Break out of loop
        if not line:
            break

        #| - Searching for Atomic Magmoms
        if "total magnetization" in line:
            line_list = line.strip().split(" ")
            line = [i for i in line_list if i != ""]
            tot_mag = float(line[3])
            tot_mag_list.append(tot_mag)

        if "absolute magnetization" in line:
            line_list = line.strip().split(" ")
            line = [i for i in line_list if i != ""]
            abs_mag = float(line[3])
            abs_mag_list.append(abs_mag)

        #__|

    df_tot = pd.DataFrame(tot_mag_list, columns=["tot_mag"])
    df_abs = pd.DataFrame(abs_mag_list, columns=["abs_mag"])

    df = pd.concat([df_tot, df_abs], axis=1)

    return(df)
    #__|

def element_index_dict(path_i=".", log="log"):
    """Return index: element dictionary.

    Format: {0: 'Fe', 1: 'Fe', 2: 'Fe'}

    Args:
        path_i
        log
    """
    #| - element_index_dict
    elem_ind_dict = {}
    file_name = path_i + "/" + log
    with open(file_name, "r") as fle:

        fle.seek(0)  # just in case
        while True:
            line = fle.readline()

            # Break out of loop
            if not line:
                break

            #| - Atom Index <---> Atom Type
            if "Cartesian axes" in line:
                fle.readline()  # Blank line
                fle.readline()  # Column header line

                # elem_ind_dict = {}
                while True:
                    line_i = fle.readline()
                    if "tau(" in line_i:

                        line_list = line_i.strip().split(" ")
                        line_list = [i_ct for i_ct in line_list if i_ct != ""]

                        ind_i = int(line_list[0]) - 1  # "0" indexed

                        elem_i = line_list[1]
                        elem_i = ''.join([i for i in elem_i if not i.isdigit()])

                        elem_ind_dict[ind_i] = elem_i

                    else:
                        break

            #__|

    return(elem_ind_dict)
    #__|

def magmom_charge_data(path_i=".", log="log"):
    """Return charge and magmom data per atom for all SCF iterations.

    Args:
        path_i
        log
    """
    #| - magmom_charge_data

    #| - Reading Log File
    file_name = path_i + "/" + log

    fle = open(file_name, "r")
    fle.seek(0)  # just in case

    master_list = []
    while True:
        line = fle.readline()

        # Break out of loop
        if not line:
            break

        #| - Searching for Atomic Magmoms
        if "Magnetic moment per site" in line:

            list_i = []
            while True:
                # print("ksjfksjfks - 3")
                line_i = fle.readline()
                if "atom:" in line_i:

                    line_list = line_i.strip().split(" ")
                    line_list = [i_cnt for i_cnt in line_list if i_cnt != ""]

                    atom_num = int(line_list[1]) - 1  # "0" indexed
                    charge = float(line_list[3])
                    magmom = float(line_list[5])

                    entry_dict = {
                        "atom_num": atom_num,
                        "charge": charge,
                        "magmom": magmom,
                        }

                    list_i.append(entry_dict)
                else:
                    break

            master_list.append(list_i)
        #__|

    #__|

    #| - Creating Pandas DataFrame
    if master_list == []:
        print("Magmom/charge data not found, ",
            "calculation probably not spin-polarized"
            )
        return(None)

    elem_ind_dict = element_index_dict(path_i=path_i, log=log)

    df_list = []
    for i_cnt, iter_i in enumerate(master_list):
        df_i = pd.DataFrame(iter_i)

        df_list.append(df_i)
        df_i["iteration"] = i_cnt

    # Combining data frames
    df = pd.concat(df_list)

    # Creating element type column (mapping index <--> element)
    df["element"] = df["atom_num"].map(elem_ind_dict)

    # Resetting index
    df.reset_index()
    #__|

    return(df)
    #__|

def estimate_magmom(
    path_i=".",
    atoms=None,
    log="log",
    ):
    """Estimage magnetic moments from QE log file.

    TODO Also read charges data for last iteration
    TODO Modify to not need atoms object to function

    Author: Colin Dickens

    Estimate magmom from log file (based on charge spheres centered on atoms)
    and assign to atoms object to assist with calculation restart upon
    unexpected interruption.
    """
    #| - estimate_magmom
    num_atoms = number_of_atoms(path_i=path_i, log=log)
    # print(num_atoms)

    file_name = path_i + "/" + log

    with open(file_name, "r") as fle:
        lines = fle.readlines()

    i = len(lines) - 1
    while True:

        #| - If magmom/charge data is not found in log file
        # The calculation is probably not spin-polarized
        if i == 0:
            print("estimate_magmom - Could not find magmom/charge data \n",
                "Calculation is probably not spin polarized"
                )
            return(None)
            break
            # raise IOError("Could not identify espresso magmoms")

        #__|

        line = lines[i].split()
        if len(line) > 3:
            if line[0] == "absolute":
                abs_magmom = float(line[3])
        if len(line) > 6:
            if line[4] == "magn:":
                i -= num_atoms - 1
                break
        i -= 1

    magmom_list = []
    charge_list = []

    if abs_magmom < 1e-3:
        print("estimate_magmom | Absolute magnetism is near 0, setting "
            "initial atomic magmoms to 0"
            )

        # for atom in atoms:
        #     atom.magmom = 0

        for atom in range(num_atoms):
            magmom_list.append(0.)

    else:
        total_esp_magmom = 0
        # for j in range(len(atoms)):
        for j in range(num_atoms):
            total_esp_magmom += np.abs(float(lines[i + j].split()[5]))

        magmom_list = []
        charge_list = []
        # for j in range(len(atoms)):
        for j in range(num_atoms):
            charge_i = float(lines[i + j].split()[3])
            magmom_i = float(lines[i + j].split()[5])
            new_magmom_i = magmom_i * abs_magmom / total_esp_magmom
            # atoms[j].magmom = new_magmom_i
            magmom_list.append(new_magmom_i)
            charge_list.append(charge_i)

    magmom_list = np.array(magmom_list)
    charge_list = np.array(charge_list)

    if atoms is not None:
        atoms.info.update({"qe_log_magmoms": magmom_list})
        atoms.info.update({"qe_log_charges": charge_list})

    return(magmom_list, charge_list)
    #__|

def scf_convergence(path_i=".", log="log"):
    """Return SCF convergence vs iteration.

    Author: ???
    I didn't write this.

    TODO: Not working well with python3.6, probably best to rewrite this

    Args:
        log
    """
    #| - scf_convergence
    rydberg = 13.6057  # rydberg to eV conversion

    filename = os.path.join(
        os.environ["HOME"],
        "scf_temp.txt",
        )

    os.system("grep scf %s > %s" % (log, filename))

    file = open(filename)
    file.seek(0)  # just in case

    lines = file.read().split("\n")
    lines = filter(None, lines)

    scf = []
    iter = []
    n = 0
    for line in lines:
        items = filter(None, line.split(" "))
        try:
            scf.append(float(items[4]) * rydberg)
            n += 1
            iter.append(n)
        except:
            continue

    os.system("grep cpu log > %s" % filename)
    file = open(filename)
    lines = file.read().split("\n")
    lines = filter(None, lines)
    #time = float(lines[-1][40:50])/3600/48
    file.close()

    os.system('grep "kinetic-energy cutoff" %s > %s' % (log, filename))
    file = open(filename)

    # FIXME This is breaking in python3.6
    # items = filter(None, file.read().split(" "))
    items = list(filter(None, file.read().split(" ")))

    pw = float(items[3]) * 13.606
    os.system("rm %s" % filename)

    return(scf, iter, pw)
    #__|

#__|
