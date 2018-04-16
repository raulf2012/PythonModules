#!/usr/bin/env python

"""Quantum Espresso methods and code.

Author: Raul A. Flores
"""
#| - Import Modules
import os
import pandas as pd
#__|

#| - Log File Methods

def tot_abs_magnetization(log="log"):
    """Return total and absolute magnetization vs SCF iteration.

    Abs:
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

def element_index_dict(log="log"):
    """Returns index: element dictionary.

    Format: {0: 'Fe', 1: 'Fe', 2: 'Fe'}

    Args:
        log
    """
    #| - element_index_dict
    elem_ind_dict = {}
    with open(log, "r") as fle:

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
                        line_list = [i_cnt for i_cnt in line_list if i_cnt != ""]

                        ind_i = int(line_list[0]) - 1  # "0" indexed

                        elem_i = line_list[1]
                        elem_i = ''.join([i for i in elem_i if not i.isdigit()])

                        elem_ind_dict[ind_i] = elem_i

                    else:
                        break

            #__|

    return(elem_ind_dict)
    #__|

def magmom_charge_data(log="log"):
    """Returns charge and magmom data per atom for all SCF iterations.
    # [{'magmom': 2.7931, 'charge': 7.3942, 'atom_num': 0},
    # {'magmom': 2.7931, 'charge': 7.3942, 'atom_num': 1}, ...],
    # ]

    Args:
        log
    """
    #| - magmom_charge_data

    #| - Reading Log File
    fle = open(log, "r")

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
    elem_ind_dict = element_index_dict(log=log)

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

def scf_convergence(log="log"):
    """Return SCF convergence vs iteration.

    Author: ???
    I didn't write this.

    Args:
        log
    """
    #| - scf_convergence
    rydberg = 13.6057  # rydberg to eV conversion

    filename = "scf_temp.txt"
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
    items = filter(None, file.read().split(" "))
    pw = float(items[3]) * 13.606
    os.system("rm %s" % filename)

    return(scf, iter, pw)
    #__|

#__|
