#!/usr/bin/env python

"""Wrapper code around ASAP & ASE LJ implementation.

Author: Raul A. Flores
"""

#| - Import Modules
# import os
# import copy
# import pickle
from scipy import stats
from scipy.optimize import fmin

import math

import numpy as np
import pandas as pd

import collections

import plotly.plotly as py
import plotly.graph_objs as go

from ase.visualize import view
from ase import io
from asap3 import LennardJones
#__|


def lennard_jones_sp(epsilon, sigma, atoms):
    """Calculate single-point energy and forces with Lennard Jones force field.

    Args:
        pars:
        atoms:
    """
    #| - lennard_jones_sp

    #| - __old__
    # ep_11 = pars[0]
    # ep_22 = pars[1]
    # ep_12 = pars[2]
    # si_11 = pars[3]
    # si_22 = pars[4]
    # si_12 = pars[5]
    #
    # epsilon = [
    #     ep_11, 0.0,
    #     ep_12, ep_22,
    #     ]
    #
    # sigma = [
    #     si_11, 0.0,
    #     si_12, si_22,
    #     ]
    #__|

    atomic_num_list = atoms.get_atomic_numbers()
    atomic_type_num_dict = collections.Counter(atomic_num_list)



    #| - Filter Relevant LJ Parameters
    # row_col_to_keep =

    epsilon = epsilon.loc[row_col_to_keep]
    epsilon = epsilon[row_col_to_keep]

    sigma = sigma.loc[row_col_to_keep]
    sigma = sigma[row_col_to_keep]

    #__|


    calc = LennardJones(
        list(atomic_type_num_dict),
        epsilon,
        sigma,
        )

    #| - Repeat Unit Cell
    # Repeat until cell so that ASAP doesn't freak out
    cut_off_length = 6

    repeat_unit_cell = []
    for lattice_vect_i in atoms.cell:
        lattice_vect_i_magn = np.linalg.norm(lattice_vect_i)

        if lattice_vect_i_magn < cut_off_length:
            cell_repeat_fact = math.ceil(cut_off_length / lattice_vect_i_magn)
            cell_repeat_fact = int(cell_repeat_fact)
            repeat_unit_cell.append(cell_repeat_fact)
        else:
            repeat_unit_cell.append(1)

    atoms = atoms.repeat(repeat_unit_cell)
    #__|

    atoms.set_calculator(calc)
    lj_energy = atoms.get_potential_energy()

    lj_energy = lj_energy / atoms.get_number_of_atoms()

    #| - __old__
    # atoms_list = [atoms_i.repeat((3,3,3)) for atoms_i in atoms_list]
    # [atoms_i.set_calculator(calc) for atoms_i in atoms_list]
    # predicted_energies = np.array([atoms.get_potential_energy() for atoms in atoms_list])
    # return(predicted_energies)
    #__|

    return(lj_energy)
    #__|
