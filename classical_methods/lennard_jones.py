#!/usr/bin/env python

"""Wrapper code around ASAP & ASE LJ implementation.

Author: Raul A. Flores
"""

#| - Import Modules
import gc
import math
import numpy as np  # import pandas as pd
import collections
from pymatgen.core.periodic_table import Element
from asap3 import LennardJones
#__|

def lennard_jones_sp(
    epsilon,
    sigma,
    atoms,
    modified_lj=True,
    normalize_per_atom=True,
    ):
    """Calculate single-point energy and forces with Lennard Jones force field.

    Because the order of the elements must be internally consistent, the
    convention will be that the element list goes from smalles to largest
    atomic number.

    Args:
        epsilon:
        sigma:
        atoms:
        modified_lj:
    """
    #| - lennard_jones_sp
    atomic_num_list = atoms.get_atomic_numbers()
    atomic_num_list = list(set(atomic_num_list))
    atomic_num_list.sort()

    atomic_type_num_dict = collections.Counter(atomic_num_list)

    orig_num_of_atoms = atoms.get_number_of_atoms()

    #| - Filter Relevant LJ Parameters
    row_col_to_keep = [
        Element.from_Z(atomic_num).name
        for
        atomic_num
        in
        atomic_num_list
        ]

    epsilon = epsilon.loc[row_col_to_keep]
    epsilon = epsilon[row_col_to_keep]

    sigma = sigma.loc[row_col_to_keep]
    sigma = sigma[row_col_to_keep]

    epsilon = epsilon.values
    sigma = sigma.values
    #__|

    # TEMP_PRINT
    # print(epsilon)
    # print(sigma)

    calc = LennardJones(
        list(atomic_type_num_dict),
        epsilon,
        sigma,
        rCut=-1,
        modified=modified_lj,
        )

    #| - Repeat Unit Cell
    def calc_cell_heights(unit_cell):
        """Calculate heights of cell.

        Obtained code from ASAP NormalAtoms.invert_cell

        Args:
            unit_cell:
        """
        #| - calc_cell_heights
        determinant = np.cross(
            unit_cell[0],
            unit_cell[1],
            )
        determinant = abs(determinant.dot(unit_cell[2]))

        heights = []
        for i_ind, unit_vect_i in enumerate(unit_cell):
            inv = np.cross(
                unit_cell[(i_ind + 1) % 3],
                unit_cell[(i_ind + 2) % 3],
                )

            den = math.sqrt(np.dot(inv, inv))
            height_i = determinant / den
            heights.append(height_i)

        return(heights)
        #__|

    heights = calc_cell_heights(atoms.cell)

    # The cut-off length is internally set by ASAP to 3 * 2 * (max_sigma_value)
    max_sigma = sigma.flatten().max()
    cut_off_length = 3 * 2 * max_sigma

    #| - Repeat Unit Cell
    repeat_unit_cell = []
    for i_ind, height_i in enumerate(heights):
        if height_i < cut_off_length:
            cell_repeat_fact = math.ceil(cut_off_length / height_i)
            cell_repeat_fact = int(cell_repeat_fact)
            repeat_unit_cell.append(cell_repeat_fact)
        else:
            repeat_unit_cell.append(1)
    #__|

    # print("KJIJIJKJUKJIJ")
    # print(repeat_unit_cell)

    atoms = atoms.repeat(repeat_unit_cell)
    #__|

    atoms.set_calculator(calc)

    lj_energy = atoms.get_potential_energy()

    # print("Unnormalized LJ energy:")
    # print(lj_energy)

    lj_energy_per_atom = lj_energy / atoms.get_number_of_atoms()

    # This is total energy w.r.t. the original number of atoms in the
    # computational cell (since the cell was repeated)
    lj_total_energy = lj_energy_per_atom * orig_num_of_atoms

    if normalize_per_atom:
        lj_energy = lj_energy_per_atom

    else:
        lj_energy = lj_total_energy

    # print("lj_energy - lennard_jones.py")
    # print(lj_energy)

    # TEMP
    # print("TMP - running del atoms and gc.collect()")
    # del atoms
    # del calc
    # gc.collect()

    return(lj_energy)
    #__|



# def calc_lennard_jones_formation_energy(
#     atoms,
#
#     ):
#     """
#     """
#     #| - calc_lennard_jones_formation_energy
#
#
#     #__|
