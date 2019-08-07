
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
    return_quantity="energy",  # 'energy', 'forces', 'both', 'energy&forces'
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
        return:
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

    # epsilon = epsilon.values
    # sigma = sigma.values
    #__|

    calc = LennardJones(
        list(atomic_type_num_dict),
        epsilon.values,
        sigma.values,
        rCut=-1,
        modified=modified_lj,
        )

    #| - Repeat Unit Cell
    repeat_unit_cell = repeat_unit_cell_ASAP(atoms, sigma)

    atoms = atoms.repeat(repeat_unit_cell)
    #__|

    atoms.set_calculator(calc)

    #| - Energy
    lj_energy = atoms.get_potential_energy()
    lj_energy_per_atom = lj_energy / atoms.get_number_of_atoms()

    # This is total energy w.r.t. the original number of atoms in the
    # computational cell (since the cell was repeated)
    lj_total_energy = lj_energy_per_atom * orig_num_of_atoms

    if normalize_per_atom:
        lj_energy = lj_energy_per_atom

    else:
        lj_energy = lj_total_energy
    #__|

    #| - Forces -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    # COMBAK
    # Find a way to decrease the number of atoms again
    lj_forces = atoms.get_forces()

    #__|

    if return_quantity == "energy":
        out = lj_energy
    elif return_quantity == "forces":
        out = lj_forces
    elif return_quantity == "both" or return_quantity == "both":
        out = (lj_energy, lj_forces)

    return(out)
    #__|

def repeat_unit_cell_ASAP(atoms, sigma):
    """Return repeat array such that ASAP doesn't complain about r_cut.

    Args:
        atoms
    """
    #| - repeat_unit_cell_ASAP

    def calc_cell_heights(unit_cell):
        """Calculate heights of cell.

        Obtained code from ASAP NormalAtoms.invert_cell

        Args:
            unit_cell:
            sigma:
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

    # max_sigma = sigma.flatten().max()
    max_sigma = sigma.values.flatten().max()

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

    return(repeat_unit_cell)
    #__|
