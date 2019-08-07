#!/usr/bin/env python

"""General method to calculate formation energies.

Author: Raul A. Flores
"""

#| - Import Modules
import numpy as np

from ase_modules.ase_methods import create_species_element_dict
#__|

def calc_formation_energy(
    atoms,
    energy,
    reference_states,
    normalize_per_atom=True,
    ):
    """Calculate formation energy of structure on a per atom basis.

        reference states = [
            {
                "energy": ,
                "element_dict": {"H": 1, "Ir": 2},
                },
            {

                },
            ...
            ]

    Args:
        atoms:
        energy:
        reference_dict:
        normalize_per_atom:
    """
    #| - calc_formation_energy
    atoms.info["element_dict"] = create_species_element_dict(atoms)

    num_atoms = atoms.get_number_of_atoms()

    # ordered_elems = list(atoms.info["element_dict"])
    # ordered_elems.sort()


    #| - Finding List of Unique Elements Defined by reference_states list
    ref_state_elem_list = []
    for i in reference_states:
        ref_i_elems = list(i["element_dict"])
        ref_state_elem_list.extend(ref_i_elems)

    ordered_elems = list(set(ref_state_elem_list))
    ordered_elems.sort()
    #__|


    #| - Constructing the A matrix
    a_matrix = []
    for ref_state_i in reference_states:
        ref_i_comp_vect = []
        for elem_i in ordered_elems:
            if elem_i in list(ref_state_i["element_dict"]):
                elem_i_num = ref_state_i["element_dict"][elem_i]
                ref_i_comp_vect.append(elem_i_num)
            else:
                ref_i_comp_vect.append(0.)

        a_matrix.append(np.array(ref_i_comp_vect))
    a_matrix = np.array(a_matrix).transpose()
    #__|

    #| - Constructing the b vector
    # phi = 1.
    b_vect = []
    for elem_i in ordered_elems:
        if elem_i in list(atoms.info["element_dict"]):
            elem_i_num = atoms.info["element_dict"][elem_i]
            b_vect.append(elem_i_num)
        else:
            b_vect.append(0.)

    b_vect = np.array(b_vect)
    #__|

    #| - Solve linear system of equations
    x = np.linalg.solve(
        a_matrix,
        b_vect,
        )
    #__|

    #| - Calculate Formation Energy
    ref_e_sum = 0.
    for coeff_i, ref_i in zip(x, reference_states):
        ref_i_contribution = coeff_i * ref_i["elec_e"]
        ref_e_sum += ref_i_contribution

    form_e = energy - ref_e_sum

    if normalize_per_atom:
        form_e = form_e / num_atoms
    #__|


    return(form_e)
    #__|
