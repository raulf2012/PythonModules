#!/usr/bin/env python

"""Methods for the regression of Lennard Jones parameters to DFT energies.

Author: Raul A. Flores
"""

#| - Import Modules
import gc
import pickle

import itertools

import numpy as np
import pandas as pd
from scipy.optimize import minimize

#My Modules
from energetics.formation_energy import calc_formation_energy
from classical_methods.lennard_jones import lennard_jones_sp

from IPython.display import display, clear_output

#| - __old__
# import os
# import copy
# from scipy import stats
# from scipy.optimize import fmin
# import collections
# import plotly.plotly as py
# import plotly.graph_objs as go
# from ase.visualize import view
# from ase import io
# from asap3 import LennardJones
# import math
# from an_data_processing import load_df
# from ase_modules.ase_methods import create_species_element_dict
#__|

#__|

def calc_lennard_jones_form_e(
    atoms_i=None,

    atoms_H2=None,
    atoms_Ir=None,
    atoms_O2=None,

    epsilon=None,
    sigma=None,

    ):
    """Calculate the Lennard Jones formation energy of atoms object.

    Args:
        atoms_i:
        atoms_H2:
        atoms_Ir:
        atoms_O2:
        epsilon:
        sigma:
    """
    #| - calc_lennard_jones_form_e
# atoms_H2 = io.read("reference_states/H2_mp-632172_conventional_standard.cif")
# atoms_Ir = io.read("reference_states/Ir_mp-101_conventional_standard.cif")
# atoms_O2 = io.read("reference_states/O2_mp-12957_conventional_standard.cif")

    E_H = lennard_jones_sp(epsilon, sigma, atoms_H2)
    E_O = lennard_jones_sp(epsilon, sigma, atoms_O2)
    E_Ir = lennard_jones_sp(epsilon, sigma, atoms_Ir)

    reference_states = [
        {
            "elec_e": E_H,
            "element_dict": {"H": 1},
            },

        {
            "elec_e": E_O,
            "element_dict": {"O": 1},
            },

        {
            "elec_e": E_Ir,
            "element_dict": {"Ir": 1},
            },

        ]

    # E_IrOxHy = lennard_jones_sp(
    #     epsilon,
    #     sigma,
    #     atoms_i,
    #     normalize_per_atom=False,
    #     )

    E_form_i = calc_formation_energy(
        atoms_i,
        lennard_jones_sp(epsilon, sigma, atoms_i, normalize_per_atom=False),
        reference_states,
        normalize_per_atom=False,
        )

    E_form_per_atom_i = E_form_i / atoms_i.get_number_of_atoms()

    E_out = E_form_per_atom_i

    return(E_out)
    #__|

def calc_lennard_jones_all_atoms(
    pars,
    atoms_list,
    reference_atoms,
    ):
    """Calculate the Lennard Jones formation energies of list of atoms.

    Args:
        pars:
        atoms_list:
        reference_atoms:
    """
    #| - calc_lennard_jones_all_atoms
    epsilon = pars[0]
    sigma = pars[1]

    atoms_H2 = reference_atoms[0]
    atoms_O2 = reference_atoms[1]
    atoms_Ir = reference_atoms[2]

    predicted_energies = []
    for atoms_i in atoms_list:

        lj_energy_i = calc_lennard_jones_form_e(
            atoms_i=atoms_i,

            atoms_H2=atoms_H2,
            atoms_Ir=atoms_Ir,
            atoms_O2=atoms_O2,

            epsilon=epsilon,
            sigma=sigma,
            )

#         lj_energy_i = lennard_jones_sp(epsilon, sigma, atoms_i)

        predicted_energies.append(lj_energy_i)

    predicted_energies = np.array(predicted_energies)

    return(predicted_energies)
    #__|

def objective(
    pars,
    known_energies,
    atoms_list,
    eps_shape,
    sig_shape,
    elem_list,
    reference_atoms,
    info,
    ):
    """Objective function to be minimized for LJ parameter fitting.

    Args:
        pars:
        known_energies:
        atoms_list:
        eps_shape:
        sig_shape:
        reference_atoms:
    """
    #| - objective
    epsilon, sigma = unflatten_eps_sig_array(
        pars,
        eps_shape,
        sig_shape,
        elem_list,
        )

    err = known_energies - \
        calc_lennard_jones_all_atoms(
            (epsilon, sigma),
            atoms_list,
            reference_atoms,
            )

    MSE = np.mean(err ** 2)

    # gc.collect()

    clear_output(wait=True)
    display("Iter: " + str(info["Nfeval"]))
    display("MSE: " + str(MSE))
    display("Epsilon Matrix: ")
    display(epsilon)
    display("Sigma Matrix: ")
    display(sigma)
    display("__________________________")

    display(epsilon.values)
    display(sigma.values)

    info["Nfeval"] += 1

    return(MSE)
    #__|

def flatten_eps_sig_triangular_matrices(
    epsilon,
    sigma,
    mode="traingular",  # 'triangular' or 'diagonal'
    ):
    """Flatten triangular epsilon and sigma matrices into a 1D array.

    Args:
        epsilon:
        sigma:
    """
    #| - flatten_eps_sig_triangular_matrices
    epsilon = epsilon.values
    sigma = sigma.values

    if mode == "triangular":
        flat_eps_sig = np.hstack([
            epsilon.flatten(),
            sigma.flatten(),
            ])

        # Remove 0s
        flat_eps_sig_no_0s = flat_eps_sig[flat_eps_sig != 0.]

        return(flat_eps_sig_no_0s)

    elif mode == "diagonal":
        flat_eps_sig_diag = np.hstack(
            [
                epsilon.diagonal(),
                sigma.diagonal(),
                ]
            )

        return(flat_eps_sig_diag)
    #__|

def unflatten_eps_sig_array(
    flat_eps_sig,
    eps_shape,
    sig_shape,
    elem_list,
    ):
    """Unflatten a 1D array into 2 triangular epsilon and sigma matrices.

    Args:
        flat_eps_sig:
        eps_shape:
        sig_shape:
    """
    #| - unflatten_eps_sig_array

    #| - Array Dimension Check
    assert eps_shape[0] == eps_shape[1]
    assert sig_shape[0] == sig_shape[1]

    N_eps = eps_shape[0]
    N_sig = sig_shape[0]

    assert N_eps == N_sig

    N = N_eps
    #__|

    len_pars = len(flat_eps_sig)
    half_way = int(len_pars / 2)

    epsilon_short = flat_eps_sig[:half_way]
    sigma_short = flat_eps_sig[half_way:]

    if len(epsilon_short) == N and len(sigma_short) == N:
        pars_mode = "diagonal"
    else:
        pars_mode = "triangular"

    #| - Methods
    def unflatten_tri_matrix_with_defined_cross_terms(
        flat_array,
        N,
        cross_terms_mode="geo"  # "geo" or "ave"
        ):
        """Convert array into a diagonal matrix with defined lower cross-terms.

        The lower-half cross terms are calculated from the diagonal terms.

        Args:
            flat_array:
            N:
            cross_terms_mode:
                "geo" or "ave"
        """
        #| - unflatten_tri_matrix_with_defined_cross_terms
        matrix = np.diag(flat_array)

        # Return list of i, j indices corresponding to the off diagonal
        # cross-terms in the lower section of a triangular matrix
        all_comb_indices = list(itertools.product(
            np.array(range(N)),
            np.array(range(N)),
            ))

        unique_ij_list = []
        for i_j_pair in all_comb_indices:
            i_ind = i_j_pair[0]
            j_ind = i_j_pair[1]
            if i_ind == j_ind:
                continue
            unique_ij_list.append(set(i_j_pair))

        unique_ij_list = list(set([frozenset(item) for item in unique_ij_list]))
        unique_ij_list = [list(item) for item in unique_ij_list]

        for i in unique_ij_list:
            i.sort(reverse=True)

        for ij_ind in unique_ij_list:
            i_ind = ij_ind[0]
            j_ind = ij_ind[1]

            i_term = matrix[i_ind][i_ind]
            j_term = matrix[j_ind][j_ind]

            average_ij = (i_term + j_term) / 2.
            geo_ave_ij = (i_term * j_term) ** (1. / 2.)

            if cross_terms_mode == "geo":
                matrix[i_ind][j_ind] = geo_ave_ij
            elif cross_terms_mode == "ave":
                matrix[i_ind][j_ind] = average_ij

        return(matrix)
        #__|

    def unflatten_single_triangular_matrix(flat_array, N):
        """Unflatten a single triangular matrix.

        Args:
            flat_array:
            N:
        """
        #| - unflatten_single_traingular_matrix
        start_index_list = []
        stop_index_list = []
        j_cnt = 0
        for i in range(N):
            start_index_list.append(j_cnt)
            increment = (i + 1)
            j_cnt += increment
            stop_index_list.append(j_cnt)

        rebuilt_matrix = []
        for i_ind, (start_i, stop_i) in enumerate(zip(
            start_index_list,
            stop_index_list)):

            num_of_0_to_add = N - i_ind - 1

            final_row = np.append(
                flat_array[start_i:stop_i],
                num_of_0_to_add * [0.],
                )

            rebuilt_matrix.append(final_row)

        rebuilt_matrix = np.array(rebuilt_matrix)

        return(rebuilt_matrix)
        #__|

    #__|

    if pars_mode == "triangular":
        epsilon = unflatten_single_triangular_matrix(epsilon_short, N)
        sigma = unflatten_single_triangular_matrix(sigma_short, N)

    elif pars_mode == "diagonal":

        epsilon = unflatten_tri_matrix_with_defined_cross_terms(
            epsilon_short,
            N,
            cross_terms_mode="ave",
            )

        sigma = unflatten_tri_matrix_with_defined_cross_terms(
            sigma_short,
            N,
            cross_terms_mode="geo",
            )

    epsilon = pd.DataFrame(
        epsilon,
        index=elem_list,
        columns=elem_list,
        )

    sigma = pd.DataFrame(
        sigma,
        index=elem_list,
        columns=elem_list,
        )

    return(epsilon, sigma)

    #__|

def fit_LJ_to_DFT(
    objective=None,
    known_energies=None,
    atoms_list=None,
    elem_list=None,
    reference_atoms=None,
    epsilon_0=None,
    sigma_0=None,
    tol=1e-4,
    maxiter=50,
    maxfun=20,
    params_mode="triangular",  # "triangular" or "diagonal"
    ):
    """Fit LJ parameters to DFT energies.

    Args:
        objective:
        known_energies:
        atoms_list:
        reference_atoms:
        epsilon_0:
        sigma_0:
        tol:
        maxiter:
        maxfun:
    """
    #| - fit_LJ_to_DFT
    atoms_H2 = reference_atoms[0]
    atoms_O2 = reference_atoms[1]
    atoms_Ir = reference_atoms[2]

    epsilon = epsilon_0
    sigma = sigma_0

    eps_shape = epsilon.shape
    sig_shape = sigma.shape


    # if params_mode == "triangular":
    pars = flatten_eps_sig_triangular_matrices(
        epsilon,
        sigma,
        mode=params_mode,  # 'triangular' or 'diagonal'
        )

    # elif params_mode == "diagonal":
    # pars = flatten_eps_sig_triangular_matrices(epsilon, sigma)


    #| - Minimize Method
    opt_out = minimize(
        objective,
        pars,

        args=(
            known_energies,
            atoms_list,
            eps_shape,
            sig_shape,
            elem_list,
            [
                atoms_H2,
                atoms_O2,
                atoms_Ir,
                ],

            {'Nfeval': 0},
            ),

        # method=None,
        method='L-BFGS-B',

        jac=None,
        hess=None,
        hessp=None,

        # I'm adding a bit to the lower bound to avoid 0. values
        # bounds=[(0, None) for i in pars],
        bounds=[(0.0001, None) for i in pars],

        constraints=(),
        tol=tol,
        callback=None,
        options={
            "maxiter": maxiter,
            "maxfun": maxfun,

            "disp": True,
            },
        )
    #__|

    LJ_pars = opt_out.x

    epsilon_out, sigma_out = unflatten_eps_sig_array(
        LJ_pars,
        eps_shape,
        sig_shape,
        elem_list,
        )

    with open("opt_params.pickle", "wb") as fle:
        pickle.dump((epsilon_out, sigma_out), fle)

    #| - __old__
    # new_energies = calc_lennard_jones_all_atoms(
    #     (epsilon_out, sigma_out),
    #     atoms_list,
    #     [atoms_H2, atoms_O2, atoms_Ir],
    #     )
    #__|

    return(epsilon_out, sigma_out)
    #__|

def calc_MSE(
    pars,
    df_i,
    DFT_energies_col,
    ref_atoms_list,
    ):
    """Calculate the mean-squared-error of model on dataset.

    Args:
        pars:
        df_i:
        DFT_energies_col:
        ref_atoms_list:
    """
    #| - calc_MSE
    # df_i = df_orig.iloc[cv_folds_indices[0]["testing"]]
    epsilon = pars[0]
    sigma = pars[1]

    known_energies = np.array(df_i[DFT_energies_col].tolist())
    atoms_list = df_i["final_atoms"].tolist()

    new_energies_test = calc_lennard_jones_all_atoms(
        (epsilon, sigma),
        atoms_list,
        ref_atoms_list,
        # [atoms_H2, atoms_O2, atoms_Ir],
        )

    err = known_energies - new_energies_test
    MSE = np.mean(err ** 2)

    return(MSE)
    #__|

def k_fold_cross_validation(data, k=5):
    """k-fold cross-validation indices list.

    Args:
        data:
        k:
    """
    #| - k_fold_cross_validation
    folds = np.array_split(data, k)

    cv_data = []
    for i_cnt in range(k):

        training_data_i = []
        for j_cnt in np.delete(np.arange(k), i_cnt):
            training_data_i = np.concatenate((training_data_i, folds[j_cnt]))
        testing_data_i = folds[i_cnt]

        cv_data.append({
            "training": training_data_i,
            "testing": testing_data_i,
            })

    return(cv_data)
    #__|
