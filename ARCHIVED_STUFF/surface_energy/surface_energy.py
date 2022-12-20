#!/usr/bin/env python

"""Module to calclate surface energies from symmetric slabs.
TEMP
Author: Raul A. Flores
"""

# | - IMPORT MODULES
import numpy as np
import pandas as pd

import plotly.graph_objs as go

from ase_modules.ase_methods import create_species_element_dict

from pymatgen.core.composition import Composition
# __|


class SurfaceEnergy:
    """Calculates surface energy of slab relative to a given bulk structure.

    TODO:
      Currently assumes that the slab is symmetric (each surface is the same),
      try to add check for this
    """

    # | - SurfaceEnergy ********************************************************
    _TEMP = "TEMP"

    def __init__(self,
        # | - TEMP | ARGS ------------------------------------------------------
        atoms=None,
        electronic_energy=None,  # Only needed if not within atoms object

        bulk_atoms=None,
        bulk_electronic_energy_per_atom=None,

        H_ref_electronic_energy=None,
        O_ref_electronic_energy=None,

        special_surface_species_corrections=None,
        apply_special_species_corrections=False,
        num_surface_atoms=None,

        bias=0.,
        pH=0.,

        verbose=True,
        # __| ------------------------------------------------------------------
        ):
        """
        """
        # | - __init__

        # | - Setting Argument Instance Attributes
        self.atoms = atoms
        self.electronic_energy = electronic_energy

        self.bulk_atoms = bulk_atoms
        self.bulk_electronic_energy_per_atom = bulk_electronic_energy_per_atom

        self.H_ref_electronic_energy = H_ref_electronic_energy
        self.O_ref_electronic_energy = O_ref_electronic_energy
        self.special_surface_species_corrections = \
            special_surface_species_corrections
        self.apply_special_species_corrections = \
            apply_special_species_corrections

        self.num_surface_atoms = num_surface_atoms

        self.bias = bias
        self.pH = pH

        self.verbose = verbose
        # __|

        # | - Initializing Internal Instance Attributes
        self.__num_atoms_reduced_bulk__ = None

        self.non_stoich_comp = None

        self.slab_thickness = None

        self.__bulk_energy_per_atom__ = None

        self.non_stoich_comp_new = dict()
        # __|

        self.surface_area = self.__calc_surface_area__()

        self.electronic_energy = self.__get_electronic_energy__(
            self.atoms, self.electronic_energy)

        self.bulk_electronic_energy = self.__get_electronic_energy__(
            self.bulk_atoms, self.bulk_electronic_energy_per_atom,
            per_atom=True)

        self.__bulk_formula_units_in_slab__ = \
            self.__calc_units_of_bulk_in_slab__()

        self.__bulk_formula_reduction__ = \
            self.__calc_units_of_reduced_bulk_in_bulk__()

        self.__bulk_energy_per_formula_unit__ = \
            self.__calc_bulk_energy_per_formula_unit__()

        self.special_surface_species = \
            self.__count_special_surface_species__()

        self.std_surface_e_per_side = self.calc_std_surface_energy()

        self.std_surface_e_per_area = \
            self.__calc_surface_energy_per_area__(
                unnorm_surface_e=self.std_surface_e_per_side)


        if self.num_surface_atoms is not None:
            self.surface_e_per_surface_atom = \
                self.__calc_std_surface_energy_per_surface_atom__()

        self.slab_thickness = self.__calc_slab_thickness__()
        # __|

    def calc_surface_energy(self, bias, pH, norm_type="area"):
        """

        Args:
            norm_type: 'area', 'surface_atom', None
        """
        # | - calc_surface_energy
        surface_e_per_side = self.std_surface_e_per_side
        non_stoich_comp = self.non_stoich_comp

        # print(surface_e_per_side)

        nonstoich_O = non_stoich_comp.get("O", 0)
        nonstoich_H = non_stoich_comp.get("H", 0)

        slope = 2 * nonstoich_O - nonstoich_H

        surf_e_V_ph = 0. + \
            +surface_e_per_side + \
            +(slope * 0.0591 * pH) + \
            -(slope * bias) + \
            +0.

        # print(surf_e_V_ph)


        if norm_type == "area":
            surf_e_V_ph__norm = self.__calc_surface_energy_per_area__(
                unnorm_surface_e=surf_e_V_ph)
        elif norm_type == "surface_atom":
            print("NOT SUPPORTEE")
            assert False
            surf_e_V_ph__norm = "TEMP"
        elif norm_type is None:
            surf_e_V_ph__norm = surf_e_V_ph
        else:
            assert False
            surf_e_V_ph__norm = "TEMP"


        return(surf_e_V_ph__norm)
        # __|

    def calc_std_surface_energy(self):
        """
        """
        # | - calc_std_surface_energy

        electronic_energy = self.electronic_energy
        bulk_formula_units_in_slab = self.__bulk_formula_units_in_slab__
        bulk_energy_per_formula_unit = self.__bulk_energy_per_formula_unit__
        H_ref_electronic_energy = self.H_ref_electronic_energy
        O_ref_electronic_energy = self.O_ref_electronic_energy


        non_stoich_comp = self.non_stoich_comp
        # non_stoich_comp_new = self.non_stoich_comp_new
        special_surface_species = self.special_surface_species


        # TODO add this as class input
        apply_special_species_corrections = \
            self.apply_special_species_corrections

        special_surface_species_corrections = \
            self.special_surface_species_corrections
        # special_surface_species_corrections = {
        #     "*OH": 0.2945, "*O": 0.044, "*OOH": 0.3765}
        sssc = special_surface_species_corrections


        # # TEMP
        # print(non_stoich_comp)
        # print(H_ref_electronic_energy)

        # TODO Make the referencing more robust, take arbitary dict of
        # referencec atoms
        surf_e_0 = 0. + \
            +electronic_energy + \
            -bulk_formula_units_in_slab * bulk_energy_per_formula_unit + \
            -non_stoich_comp.get("O", 0.) * (O_ref_electronic_energy) + \
            -non_stoich_comp.get("H", 0.) * (H_ref_electronic_energy) + \
            +0.

        print("TEMP | ksjdfidsihjgisdfjgf8sgsd | TEMP")
        print("electronic_energy:", electronic_energy)
        print("bulk_formula_units_in_slab:", bulk_formula_units_in_slab)
        print("bulk_energy_per_formula_unit", bulk_energy_per_formula_unit)
        print("std_surface_e_per_side:", surf_e_0)

        if apply_special_species_corrections:
            for spec_i, num_spec_i in special_surface_species.items():
                corr_i = num_spec_i * sssc[spec_i]
                print("Special spec corr:", spec_i, corr_i)
                # surf_e_0 += corr_i
                surf_e_0 -= corr_i

        print("std_surface_e_per_side__after_corr:", surf_e_0)
        print("TEMP | ksjdfidsihjgisdfjgf8sgsd | TEMP")


        # Divide surface energy across two sides of slab
        surf_e_0 /= 2

        return(surf_e_0)
        # __|


    def __count_special_surface_species__(self):
        """Count the number of *OH and *O surface species.

        Taks the non-stoich dict and pairs the max number of O and H pairs into
        *OH and the remaining O's into *O
        """
        # | - __count_special_surface_species__
        import copy

        non_stoich_comp = self.non_stoich_comp

        non_stoich_comp_new = copy.copy(non_stoich_comp)

        # print(non_stoich_comp)

        special_species_dict = dict()
        if "O" in non_stoich_comp.keys():
            #| - If *O present
            num_Os = non_stoich_comp.get("O")

            if "H" in non_stoich_comp.keys():
                num_Hs = non_stoich_comp.get("H")

                min_num = min([num_Os, num_Hs])

                num_OHs = min_num

                left_over_Hs = num_Hs - min_num
                left_over_Os = num_Os - min_num

                special_species_dict["*OH"] = num_OHs
                special_species_dict["*O"] = left_over_Os
                special_species_dict["*H"] = left_over_Hs


                # All nonstoich Os will be *O species
                non_stoich_comp_new["O"] = 0
                non_stoich_comp_new["H"] = left_over_Hs
            else:
                num_OHs = 0
                special_species_dict["*OH"] = num_OHs

                left_over_Hs = 0
                left_over_Os = num_Os

                special_species_dict["*O"] = left_over_Os
                special_species_dict["*OH"] = 0

                # All nonstoich Os will be *O species
                non_stoich_comp_new["O"] = 0
                non_stoich_comp_new["H"] = left_over_Hs
            #__|



        #| - __old__
        # else:
        #     num_OHs = 0
        #     left_over_Os = num_Os
        #     left_over_Hs = 0
        #
        #     if "H" in non_stoich_comp.keys():
        #         if non_stoich_comp.get("H") > 0:
        #             raise ValueError(
        #                 "NOT GOOD HERE, THERE IS AN *H WITHOUT and *OH")


        # print("----")
        # print(non_stoich_comp_new)
        # print(special_species_dict)
        #__|

        self.non_stoich_comp_new = non_stoich_comp_new

        return(special_species_dict)
        # __|

    def __calc_surface_energy_per_area__(self, unnorm_surface_e=None):
        """Normalize the surface energy to surface area (A^2)."""
        # | - __calc_surface_energy_per_area__
        surface_area = self.surface_area
        # surface_e_per_side = self.std_surface_e_per_side
        surface_e_per_side = unnorm_surface_e

        surface_e_per_area = surface_e_per_side / surface_area

        return(surface_e_per_area)
        # __|


    def __calc_std_surface_energy_per_surface_atom__(self):
        """Normalize the surface area to a per surface atom basis."""
        # | - calc_std_surface_energy_per_area
        num_surface_atoms = self.num_surface_atoms
        surface_e_per_side = self.std_surface_e_per_side

        surface_e_per_surface_atoms = surface_e_per_side / num_surface_atoms

        return(surface_e_per_surface_atoms)
        # __|


    def __calc_surface_area__(self):
        """
        """
        # | - __calc_surface_area__
        atoms = self.atoms

        cell = atoms.cell

        cross_prod_i = np.cross(cell[0], cell[1])
        area_i = np.linalg.norm(cross_prod_i)

        return(area_i)
        # __|


    def __get_electronic_energy__(self,
        atoms,
        electronic_energy,
        per_atom=False):
        """
        """
        # | - __get_electronic_energy__
        energy_out = None
        if electronic_energy is not None:

            if self.verbose:
                print("Using energy provided instead of energy attached to atoms")

            if per_atom:
                energy_out = atoms.get_number_of_atoms() * electronic_energy
            else:
                energy_out = electronic_energy
        else:
            try:
                energy_out = atoms.get_potential_energy()
            except RuntimeError as err:
                print(err, "AHHHHHHHH!!!!")
                raise ValueError('No where to get energy from!!!')
            except:
                print("Uncaught error type")
                raise ValueError('No where to get energy from!!!')

        return(energy_out)
        # __|

    def __calc_units_of_bulk_in_slab__(self):
        """
        """
        # | - __calc_units_of_bulk_in_slab__
        # TODO
        main_atom = "Ir"  # Make this a class attribute

        # 'main_atom' or 'gcm' (greatest common multiple)
        find_bulk_form_units_method = "main_atom"

        bulk_atoms = self.bulk_atoms
        atoms = self.atoms


        comp0 = Composition(bulk_atoms.get_chemical_formula())

        df = pd.DataFrame([
            create_species_element_dict(atoms,
                elems_to_always_include=["O", "H"]),
            dict(comp0.to_data_dict["reduced_cell_composition"])],
            index=["slab", "bulk"])

        # Replace NaNs with 0.
        df = df.replace(np.nan, 0.0, regex=True)

        # Removingg columns with 0
        df = df.loc[:, (df != 0).any(axis=0)]

        # slab_comp_array = np.array(list(df.loc["slab"]))
        # bulk_comp_array = np.array(list(df.loc["bulk"]))
        # df.loc["slab"].to_numpy()
        # df.loc["bulk"].to_numpy()

        # Number of unit of the bulk's reduced formula that fit into the slab
        if find_bulk_form_units_method == "main_atom":
            bulk_formula_units_in_slab = int(
                df.loc["slab"][main_atom] / df.loc["bulk"][main_atom]
                )

        elif find_bulk_form_units_method == "gcm":
            bulk_formula_units_in_slab = int(min(
                df.loc["slab"].to_numpy() / df.loc["bulk"].to_numpy()
                ))
        bfuis = bulk_formula_units_in_slab

        # #####################################################################
        # Getting the non-stoicheometric atoms composition
        df.loc["nonstoich"] = df.loc["slab"] - bfuis * df.loc["bulk"]
        non_stoich_comp = df.loc["nonstoich"].to_dict()
        self.non_stoich_comp = non_stoich_comp

        # print(bulk_formula_units_in_slab)
        # print(non_stoich_comp)

        return(bulk_formula_units_in_slab)
        # __|

    def __calc_units_of_reduced_bulk_in_bulk__(self):
        """
        """
        # | - __calc_units_of_reduced_bulk_in_bulk__
        bulk_atoms = self.bulk_atoms

        bulk_atoms.get_chemical_formula()

        comp0 = Composition(bulk_atoms.get_chemical_formula())


        reduced_comp = dict(comp0.to_data_dict["reduced_cell_composition"])
        num_atoms_in_reduc_form = int(sum(list(reduced_comp.values())))

        self.__num_atoms_reduced_bulk__ = num_atoms_in_reduc_form


        df = pd.DataFrame([
            dict(comp0.to_data_dict["reduced_cell_composition"]),
            dict(comp0.to_data_dict["unit_cell_composition"])],
            index=["reduced", "original"])

        # Replace NaNs with 0.
        df = df.replace(np.nan, 0.0, regex=True)

        # Removingg columns with 0
        df = df.loc[:, (df != 0).any(axis=0)]

        reduced_comp_array = np.array(list(df.loc["reduced"]))
        orig_comp_array = np.array(list(df.loc["original"]))

        # Number of unit of the bulk's reduced formula that fit into the slab
        bulk_reduction_factor = int(min(orig_comp_array / reduced_comp_array))
        # print(bulk_reduction_factor)

        return(bulk_reduction_factor)
        # __|

    def __calc_bulk_energy_per_formula_unit__(self):
        """
        """
        # | - __calc_bulk_energy_per_formula_unit__
        # bulk_formula_reduction = self.__bulk_formula_reduction__
        bulk_electronic_energy = self.bulk_electronic_energy
        num_atoms_reduced_bulk = self.__num_atoms_reduced_bulk__
        bulk_atoms = self.bulk_atoms

        num_atoms_bulk = bulk_atoms.get_number_of_atoms()

        # norm_fact = (E_per_atom) * (#_atoms_in_reduced_formula)
        norm_fact = (num_atoms_reduced_bulk / num_atoms_bulk)
        bulk_electronic_energy_per_formula = bulk_electronic_energy * norm_fact

        self.__bulk_energy_per_atom__ = \
            bulk_electronic_energy / num_atoms_bulk

        return(bulk_electronic_energy_per_formula)
        # __|

    def __calc_slab_thickness__(self):
        """Calculate the thickness of the atoms slab."""
        # | - __calc_slab_thickness__
        atoms = self.atoms

        positions = atoms.positions

        slab_thickness = max(positions[:, 2]) - min(positions[:, 2])

        return(slab_thickness)
        # __|



    # __| **********************************************************************



class SurfaceEnergyConvergence:
    """Compute properties from series of slabs of increasing thickness.

    Regress bulk energy from total E's vs N
        Extracting convergent surface energies from slab calculations
        https://iopscience.iop.org/article/10.1088/0953-8984/8/36/005
    """

    # | - SurfaceEnergyConvergence *********************************************
    _TEMP = "TEMP"

    def __init__(self,
        SurfaceEnergy_instances=None,
        bulk_electronic_energy_per_atom=None,
        num_points_to_exclude=1,
        verbose=True,
        ):
        """

        Args:
            SurfaceEnergy_instances:
            num_points_to_exclude: int
                Number of points to exclude for the bulk energy regression
                Thin slabs will suffer from finite size effects, so they should
                not be included in regression
        """
        # | - __init__

        # | - Setting Argument Instance Attributes
        self.SurfaceEnergy_instances = SurfaceEnergy_instances
        self.num_points_to_exclude = num_points_to_exclude
        self.bulk_electronic_energy_per_atom = bulk_electronic_energy_per_atom
        self.verbose = verbose
        # __|

        # | - Initializing Internal Instance Attributes
        self.sufficient_data_to_fit_bulk = True
        self.fitted_bulk_energy = None
        self.new_SurfaceEnergy_instances = []
        self.new_ave_surface_energy_per_area = None
        # __|

        self.df = self.__init_dataframe__()

        self.ave_surface_energy_per_area = \
            self.__calc_ave_surface_energy__(
                self.SurfaceEnergy_instances)

        # Enough data to fit?
        self.__sufficient_data_to_fit_bulk()

        # __|


    def calculate_surface_energies(self, bulk_energy=None):
        """
        """
        # | - calculate_surface_energies
        self.new_SurfaceEnergy_instances = \
            self.__recalc_SurfaceEnergy_w_new_bulk__(
                bulk_energy)

        # self.fitted_bulk_energy)

        self.new_ave_surface_energy_per_area = \
            self.__calc_ave_surface_energy__(
                self.new_SurfaceEnergy_instances)
        # __|

    def fit_bulk_energy(self):
        """
        """
        # | - fit_bulk_energy
        if self.sufficient_data_to_fit_bulk:
            self.fitted_bulk_energy = self.__calc_regressed_bulk_energy__()
        # __|

    def inst_surf_e_with_fitted_bulk(self):
        """
        COMBAK: This isn't being used now
        """
        # | - tmp_meth
        # fitted_bulk_energy = self.fitted_bulk_energy
        # SurfaceEnergy_instances = self.SurfaceEnergy_instances
        #
        # SE_old_i = SurfaceEnergy_instances[0]
        #
        # atoms = SE_old_i.atoms
        #
        #
        # SurfaceEnergy(
        #     atoms=atoms,
        #     # electronic_energy=None,
        #
        #     bulk_atoms=None,
        #     bulk_electronic_energy=None,
        #
        #     H_ref_electronic_energy=None,
        #     O_ref_electronic_energy=None,
        #
        #     # num_surface_atoms=None,
        #     # bias=0.,
        #     # pH=0.,
        #     )

        # __|


    def __init_dataframe__(self):
        """
        """
        # | - __init_dataframe__
        SurfaceEnergy_instances = self.SurfaceEnergy_instances

        atoms_list = [i.atoms for i in SurfaceEnergy_instances]
        number_of_atoms = [i.get_number_of_atoms() for i in atoms_list]
        potential_energy_list = [i.get_potential_energy() for i in atoms_list]

        df = pd.DataFrame()

        df["number_of_atoms"] = number_of_atoms
        df["potential_energy"] = potential_energy_list

        df = df.sort_values("number_of_atoms")

        return(df)
        # __|

    def __sufficient_data_to_fit_bulk(self):
        """
        """
        # | - __sufficient_data_to_fit_bulk
        df = self.df
        num_points_to_exclude = self.num_points_to_exclude

        # If there are too few points excluding points may not be the best
        # Need at least 3 points after removing points to do fit
        if len(df) - num_points_to_exclude < 3:

            # | - Less than 3 data points, don't fit
            if len(df) <= 2:
                self.num_points_to_exclude = 0
                self.sufficient_data_to_fit_bulk = False

                if self.verbose:
                    print("Only ", str(len(df)), " data points, in dataframe")
                    print("Will not fit bulk energy")
            # __|

            # | - Modifiying num_points_to_exclude to have at least 3 points
            else:
                num_points_to_exclude = int(len(df) - 3)
                self.num_points_to_exclude = num_points_to_exclude

                if self.verbose:
                    print(
                        "Changed num_points_to_exclude to ",
                        num_points_to_exclude)
                    # print("")
            # __|

        # __|

    def __calc_regressed_bulk_energy__(self):
        """
        """
        # | - __calc_regressed_bulk_energy__
        df = self.df
        num_points_to_exclude = self.num_points_to_exclude
        nbpe = num_points_to_exclude

        number_of_atoms = df["number_of_atoms"]
        potential_energy = df["potential_energy"]

        x_i = number_of_atoms
        y_i = potential_energy

        z = np.polyfit(x_i[nbpe:], y_i[nbpe:], 1)

        bulk_energy = z[0]

        return(bulk_energy)
        # __|


    def __recalc_SurfaceEnergy_w_new_bulk__(self, new_bulk_energy):
        """
        """
        # | - __recalc_SurfaceEnergy_w_new_bulk__
        SurfaceEnergy_instances = self.SurfaceEnergy_instances
        # new_bulk_energy = self.new_bulk_energy
        verbose = self.verbose

        new_SurfaceEnergy_instances = []
        for SurfaceEnergy_instance_i in SurfaceEnergy_instances:
            SE_old_i = SurfaceEnergy_instance_i

            atoms = SE_old_i.atoms
            bulk_atoms = SE_old_i.bulk_atoms
            H_ref_electronic_energy = SE_old_i.H_ref_electronic_energy
            O_ref_electronic_energy = SE_old_i.O_ref_electronic_energy

            SE_new = SurfaceEnergy(
                atoms=atoms,
                bulk_atoms=bulk_atoms,
                bulk_electronic_energy_per_atom=new_bulk_energy,
                H_ref_electronic_energy=H_ref_electronic_energy,
                O_ref_electronic_energy=O_ref_electronic_energy,
                verbose=verbose,
                )

            new_SurfaceEnergy_instances.append(SE_new)

        return(new_SurfaceEnergy_instances)
        # __|


    def __calc_ave_surface_energy__(self, SurfaceEnergy_instances):
        """
        """
        # | - __calc_ave_surface_energy__
        # new_SurfaceEnergy_instances = self.new_SurfaceEnergy_instances
        # SurfaceEnergy_instances = self.SurfaceEnergy_instances
        # if

        y_surface_e = []
        x_slab_thickness = []
        for SE_inst_i in SurfaceEnergy_instances:
            y_surface_e.append(SE_inst_i.std_surface_e_per_area)
            x_slab_thickness.append(SE_inst_i.slab_thickness)

        df = pd.DataFrame()

        df["surface_energy_per_area"] = y_surface_e
        df["slab_thickness"] = x_slab_thickness


        # Number of points in the Energy vs  Slab_thickness plot to include
        # in the average for the new surface energy
        num_points_to_average = 3
        nptoa = num_points_to_average

        df = df.sort_values("slab_thickness")
        ave_surface_e = df.iloc[-nptoa:]["surface_energy_per_area"].mean()

        return(ave_surface_e)
        # __|

    def plot_E_vs_N_convergence(self):
        """Plot raw slab energy vs number of atoms.

        The slope of this plot is the fitted bulk energy
        """
        # | - plot_E_vs_N_convergence
        df = self.df

        number_of_atoms = df["number_of_atoms"]
        potential_energy_list = df["potential_energy"]

        trace_i = go.Scatter(
            x=number_of_atoms,
            y=potential_energy_list,
            mode='markers+lines',
            # name=name_i,

            marker=dict(
                symbol="square",
                size=10,
                # color=color,
                line=dict(
                    width=1,
                    color='rgb(0, 0, 0)',
                    ),
                ),

            )

        return(trace_i)
        # __|


    def plot_surface_energy(self,
        name_i="TEMP",
        color_i="pink"
        ):
        """
        """
        # | - plot_surface_energy
        data = []

        # | - Surface Energy (DFT Bulk) ########################################
        y_surface_e = []; x_slab_thickness = []
        for SE_inst_i in self.SurfaceEnergy_instances:
            y_surface_e.append(SE_inst_i.std_surface_e_per_area)
            x_slab_thickness.append(SE_inst_i.slab_thickness)

        trace_i = go.Scatter(
            x=x_slab_thickness,
            y=y_surface_e,
            mode='markers+lines',
            name=name_i,
            legendgroup=name_i,
            showlegend=True,
            line=dict(
                width=1.5,
                color=color_i,
                dash='dash',
                ),
            marker=dict(
                symbol="square",
                size=8,
                color=color_i,
                line=dict(
                    width=1.0,
                    color="black",
                    ),
                ),

            )
        data.append(trace_i)
        # __|

        # | - Surface Energy (Fitted Bulk) #####################################
        y_surface_e = []; x_slab_thickness = []
        for SE_inst_i in self.new_SurfaceEnergy_instances:
            y_surface_e.append(SE_inst_i.std_surface_e_per_area)
            x_slab_thickness.append(SE_inst_i.slab_thickness)

        trace_i = go.Scatter(
            x=x_slab_thickness,
            y=y_surface_e,
            mode='markers+lines',
            name=name_i,
            # legendgroup=name_i,
            showlegend=False,
            line=dict(
                width=1.5,
                color=color_i,
                dash='solid',
                ),
            marker=dict(
                symbol="circle",
                size=10,
                color=color_i,
                line=dict(
                    width=1.5,
                    color="black",
                    # dash='solid',
                    ),
                ),
            )
        data.append(trace_i)
        # __|

        # | - Average Surface Energy (DFT Bulk) ################################
        ave_surface_energy = self.ave_surface_energy_per_area
        trace_i = go.Scatter(
            x=[0, 30],
            y=[ave_surface_energy, ave_surface_energy],
            name=name_i,
            # legendgroup=name_i,
            mode='lines',
            showlegend=False,
            line=dict(
                width=1.5,
                color=color_i,
                dash='dash',
                ),
            marker=dict(
                symbol="square",
                size=10,
                color=color_i,
                line=dict(
                    width=1,
                    # color='rgb(0, 0, 0)',
                    color=color_i,
                    ),
                ),
            )
        data.append(trace_i)
        # __|

        # | - Average Surface Energy (Fitted Bulk) #############################
        ave_surface_energy = self.new_ave_surface_energy_per_area
        trace_i = go.Scatter(
            x=[0, 30],
            y=[ave_surface_energy, ave_surface_energy],
            name=name_i,
            # legendgroup=name_i,
            mode='lines',
            showlegend=False,
            line=dict(
                width=1.5,
                color=color_i,
                dash='solid',
                ),
            marker=dict(
                symbol="square",
                size=10,
                color=color_i,
                line=dict(
                    width=1,
                    color="black",
                    ),
                ),
            )
        data.append(trace_i)
        # __|

        return(data)
        # __|


    # __| **********************************************************************



# | - METHODS

def surface_energy_2(
    df_i,
    num_points_to_exclude=0,
    trace_title="surface_energy_2",
    color="blue",
    ):
    """
    Calculates the "average" surface energy, but with a bulk_e_per_atom that
    is fitted to a range of slab thicknesses.
    """
    # | - surface_energy_2
    nbpe = num_points_to_exclude

    y_i = df_i["elec_energy"].tolist()
    x_i = df_i["N_atoms"].tolist()

    z = np.polyfit(x_i[nbpe:], y_i[nbpe:], 1)
    bulk_e_per_atom_i = z[0]

    surf_e_2 = df_i["elec_energy"] - df_i["N_atoms"] * bulk_e_per_atom_i
    surf_e_2 = surf_e_2 / df_i.iloc[0]["slab_area"]

    trace = go.Scatter(
        x=df_i["layers"],
        y=surf_e_2.tolist(),
        # mode='markers',
        mode='lines+markers',
        name=trace_title,
        marker=dict(
            symbol="circle",
            size=12,
            color=color,
            line=dict(
                width=2,
                color='rgb(0, 0, 0)'
                )
            ),
        )

    return(trace)
    # __|


def surf_e_4(
    row_i,
    G_H2=0.,
    G_H2O=0.,
    bias=0.,
    pH=0.,
    bulk_e_per_atom=None,
    get_e0_from_row=False,
    norm_mode="area",  # 'area' or 'atoms'
    num_atoms=None,
    units="eV/A^2",  # 'eV/A^2' or 'J/m^2'
    ):
    """
    DEPRECATED SWITCH TO USING CLASS

    Calculate surface energy assuming a water reference state
    and using the computational hydrogen electrode.
    """
    # | - surf_e_4

    # | - Read info from row_i
    # COMBAK This shouldn't be hard coded in
    metal = "Ir"

    # atoms_i = row_i.get("init_atoms", default=None)
    elec_energy = row_i.get("elec_energy", default=0.)
    nonstoich_Os = row_i.get("nonstoich_Os", default=0)
    elems_dict = row_i.get("elem_num_dict", default={})

    if bulk_e_per_atom is None:
        bulk_e_per_atom = row_i.get("bulk_e_per_atom_DFT", default=0.)

    # print("bulk_e_per_atom: ", bulk_e_per_atom)

    N_stoich_in_slab = elems_dict[metal] + \
        elems_dict.get("O", 0) - nonstoich_Os

    nonstoich_Hs = elems_dict.get("H", 0)
    # __|

    # | - Calculate Standard State Surface Energy
    if "surf_e_0" in row_i:
        surf_e_0 = row_i.get("surf_e_0", default=0.)
    else:
        surf_e_0 = 0. + \
            +elec_energy + \
            -N_stoich_in_slab * bulk_e_per_atom + \
            -nonstoich_Os * (G_H2O - G_H2) + \
            -nonstoich_Hs * (G_H2 / 2) + \
            +0.
            # -nonstoich_Os * (G_H2O - G_H2 - 2.518583065) + \
            # -nonstoich_Hs * (G_H2) + \


        if norm_mode == "area":
            norm_term = 2 * row_i["slab_area"]
            surf_e_0 = surf_e_0 / norm_term
        elif norm_mode == "atoms" and num_atoms is not None:
            norm_term = 2 * num_atoms
            surf_e_0 = surf_e_0 / norm_term
    # __|

    # | - Calculate V, pH Dependant Surface Energy
    slope = 2 * nonstoich_Os - nonstoich_Hs

    surf_e = 0. + \
        +surf_e_0 + \
        +(slope * 0.0591 * pH) / norm_term + \
        -(slope * bias) / norm_term + \
        +0.

        # +(slope * 0.0591 * pH) / (2 * row_i["slab_area"]) + \
        # -(slope * bias) / (2 * row_i["slab_area"]) + \

#     surf_e = surf_e / (2 * row_i["slab_area"])
    # __|


    # | - Unit conversion
    units="eV/A^2",  # 'eV/A^2' or 'J/m^2'

    if norm_mode == "area":
        if units == "eV/A^2":
            pass
        elif units == "J/m^2":
            # Convert eV/A^2 to J/m^2
            # (1E10 A/m) ^ 2 * (1.6022E-19 J/eV) = 16.022
            ev_A2__to__J_m2 = 16.022
            surf_e = surf_e * ev_A2__to__J_m2
    # __|

    return(surf_e)
    # __|

# __|
