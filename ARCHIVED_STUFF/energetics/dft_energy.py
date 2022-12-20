#!/usr/bin/env python

"""Encapsulate DFT energy methods.

Author: Raul A. Flores

Development:
    NOTE Take into account BEEF ensemble of energies
"""

# | - Import Modules
import os
import sys
# __|


class Energy(object):
    """
    """

    # | - Energy ***************************************************************
    def __init__(self,
        gibbs_e=None,
        internal_e=None,
        enthalpy_e=None,
        helmholtz_e=None,
        electronic_e=None,
        zero_point_e=None,
        Cv_trans_term=None,
        Cv_rot_term=None,
        Cv_vib_term=None,
        Cv_to_Cp=None,
        entropy_term=None,
        PV_term=None,
        # main_energy="gibbs",
        ):
        """Initialize system energy quantities.

        TODO: Wrap all math operations around my custom methods that check for
        None exeptions

        Args:
            electronic_e:
            zero_point_e:
            trans_e:
            rot_e:
            vib_e:
        """
        # | - __init__
        self.gibbs_e = gibbs_e
        self.internal_e = internal_e
        self.enthalpy_e = enthalpy_e
        self.helmholtz_e = helmholtz_e

        self.electronic_e = electronic_e
        self.zero_point_e = zero_point_e

        self.Cv_trans_term = Cv_trans_term
        self.Cv_rot_term = Cv_rot_term
        self.Cv_vib_term = Cv_vib_term
        self.Cv_to_Cp = Cv_to_Cp

        self.entropy_term = entropy_term
        self.PV_term = PV_term

        if self.internal_e is None:
            self.internal_e = self.calc_internal_energy()

        if self.enthalpy_e is None:
            self.enthalpy_e = self.calc_enthalpy_energy()

        if self.gibbs_e is None:
            self.gibbs_e = self.calc_gibbs_free_energy()
        # __|


    def __str__(self):
        """
        """
        # | - __str__
        mess = "Gibbs Free Energy: " + str(self.gibbs_e) + "\n"
        mess += "Electronic Energy: " + str(self.electronic_e) + "\n"

        return(mess)
        # __|

    def __repr__(self):
        """Representation of instance when printed to stdout."""
        # | - __repr__
        mess = "Electronic Energy: " + str(self.electronic_e) + "\n"
        mess += "Enthalpy Energy: " + str(self.enthalpy_e) + "\n"
        mess += "Gibbs Free Energy: " + str(self.gibbs_e) + "\n"

        return(mess)
        # __|

    def __sub__(self, other):
        """
        """
        # | - __sub__

        def subtract_mine(a, b):
            """Return a - b.

            Return None if either of the terms are None

            Args:
                a:
                b:
            """
            # | - subtract_mine
            if a is None or b is None:
                return(None)
            else:
                out = a - b

                return(out)
            # __|


        if isinstance(other, Energy):
            # print("Divisor is Energy class instance!!!")

            electronic_e_new = subtract_mine(
                self.electronic_e,
                other.electronic_e)

            enthalpy_e_new = subtract_mine(
                self.enthalpy_e,
                other.enthalpy_e)

            gibbs_e_new = subtract_mine(
                self.gibbs_e,
                other.gibbs_e)

            # electronic_e_new = self.electronic_e - other.electronic_e
            # enthalpy_e_new = self.enthalpy_e - other.enthalpy_e
            # gibbs_e_new = self.gibbs_e - other.gibbs_e


        elif isinstance(other, float) or isinstance(other, int):


            electronic_e_new = subtract_mine(self.electronic_e, float(other))
            enthalpy_e_new = subtract_mine(self.enthalpy_e, float(other))
            gibbs_e_new = subtract_mine(self.gibbs_e, float(other))


            # electronic_e_new = self.electronic_e - float(other)
            # enthalpy_e_new = self.enthalpy_e - float(other)
            # gibbs_e_new = self.gibbs_e - float(other)

        else:
            print("type:")
            print(type(other))
            raise NotImplementedError(
                "Haven't figured out how to subract to the Energy class yet"
                )

        E_dict={
            "enthalpy_e": enthalpy_e_new,
            "gibbs_e": gibbs_e_new,
            "electronic_e": electronic_e_new,
            }

        out_Energy = Energy(**E_dict)

        return(out_Energy)
        # __|

    def __add__(self, other):
        """
        """
        # | - __add__
        return (self.gibbs_e + other.gibbs_e)
        # __|

    def __truediv__(self, other):
        """
        """
        # | - __truediv__

        if isinstance(other, Energy):
            electronic_e_new = self.electronic_e / other.electronic_e
            enthalpy_e_new = self.enthalpy_e / other.enthalpy_e
            gibbs_e_new = self.gibbs_e / other.gibbs_e

        elif isinstance(other, float) or isinstance(other, int):
            electronic_e_new = self.electronic_e / float(other)
            enthalpy_e_new = self.enthalpy_e / float(other)
            gibbs_e_new = self.gibbs_e / float(other)

        else:
            raise NotImplementedError(
                "Haven't figured out how to divide my Energy class by that yet"
                )

        E_dict={
            "electronic_e": electronic_e_new,
            "enthalpy_e": enthalpy_e_new,
            "gibbs_e": gibbs_e_new,
            }

        out_Energy = Energy(**E_dict)

        return(out_Energy)
        # __|

    def __floordiv__(self, other):
        """
        """
        # | - __floordiv__
        print("__floordiv__")
        return (self.gibbs_e / other.gibbs_e)
        # __|



    @staticmethod
    def add_entries(entries_list):
        """Adds entries in entries_list.

        Converts NoneType to 0

        Args:
            entries_list:
        """
        # | - add_entries
        sum_tot = 0.
        for entry in entries_list:
            if entry is None:
                summand = 0.
            else:
                summand = entry
            sum_tot += summand

        return(sum_tot)
        # __|

    def calc_internal_energy(self):
        """Calculate internal energy.

        Args:

        """
        # | - internal_energy
        energy_list = [
            self.electronic_e,
            self.zero_point_e,
            self.Cv_trans_term,
            self.Cv_rot_term,
            self.Cv_vib_term,
            self.Cv_to_Cp,
            ]

        internal_energy = self.add_entries(energy_list)

        return(internal_energy)
        # __|

    def calc_enthalpy_energy(self):
        """
        """
        # | - calc_enthalpy_energy
        energy_list = [
            self.internal_e,
            self.PV_term,
            ]

        enthalpy_e = self.add_entries(energy_list)

        return(enthalpy_e)
        # __|

    def calc_gibbs_free_energy(self):
        """
        """
        # | - calc_gibbs_free_energy
        if self.entropy_term is not None:
            entropy_term = -self.entropy_term
        else:
            entropy_term = None

        energy_list = [
            self.enthalpy_e,
            entropy_term,
            ]

        gibbs_e = self.add_entries(energy_list)

        return(gibbs_e)
        # __|

    # __| **********************************************************************


class Element_Refs():
    """docstring for [object Object].

    /u/if/flores12/02_ref_energies/h2/_3

    Notes:

        * Water formation energy:
            G_Liq: | -237.14 kJ/mol  --> -2.457784 eV (-4.915567)
            G_Gas: | -228.61 kJ/mol  --> -2.369376 eV (-4.738753)
            --------------------------------------------------------
            H_Liq: | −285.8 kJ/mol   --> -2.96211 eV  (-5.92422)
            H_Gas: | −241.818 kJ/mol --> -2.506268 eV (-5.012535)

    """

    # | - Element_Refs *********************************************************

    def __init__(self,
         H2O_dict={
            # "gibbs_e": -476.54410902835434,
            "gibbs_e": -476.630481653821,
            "electronic_e": -476.6430744328504,
            },

        H2_dict={
            # "gibbs_e": -32.956123446117,
            "gibbs_e": -32.956123446117,
            "electronic_e": -32.9203465588979,
            },

        O2_dict={
            "gibbs_e": -883.1904356,
            "electronic_e": -882.7384356,
            },

        # COMBAK I only have the electronic energy for now
        N2_dict={
            "gibbs_e": -553.638294974,
            "electronic_e": -553.638294974,
            },


        # * Water formation energy:
        #     G_Liq: | -237.14 kJ/mol  --> -2.457784 eV (-4.915567)
        #     G_Gas: | -228.61 kJ/mol  --> -2.369376 eV (-4.738753)
        #     --------------------------------------------------------
        #     H_Liq: | −285.8 kJ/mol   --> -2.96211 eV  (-5.92422)
        #     H_Gas: | −241.818 kJ/mol --> -2.506268 eV (-5.012535)
        #

        # H2O Formation Energy from H2 and O2 *********************************
        H2O_form_e_dict={
            # "gibbs_e": -2.4583,
            # "enthalpy_e": -2.96211,

            "gibbs_e": -2.4583,  # Liquid phase Gibbs
            "enthalpy_e": -2.506268,  # Gas phase enthalpy
            },

        # H2O_form_gibbs=-2.4583,

        oxygen_ref="H2O",  # "H2O" or "O2"
        hydrogen_ref="H2",

        # reference_states_dict={
        #     "H": "H2",
        #     "O": "H2O",
        #     }

        ):
        """Initialize instance.

        # Directories of gas phase DFT jobs (QE):
        H2O, H2, and O2 energies were obtained from:
        /scratch/users/flores12/gas_phase_molec/BEEF-vdW
        H2:
        /scratch/users/flores12/gas_phase_molec/BEEF-vdW/h2/_4
        H2O:
        /scratch/users/flores12/gas_phase_molec/BEEF-vdW/h2o/_3
        O2:
        /scratch/users/flores12/gas_phase_molec/BEEF-vdW/o2/_3

        Args:
            G_H2O:
            G_H2:
            G_O2:
            H2O_form_gibbs:
            oxygen_ref and hydrogen_ref:
                States of 0 energy, reference for each element type
                Conventionally for ORR/OER I use a H2, H2O reference state,
                such that H2O has 0 energy
        """
        # | - __init__
        self.En_H2O = Energy(**H2O_dict)
        self.En_H2 = Energy(**H2_dict)
        self.En_O2 = Energy(**O2_dict)
        self.En_N2 = Energy(**N2_dict)

        self.oxygen_ref = oxygen_ref
        self.hydrogen_ref = hydrogen_ref

        # self.H2O_form_gibbs = Energy(gibbs_e=H2O_form_gibbs)
        self.H2O_form_gibbs = Energy(**H2O_form_e_dict)

        self.E_O_ref, self.E_H_ref = self.calc_ref_energies()
        # __|

    def calc_ref_energies(self):
        """
        """
        # | - calc_ref_energies

        if self.hydrogen_ref == "H2":
            # hyd_ref = self.En_H2 / 2.
            hyd_ref = self.En_H2 / 2.
        else:
            print("Non H2 H reference state! Do more work here")


        if self.oxygen_ref == "H2O":
            oxy_ref = self.En_H2O - self.En_H2

        elif self.oxygen_ref == "O2":
            oxy_ref = self.En_H2O - self.En_H2 - self.H2O_form_gibbs

            # tmp = oxy_ref - self.H2O_form_gibbs
            # print(tmp)

        return(oxy_ref, hyd_ref)
        # __|

    # __| **********************************************************************
