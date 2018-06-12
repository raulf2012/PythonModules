"""Encapsulate DFT energy methods.

Author: Raul A. Flores

Development:
    NOTE Take into account BEEF ensemble of energies
"""

#| - Import Modules
import os
import sys
#__|


class Energy(object):
    """
    """

    #| - Energy ***************************************************************
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

        Args:
            electronic_e:
            zero_point_e:
            trans_e:
            rot_e:
            vib_e:
        """
        #| - __init__
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
        #__|


    def __str__(self):
        """
        """
        #| - __str__
        mess = "Gibbs Free Energy: " + str(self.gibbs_e) + "\n"
        mess += "Electronic Energy: " + str(self.electronic_e) + "\n"

        return(mess)
        #__|

    def __repr__(self):
        """
        """
        #| - __repr__
        mess = "Gibbs Free Energy: " + str(self.gibbs_e) + "\n"
        mess += "Electronic Energy: " + str(self.electronic_e) + "\n"

        return(mess)
        #__|

    def __sub__(self, other):
        """
        """
        #| - __sub__
        # return (self.gibbs_e - other.gibbs_e)

        if isinstance(other, Energy):
            print("Divisor is Energy class instance!!!")

            # lst = ["gibbs_e", "enthalpy", "electronic_e"]
            gibbs_e_new = self.gibbs_e - other.gibbs_e
            electronic_e_new = self.electronic_e - other.electronic_e


        elif isinstance(other, float) or isinstance(other, int):
            print("Divisor is integer or float")
            print(type(other))

            gibbs_e_new = self.gibbs_e - float(other)
            electronic_e_new = self.electronic_e - float(other)

        else:
            print("type:")
            print(type(other))
            raise NotImplementedError("Haven't figured out how to subract to the Energy class yet")

        E_dict={
            "gibbs_e": gibbs_e_new,
            "electronic_e": electronic_e_new,
            }

        out_Energy = Energy(**E_dict)

        return(out_Energy)
        #__|

    def __add__(self, other):
        """
        """
        #| - __add__
        return (self.gibbs_e + other.gibbs_e)
        #__|

    def __truediv__(self, other):
        """
        """
        #| - __truediv__
        # print("__truediv__")

        if isinstance(other, Energy):
            # print("Divisor is Energy class instance!!!")
            # lst = ["gibbs_e", "enthalpy", "electronic_e"]
            gibbs_e_new = self.gibbs_e / other.gibbs_e
            electronic_e_new = self.electronic_e / other.electronic_e


        elif isinstance(other, float) or isinstance(other, int):
            # print("Divisor is integer or float")
            # print(type(other))

            gibbs_e_new = self.gibbs_e / float(other)
            electronic_e_new = self.electronic_e / float(other)

        else:
            # print("Dividand type:")
            # print(type(other))
            raise NotImplementedError("Haven't figured out how to divide my Energy class by that yet")

        E_dict={
            "gibbs_e": gibbs_e_new,
            "electronic_e": electronic_e_new,
            }

        out_Energy = Energy(**E_dict)

        return(out_Energy)
        #__|

    def __floordiv__(self, other):
        """
        """
        #| - __floordiv__
        print("__floordiv__")
        return (self.gibbs_e / other.gibbs_e)
        #__|



    @staticmethod
    def add_entries(entries_list):
        """Adds entries in entries_list.

        Converts NoneType to 0

        Args:
            entries_list:
        """
        #| - add_entries
        sum_tot = 0.
        for entry in entries_list:
            if entry is None:
                summand = 0.
            else:
                summand = entry
            sum_tot += summand

        return(sum_tot)
        #__|

    def calc_internal_energy(self):
        """Calculate internal energy.

        Args:

        """
        #| - internal_energy
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
        #__|

    def calc_enthalpy_energy(self):
        """
        """
        #| - calc_enthalpy_energy
        energy_list = [
            self.internal_e,
            self.PV_term,
            ]

        enthalpy_e = self.add_entries(energy_list)

        return(enthalpy_e)
        #__|

    def calc_gibbs_free_energy(self):
        """
        """
        #| - calc_gibbs_free_energy
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
        #__|

    #__| **********************************************************************

class Element_Refs():
    """docstring for [object Object].

    /u/if/flores12/02_ref_energies/h2/_3

    Notes:

    """

    #| - Element_Refs *********************************************************

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
            "gibbs_e" -553.638294974:,
            "electronic_e": -553.638294974,
            },

        H2O_form_gibbs=-2.4583,

        oxygen_ref="H2O",  # "H2O" or "O2"
        hydrogen_ref="H2",
        ):
        """Initialize instance.

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
        """
        #| - __init__
        self.En_H2O = Energy(**H2O_dict)
        self.En_H2 = Energy(**H2_dict)
        self.En_O2 = Energy(**O2_dict)
        self.En_N2 = ENergy(**N2_dict)

        self.oxygen_ref = oxygen_ref

        self.H2O_form_gibbs = Energy(gibbs_e=H2O_form_gibbs)

        self.E_O_ref, self.E_H_ref = self.calc_ref_energies()
        #__|

    def calc_ref_energies(self):
        """
        """
        #| - calc_ref_energies
        # hyd_ref = self.En_H2 / 2.
        hyd_ref = self.En_H2 / 2.

        if self.oxygen_ref == "H2O":
            # print("KJFKDSJKFJDSLJFLISDJFIJSDKFJKDS -------------------------")
            oxy_ref = self.En_H2O - self.En_H2
            # print("KJFKDSJKFJDSLJFLISDJFIJSDKFJKDS -------------------------")

        return(oxy_ref, hyd_ref)
        #__|

    #__| **********************************************************************
