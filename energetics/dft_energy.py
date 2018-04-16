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
        gibbs_e = None,
        internal_e = None,
        enthalpy_e = None,
        helmholtz_e= None,

        electronic_e = None,
        zero_point_e = None,

        Cv_trans_term = None,
        Cv_rot_term = None,
        Cv_vib_term = None,
        Cv_to_Cp = None,

        entropy_term = None,
        PV_term = None,
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

    def __sub__(self, other):
        return self.gibbs_e - other.gibbs_e

    def __add__(self, other):
        return self.gibbs_e + other.gibbs_e


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
            "gibbs_e": -476.843794302748,
            "electronic_e": -476.643074433,
            },

        H2_dict={
            "gibbs_e": -32.9563981542851,
            "electronic_e": -32.920117439,
            },

        O2_dict={
            "gibbs_e": -883.190570481887,
            "electronic_e": -882.7384356,
            },

        # G_H2O=-476.843794302748,
        # G_H2=-32.9563981542851,
        # G_O2=-883.190570481887,

        H2O_form_gibbs=-2.4583,

        oxygen_ref="H2O",  # "H2O" or "O2"
        hydrogen_ref="H2",
        ):
        """

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
        # print(self.En_H2O)
        # self.En_H2O = Energy(gibbs_e=G_H2O)

        self.En_H2 = Energy(gibbs_e=G_H2)
        self.En_O2 = Energy(gibbs_e=G_O2)

        self.H2O_form_gibbs = Energy(gibbs_e=H2O_form_gibbs)

        self.E_O_ref = self.calc_ref_energies()
        #__|

    def calc_ref_energies(self):
        """
        """
        #| - tmp
        # pass

        tmp = self.En_H2O - self.En_H2

        print(tmp)
        #__|

    #__| **********************************************************************
