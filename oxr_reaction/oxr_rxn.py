#!/usr/bin/env python

"""ORR energetics classes and methods.

Author: Raul A. Flores
"""

# | - IMPORT MODULES
import numpy as np
import pandas as pd

from sklearn.linear_model import LinearRegression

import plotly.graph_objs as go

pd.options.mode.chained_assignment = None

from oxr_reaction.oxr_series import ORR_Free_E_Series
from oxr_reaction.adsorbate_scaling import lim_U_i
# from oxr_reaction.oxr_rxn import ORR_Free_E_Plot
# __|


class ORR_Free_E_Plot:
    """ORR free energy diagram class.

    ACTUALLY THIS IS GOING TO BE A GENERAL ORR/OER CLASS NOW!!!!!!!!!!!!!!!!!!!

    Development Notes:
        TODO Rename this base class to OXR_BASE or something
        TODO Should we consider the case where the bulk energy is not 0, and
        we have to normalize all of the species energies by it?
    """

    # | - ORR_Free_E_Plot ******************************************************

    def __init__(self,
        free_energy_df=None,
        ORR_Free_E_series_list=None,
        state_title="adsorbate",
        free_e_title="ads_e",
        num_states=5,
        smart_format=None,
        bias=0.,
        color_list=None,
        hover_text_col=None,
        rxn_type="ORR",  # ORR and OER

        # system_properties=None,
        # opt_name=None,
        # properties=None,
        # i_cnt=0,
        # plot_mode="all",
        # smart_format=None,
        # Plotting ************************************************************
        # show_H_e_pairs_annotations=True,
        # show_legend=True,
        ):
        """
        Input variables to class instance.

        Args:
            free_energy_df:
                Pandas dataframe containing the adsorbates as rows
                Required columns, adsorbate, free energy
            system_properties:
            state_title:
            free_e_title:

            rxn_type:
                ORR or OER
        """
        # | - __init__

        # | - Setting Instance Attributes
        self.fe_df = free_energy_df
        # self.sys_props = system_properties
        self.state_title = state_title
        self.fe_title = free_e_title
        self.num_states = num_states

        # self.rxn_mech_states = ["bulk", "ooh", "o", "oh", "bulk"]
        # self.ideal_energy = [4.92, 3.69, 2.46, 1.23, 0]  #COMBAK Is this being used anymore

        self.bias = bias
        self.color_list = color_list
        self.hover_text_col = hover_text_col
        self.smart_format = smart_format

        # self.show_H_e_pairs_annotations = show_H_e_pairs_annotations
        # self.show_legend = show_legend

        # ***********************************
        # COMBAK, moving to FED class, remove later
        self.plot_states_sep = 0.3
        self.plot_states_width = 1.

        self.rxn_type = rxn_type
        # __|


        if self.rxn_type == "ORR":
            self.rxn_mech_states = ["bulk", "ooh", "o", "oh", "bulk"]
            self.ideal_energy = [4.92, 3.69, 2.46, 1.23, 0]

        elif self.rxn_type == "OER":
            self.rxn_mech_states = ["bulk", "oh", "o", "ooh", "bulk"]
            self.ideal_energy = [0, 1.23, 2.46, 3.69, 4.92]

        # self.rxn_x_coord_array = self.create_rxn_coord_array(
        #     self.num_states,
        #     spacing=self.plot_states_sep,
        #     step_size=self.plot_states_width,
        #     )
        #
        # self.mid_state_x_array = self.create_mid_state_x_array()
        # # x_array_data = self.rxn_x_coord_array

        if ORR_Free_E_series_list is None:
            self.series_list = []
        else:
            self.series_list = ORR_Free_E_series_list

        # | - __old__
        # if free_energy_df is not None:
        #     self.add_bulk_entry()
        #     self.fill_missing_data()
        #
        #     self.num_of_states = len(self.fe_df) + 1  # bulk, OOH, O, OH, bulk
        #     self.energy_lst = self.rxn_energy_lst()
        #     self.num_of_elec = range(self.num_of_states)[::-1]
        #     self.overpotential = self.calc_overpotential()[0]
        #     self.limiting_step = self.calc_overpotential()[1]
        #     # self.ideal_energy = [4.92, 3.69, 2.46, 1.23, 0]
        #     self.energy_lst_h2o2 = self.rxn_energy_lst_h2o2()
        #     self.overpotential_h2o2 = self.calc_overpotential_h2o2()
        # __|

        # __|

    def __create_series_name__(self, series_i):
        """
        """
        # | - __create_series_name__

        if series_i.properties is not None:
            name_i = ""
            for key, value in series_i.properties.items():

                # TODO | This is an old feature, get rid of it
                if key == "coverage":
                    continue

                name_i += str(key) + ": " + str(value) + " | "

        else:
            name_i = ""

        return(name_i)
        # __|

    def __create_smart_format_dict__(self, property_dict, smart_format_dict):
        """Create smart format dictionary.

        Args:
            property_dict:
            smart_format_dict:
        """
        # | - __create_smart_format_dict__
        if property_dict is None:
            return({})

        format_dict = {}
        for key_i, value_i in property_dict.items():
            for format_i in smart_format_dict:
                if list(format_i[0])[0] == key_i:
                    if list(format_i[0].values())[0] == value_i:
                        format_dict.update(format_i[1])

        return(format_dict)
        # __|

    def add_series(self,
        fe_df,
        plot_mode="all",
        name_i=None,
        group=None,
        opt_name=None,
        smart_format=True,
        format_dict=None,
        overpotential_type="ORR",
        system_properties=None,
        property_key_list=None,
        add_overpot=True,
        color=None,
        overpotential_given=False,
        ):
        """Add ORR_Free_E_Series instance to ORR_Free_E_Plot.series_list.

        Note: It would be much better to simply take all of the
        ORR_Free_E_Series arguments as a **kwargs term.

        Args:
            TEMP
        """
        # | - add_series
        if smart_format:
            smart_format_i = self.smart_format
        else:
            smart_format_i = None

        ORR_Series = ORR_Free_E_Series(
            free_energy_df=fe_df,
            properties=system_properties,
            property_key_list=property_key_list,
            state_title=self.state_title,
            free_e_title=self.fe_title,
            group=group,
            bias=self.bias,
            # rxn_x_coord_array=self.rxn_x_coord_array,
            name_i=name_i,
            opt_name=opt_name,  # #######

            color_list=self.color_list,
            color=color,
            hover_text_col=self.hover_text_col,
            plot_mode=plot_mode,  # ##########
            smart_format=smart_format_i,
            format_dict=format_dict,
            add_overpot=add_overpot,
            rxn_type=self.rxn_type,
            overpotential_given=overpotential_given,

            # properties=opt_name,
            # overpotential_type=self.rxn_type,
            # i_cnt=0,  # ##########
            )

        self.series_list.append(ORR_Series)
        # __|

    # __| **********************************************************************
