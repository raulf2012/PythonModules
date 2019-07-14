#!/usr/bin/env python

"""ORR energetics classes and methods.

Author: Raul A. Flores
"""

#| - IMPORT MODULES
import copy

import numpy as np
import pandas as pd

from plotly.graph_objs import Scatter
pd.options.mode.chained_assignment = None
#__|

class ORR_Free_E_Series():
    """ORR free energy diagram series class.

    Still a work in progress

    Take a 2nd look at how the bulk/bare state is being taken care of

    # TODO | The color_list approach is not good
    """

    #| - ORR_Free_E_Series ****************************************************

    def __init__(self,
        free_energy_df=None,
        state_title="adsorbate",
        free_e_title="ads_e",
        bias=0.,
        group=None,
        name_i=None,
        opt_name=None,
        properties=None,
        property_key_list=None,
        color_list=None,
        color=None,
        hover_text_col=None,
        plot_mode="all",
        smart_format=None,
        format_dict=None,
        add_overpot=True,
        rxn_type="ORR",
        fill_missing_data_w_scaling=True,
        overpotential_given=False,
        ):
        """
        Input variables to class instance.

        TODO: Make sure that format_dict is used if it's supplied by user

        Args:
            free_energy_df:
                Pandas dataframe containing the adsorbates as rows
                Required columns, adsorbate, free energy
            state_title:
            free_e_title:
            bias:
            opt_name:
                Optional name
                string to append to the beginning of the plot series name
            properties:
            color_list:
            hover_text_col:
            plot_mode:
            smart_format:
            format_dict:
            overpotential_given:
                If true then the overpotential is given explicitely in the
                input dataframe with name "overpotential"
        """
        #| - __init__

        #| - Setting Instance Attributes
        self.fe_df = free_energy_df
        # self.sys_props = system_properties

        self.group = group
        self.state_title = state_title
        self.fe_title = free_e_title

        self.bias = bias

        # self.rxn_x_coord_array = rxn_x_coord_array
        # print(self.rxn_x_coord_array)  # TEMP | PRINT

        self.name_i = name_i
        self.opt_name = opt_name

        self.properties = properties

        self.color_list = color_list
        self.color = color
        self.hover_text_col = hover_text_col
        self.plot_mode = plot_mode
        self.smart_format = smart_format
        self.format_dict = format_dict

        # self.overpotential_type = overpotential_type
        self.rxn_type = rxn_type

        self.add_overpot = add_overpot
        self.overpotential_given = overpotential_given

        # self.i_cnt = i_cnt
        #__|

        if self.rxn_type == "ORR":
            self.rxn_mech_states = ["bulk", "ooh", "o", "oh", "bulk"]
            self.ideal_energy = [4.92, 3.69, 2.46, 1.23, 0]

        elif self.rxn_type == "OER":
            self.rxn_mech_states = ["bulk", "oh", "o", "ooh", "bulk"]
            self.ideal_energy = [0, 1.23, 2.46, 3.69, 4.92]

        self.property_key_list = property_key_list

        # Doing this with a class method instead of in the analysis script
        if properties is None:
            self.properties = self.__create_property_dict__()

        if free_energy_df is not None:
            self.add_bulk_entry()
            self.fill_missing_data()

            # TEMP | Active Development
            self.num_of_states = len(self.fe_df) + 1  # bulk, OOH, O, OH, bulk
            self.num_of_states_new = self.__num_of_states__()

            self.energy_lst = self.rxn_energy_lst()
            # self.energy_lst_new = self.rxn_energy_lst_new()

            self.num_of_elec = range(self.num_of_states)[::-1]
            self.overpotential = self.calc_overpotential()[0]
            self.limiting_step = self.calc_overpotential()[1]
            # self.ideal_energy = [4.92, 3.69, 2.46, 1.23, 0]
            self.energy_lst_h2o2 = self.rxn_energy_lst_h2o2()
            self.overpotential_h2o2 = self.calc_overpotential_h2o2()
            self.overpotential_OER = self.calc_overpotential_OER()[0]

            self.energy_states_dict = self.__energy_states_dict__()

            if fill_missing_data_w_scaling:
                self.__fill_nan_values__()

            # TODO | Put this outside of the if-statement
            self.series_name = self.__series_plot_name__(
                opt_name=self.opt_name,
                properties=self.properties,
                overpotential_type=self.rxn_type,
                add_overpot=self.add_overpot,
                )
        #__|

    def __fill_nan_values__(self):
        """Fill nan adsorption energy values from scaling relations."""
        #| - tmp
        energy_dict = self.energy_states_dict
        if True in list(np.isnan(list(energy_dict.values()))):

            print("There is a nan in the energy dict!!!")

            if np.isnan(energy_dict["ooh"]):
                if not np.isnan(energy_dict["oh"]):
                    print("*OOH energy set by *OH and standard scaling")
                    ooh_new = 1 * energy_dict["oh"] + 3.2
                    energy_dict["ooh"] = ooh_new

            if np.isnan(energy_dict["o"]):
                if not np.isnan(energy_dict["oh"]):
                    print("*O energy set by *OH and standard scaling")
                    o_new = 2 * energy_dict["oh"] + 0.
                    energy_dict["o"] = o_new

            if np.isnan(energy_dict["oh"]):
                if not np.isnan(energy_dict["ooh"]):
                    print("*OH energy set by *OOH and standard scaling")
                    oh_new = energy_dict["ooh"] - 3.2
                    energy_dict["oh"] = oh_new

        self.energy_states_dict = energy_dict
        #__|

    def __num_of_states__(self):
        """Return number of unique states.

        Looks at the uniqe number of entries in the 'adsorbate' column of the
        data frame. The correct number of states for the OER and/or ORR
        reaction are 4, only 2 states are needed for the peroxide reaction.
        """
        #| - __num_of_states
        df_i = self.fe_df

        num_of_states = len(set(df_i["adsorbate"].tolist()))

        err_mess = "There are not enough unique calcs (less than 4)"
        assert num_of_states >= 4, err_mess

        return(num_of_states)
        #__|

    def __energy_states_dict__(self):
        """
        """
        #| - __energy_states_dict__

        energy_lst = self.energy_lst
        rxn_mech_states = self.rxn_mech_states

        energy_states_dict = dict(zip(rxn_mech_states, energy_lst))
        energy_states_dict.pop("bulk", None)

        return(energy_states_dict)
        # energy_states_dict
        #__|

    def __create_property_dict__(self):
        """
        """
        #| - __create_property_dict__
        df_i = self.fe_df

        def all_same_val(df_i, prop_i, val_1):
            """
            """
            #| - all_same_val
            out_list = []
            for i in df_i[prop_i].tolist():
                if i == val_1:
                    out_list.append(True)
                else:
                    out_list.append(False)

            out_i = all(out_list)
            return(out_i)

            # [True if i == val_1 else False for i in
            # df_i[prop_i].tolist()],
            #__|

        if self.property_key_list is not None:
            prop_dict_i = {}
            for prop_i in self.property_key_list:
                val_1 = df_i[prop_i].tolist()[0]

                all_same_value = all_same_val(df_i, prop_i, val_1)
                # all_same_value = all(
                #     [True if i == val_1 else False for i in
                # df_i[prop_i].tolist()],
                #     )

                if all_same_value:
                    prop_dict_i[prop_i] = str(val_1)
                else:
                    prop_dict_i[prop_i] = str(None)

            return(prop_dict_i)
        else:
            return({})
        #__|

    def add_bulk_entry(self,
        bulk_e=0.0,
        ):
        """
        Append a row entry to data frame corresponding to bulk state.

        Args:
            bulk_e:
        """
        #| - add_bulk_entry
        df = self.fe_df
        bulk_df = pd.DataFrame([{
            "adsorbate": "bulk",
            "ads_e": bulk_e,
            }])

        # TEMP
        # df = df.append(bulk_df, ignore_index=True)
        df = df.append(bulk_df, ignore_index=True, sort=True)

        self.fe_df = df
        #__|

    def rxn_energy_lst_h2o2(self):
        """Construct energy list of h2o2 FED."""
        #| - rxn_energy_lst_h2o2
        # h2o2_e = 3.52

        df = self.fe_df

        free_energy_list = []
        for index, row in df.iterrows():
            if row["adsorbate"] == "bulk" or row["adsorbate"] == "ooh":
                free_energy_list.append(row["ads_e"])

        # TODO | Make this parse the reaction array instead of reaction list
        # Checking length of energy list
        if len(free_energy_list) != 2:
            # raise ValueError("Not the correct # of steps for H2O2")
            print("Not the correct # of steps for H2O2")

        free_energy_list[0] += 4.92
        free_energy_list.append(3.52)

        return(free_energy_list)
        #__|

    def property_list(self, column_name):
        """General method to create a list from a column in the dataframe.

        The length of the list will correspond to the steps in the ORR
        mechanism.

        Args:
            column_name:
        """
        #| - property_list
        df = self.fe_df

        property_list = []
        for state in self.rxn_mech_states:
            tmp = df.loc[df[self.state_title] == state]
            tmp1 = tmp.iloc[0][column_name]
            property_list.append(tmp1)

        # free_energy_list[0] += 4.92

        return(property_list)
        #__|

    def fill_missing_data(self):
        """
        """
        #| - fill_missing_data
        df = self.fe_df
        df_missing_data = pd.DataFrame()
        for state in self.rxn_mech_states:
            df_state = df.loc[df[self.state_title] == state]

            #| - If df is missing state fill in row with NaN for energy
            if df_state.empty:
                df_state = pd.DataFrame([{
                    self.state_title: state,
                    self.fe_title: np.nan,
                    }])
                df_missing_data = df_missing_data.append(df_state)
            #__|

        self.fe_df = self.fe_df.append(df_missing_data, sort=True)
        #__|

    def rxn_energy_lst(self):
        """List corresponding to the steps of ORR.

        (1. O2, 2. *OOH, 3. *O, 4. *OH, 5. 2H2O)
        """
        #| - rxn_energy_lst
        df = self.fe_df

        free_energy_list = []
        for state in self.rxn_mech_states:
            df_state = df.loc[df[self.state_title] == state]

            # print(df_state)

            #| - __old__
            # Not sure what this was trying to accomplish
            # if len(df_state) == 2:
            #     state_energy_list = []
            #     for j_cnt, row_j in df_state.iterrows():
            #         energy_j = row_j[self.fe_title]
            #         state_energy_list.append(energy_j)
            #__|

            #| - If df is missing state fill in row with NaN for energy
            if df_state.empty:
                df_state = pd.DataFrame([{
                    self.state_title: state,
                    self.fe_title: np.nan,
                    }])
            #__|

            # This just takes the first species
            # If you feed a df with more than one entry per species, then
            # this will stupidly choose the first one
            state_energy_1 = df_state.iloc[0][self.fe_title]

            #| - __old__
            # if type(state_energy_1) != float:
            #     print(type(state_energy_1))
            #     print("DSKFJKLSDJFSjk_d--_d-d-_D_D_d-d-d-d-d___D_D_D_")
            # print(
            #     "state_energy_1: ",
            #     str(state_energy_1),
            #     )
            #
            # print(
            #     "type: ",
            #     str(type(state_energy_1))
            #     )
            #
            # print(isinstance(state_energy_1, np.float))
            # print(float(state_energy_1))
            # print(type(float(state_energy_1)))
            # print(np.isnan(state_energy_1))
            # if isinstance(state_energy_1, np.float) is False:
            #     print("lkjfksjd")
            #__|

            free_energy_list.append(state_energy_1)

        if self.rxn_type == "ORR":
            free_energy_list[0] += 4.92
        elif self.rxn_type == "OER":
            free_energy_list[-1] += 4.92
        else:
            free_energy_list[0] += 4.92

        return(free_energy_list)
        #__|

    def rxn_energy_lst_new(self):
        """
        """
        #| - rxn_energy_lst_new
        df = self.fe_df
        free_energy_list = []
        for state in self.rxn_mech_states:
            df_state = df.loc[df[self.state_title] == state]

            #| - If df is missing state fill in row with NaN for energy
            if df_state.empty:
                df_state = pd.DataFrame([{
                    self.state_title: state,
                    self.fe_title: np.nan,
                    }])
            #__|

            state_energy_list = []
            for j_cnt, row_j in df_state.iterrows():
                energy_j = row_j[self.fe_title]
                state_energy_list.append(energy_j)

            free_energy_list.append(state_energy_list)

            # tmp1 = df_state.iloc[0][self.fe_title]
            # free_energy_list.append(tmp1)

        if self.rxn_type == "ORR":
            free_energy_list_0_new = [i + 4.92 for i in free_energy_list[0]]
            free_energy_list[0] = free_energy_list_0_new

            # free_energy_list[0] += 4.92

        elif self.rxn_type == "OER":
            # free_energy_list[-1] += 4.92
            free_energy_list_new = [i + 4.92 for i in free_energy_list[-1]]
            free_energy_list[-1] = free_energy_list_new

        else:
            free_energy_list[0] += 4.92

        return(free_energy_list)
        #__|

    def apply_bias(self, bias, energy_list):
        """Apply bias to free energies.

        Applies a potential to every species in the 4 and 2-electron process
        and adjusts their free energies accordingly
        """
        #| - apply_bias
        mod_free_e_lst = []        # Free energy reaction path at applied bias
        for energy, elec in zip(energy_list, range(len(energy_list))[::-1]):
            mod_free_e_lst.append(energy - elec * bias)

        return(mod_free_e_lst)
        #__|

    def calc_overpotential(self):
        """
        Calculate overpotential for 4e- process.

        Returns the limiting overpotential for the given species and the
        limiting reaction step in the form of a list, species_A -> species_B is
        [species_A, species_B]
        """
        #| - calc_overpotential
        if self.overpotential_given:
            out_list = [None, None]
            if "overpotential" in list(self.fe_df):
                overpot_i = self.fe_df["overpotential"].tolist()[0]
                out_list[0] = overpot_i
            else:
                print("No 'overpotential' column in df")

        else:
            rxn_spec = self.rxn_mech_states

            overpotential_lst = []
            for energy_i in enumerate(self.energy_lst[:-1]):
                energy_i_plus1 = self.energy_lst[energy_i[0] + 1]
                overpotential_i = 1.23 + energy_i_plus1 - energy_i[1]
                overpotential_lst.append(overpotential_i)

            overpotential = max(overpotential_lst)
            # overpotential = min(overpotential_lst)

            lim_step_index = overpotential_lst.index(overpotential)

            limiting_step = [rxn_spec[lim_step_index], rxn_spec[lim_step_index + 1]]
            out_list = [overpotential, limiting_step]

        return(out_list)
        #__|

    def calc_overpotential_OER(self):
        """Calculate the OER overpotential of a ORR series."""
        #| - calc_overpotential_OER
        rxn_spec = self.rxn_mech_states

        overpotential_lst = []
        for energy_i in enumerate(self.energy_lst[:-1]):
            energy_i_plus1 = self.energy_lst[energy_i[0] + 1]
            overpotential_i = energy_i_plus1 - energy_i[1] - 1.23
            overpotential_lst.append(overpotential_i)

        overpotential = max(overpotential_lst)
        lim_step_index = overpotential_lst.index(overpotential)

        limiting_step = [rxn_spec[lim_step_index], rxn_spec[lim_step_index + 1]]
        out_list = [overpotential, limiting_step]

        return(out_list)
        #__|

    def calc_overpotential_h2o2(self):
        """
        Calculate overpotential for 2e- process.

        The overpotential for the 2e- process depends only on the energy of the
        *OOH intermediate
        """
        #| - calc_overpotential_h2o2
        df = self.fe_df
        ooh_row = df[df["adsorbate"] == "ooh"]
        ooh_ads_e = ooh_row.iloc[0]["ads_e"]

        op_4e = ooh_ads_e - 4.22

        return(op_4e)
        #__|

    def __series_plot_name__(self,
        opt_name=None,
        properties=None,
        overpotential_type="ORR",
        add_overpot=True,
        use_key_in_name=False,
        ):
        """Create name for series.

        Args:
            bias:
            opt_name:
            properties:
            color_list:
            hover_text_col:
            plot_mode:
            smart_format:
            overpotential_type:
        """
        #| - __series_plot_name__

        #| - Getting appropriate Overpotential
        if add_overpot:
            if overpotential_type == "ORR":
                overpot_i = self.overpotential
            elif overpotential_type == "OER":
                overpot_i = self.overpotential_OER
            elif overpotential_type == "H2O2":
                overpot_i = self.overpotential_h2o2
            else:
                overpot_i = self.overpotential
        else:
            overpot_i = ""
        #__|

        #| - Connecting properties key: values into string
        if properties is not None:
            properties_string_name = ""
            for key_i, value_i in properties.items():

                if use_key_in_name is True:
                    properties_string_name += str(key_i)
                    properties_string_name += "_"
                properties_string_name += str(value_i)

                properties_string_name += " | "

            # Removing the trailig ' | '
            properties_string_name = properties_string_name[0:-3]
        else:
            properties_string_name = ""
        #__|

        #| - Data Series Name
        if opt_name is not None:
            name_i = opt_name + ": " + properties_string_name

        else:
            name_i = properties_string_name

        if add_overpot:
            name_i += " (OP: " + str(round(overpot_i, 2)) + ")"

        #__|

        # NEW | If name_i given, then just use that
        if self.name_i is not None:
            return(self.name_i)

        return(name_i)
        #__|

    #__|
