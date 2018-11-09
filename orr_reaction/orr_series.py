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
        # system_properties=None,
        state_title="adsorbate",
        free_e_title="ads_e",
        bias=0.,
        rxn_x_coord_array=None,
        group=None,
        name_i=None,
        opt_name=None,

        properties=None,
        property_key_list=None,

        color_list=None,
        color=None,
        i_cnt=0,
        hover_text_col=None,
        plot_mode="all",
        smart_format=None,
        add_overpot=True,
        # overpotential_type="ORR",
        rxn_type="ORR",
        ):
        """
        Input variables to class instance.

        Args:
            free_energy_df:
                Pandas dataframe containing the adsorbates as rows
                Required columns, adsorbate, free energy
            state_title:
            free_e_title:
            bias:
            rxn_x_coord_array:
            opt_name:
                Optional name
                string to append to the beginning of the plot series name
            properties:
            color_list:
            i_cnt:
            hover_text_col:
            plot_mode:
            smart_format:

        """
        #| - __init__

        #| - Setting Instance Attributes
        self.fe_df = free_energy_df
        # self.sys_props = system_properties

        self.group = group
        self.state_title = state_title
        self.fe_title = free_e_title

        self.bias = bias
        self.rxn_x_coord_array = rxn_x_coord_array
        self.name_i = name_i
        self.opt_name = opt_name

        self.properties = properties

        self.color_list = color_list
        self.color = color
        self.i_cnt = i_cnt
        self.hover_text_col = hover_text_col
        self.plot_mode = plot_mode
        self.smart_format = smart_format
        # self.overpotential_type = overpotential_type
        self.rxn_type = rxn_type

        self.add_overpot = add_overpot
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

            # Put this outside of the if-statement
            self.series_name = self.__series_plot_name__(
                opt_name=self.opt_name,
                properties=self.properties,
                overpotential_type=self.rxn_type,
                add_overpot=self.add_overpot,
                )

            # print("__-____*9dfs")
            # print(self.series_name)

            self.series_plot = self.plot_fed_series(
                bias=self.bias,
                opt_name=self.opt_name,
                properties=self.properties,
                color_list=self.color_list,
                i_cnt=self.i_cnt,
                hover_text_col=self.hover_text_col,
                plot_mode=self.plot_mode,
                smart_format=self.smart_format,
                overpotential_type=self.rxn_type,
                )

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

            # Not sure what this was trying to accomplish
            # if len(df_state) == 2:
            #     state_energy_list = []
            #     for j_cnt, row_j in df_state.iterrows():
            #         energy_j = row_j[self.fe_title]
            #         state_energy_list.append(energy_j)

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
            tmp1 = df_state.iloc[0][self.fe_title]
            free_energy_list.append(tmp1)

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
        rxn_spec = self.rxn_mech_states

        overpotential_lst = []
        for energy_i in enumerate(self.energy_lst[:-1]):
            energy_i_plus1 = self.energy_lst[energy_i[0] + 1]
            overpotential_i = 1.23 + energy_i_plus1 - energy_i[1]
            overpotential_lst.append(overpotential_i)

        overpotential = max(overpotential_lst)
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


    #| - Plotting @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    def __convert_to_plotting_list__(self,
        energy_lst,
        spacing=0.5,
        step_size=1,
        ):
        """Repeat entries in energy list to conform to FED plot.

        Modifies an energy list for plotting by repeating each entry
        Ex. [4.92, 3.69, ... ] -> [4.92, 4.92, 3.69, 3.69, ... ]

        Args:
            energy_lst: <type 'list'>
            spacing:
            step_size:
        """
        #| - __convert_to_plotting_list__
        tmp_list = range(len(energy_lst) * 2)
        energy_dupl_lst = [energy_lst[i // 2] for i in tmp_list]

        # rxn_coord_steps = self.create_rxn_coord_array(
        #     len(energy_lst),
        #     spacing=spacing,
        #     step_size=step_size,
        #     )
        # out_list = [rxn_coord_steps, energy_dupl_lst]

        return(energy_dupl_lst)
        #__|

    def __series_plot_name__(self,
        opt_name=None,
        properties=None,
        overpotential_type="ORR",
        add_overpot=True,
        ):
        """Create name for series.

        Args:
            bias:
            opt_name:
            properties:
            color_list:
            i_cnt:
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

                properties_string_name += str(key_i)
                properties_string_name += "_"
                properties_string_name += str(value_i)

                properties_string_name += " | "

            properties_string_name = properties_string_name[0:-3]
        else:
            properties_string_name = ""
        #__|

        #| - Data Series Name
        if opt_name is not None:
            # name_i = opt_name + ": " + properties_string_name + \
            #     " (OP: " + str(round(overpot_i, 2)) + ")"

            name_i = opt_name + ": " + properties_string_name

        else:
            # name_i = properties_string_name + \
            #     " (OP: " + str(round(overpot_i, 2)) + ")"

            name_i = properties_string_name

        if add_overpot:
            name_i += " (OP: " + str(round(overpot_i, 2)) + ")"

        #__|

        # NEW | If name_i given, then just use that
        if self.name_i is not None:
            return(self.name_i)

        return(name_i)
        #__|

    def plot_fed_series(self,
        bias=0.,
        opt_name=None,
        properties=None,
        color_list=None,
        i_cnt=0,
        hover_text_col=None,
        plot_mode="all",
        smart_format=None,
        overpotential_type="ORR",
        ):
        """
        Process data for FED plot.

        Args:
            bias:
            opt_name
            properties:
            color_list:
            i_cnt:
            hover_text_col:
                Dataframe column name to be used for hover text

        #FIXME | This is  fairly rough as of right now
        """
        #| - plot_fed_series


        e_list = self.energy_lst
        e_list = self.apply_bias(bias, e_list)


        for n, i in enumerate(e_list):
            if np.isnan(i) is True:
                e_list[n] = None

        if color_list is None:
            color_list = ["red"]

        name_i = self.series_name

        #| - Hover Text
        if hover_text_col is not None:
            if type(hover_text_col) is not list:
                hover_text_list = self.property_list(hover_text_col)

            else:
                hover_text_lists = []
                for col_i in hover_text_col:

                    # replacing nan with ""
                    tmp = self.property_list(col_i)
                    hover_text_i = ["" if x is np.nan else x for x in tmp]
                    hover_text_lists.append(hover_text_i)

                # TODO Removed last trailing " | "
                hover_text_list = []
                for items in zip(*hover_text_lists):

                    if all([True if i == "" else False for i in items]) is True:
                        hover_col_state_i = ""
                    else:
                        hover_col_state_i = " | ".join(items)

                    hover_text_list.append(hover_col_state_i)

        else:
            hover_text_list = [np.nan for j_cnt in list(range(5))]
        #__|

        #| - TEMP Picking color from "color_list" or "color" variable
        if self.color is not None:
            color_i = self.color
        else:
            color_i = color_list[i_cnt - 1]
        #__|

        dat_lst = self.__create_plotly_series__(
            e_list,
            group=name_i,
            name=name_i,
            hover_text=hover_text_list,
            # color=color_list[i_cnt - 1],
            color=color_i,
            plot_mode=plot_mode,
            smart_format=smart_format,
            )

        return(dat_lst)

        #| - Delete this

        # 181101 Delete this
        # key = properties
        # if type(key) == tuple:
        #     pass
        #
        # elif key is None:
        #     key = None
        # else:
        #     key = (key,)
        # 181101 Delete this
        # if overpotential_type == "ORR":
        #     overpot_i = self.overpotential
        # elif overpotential_type == "OER":
        #     overpot_i = self.overpotential_OER
        # elif overpotential_type == "H2O2":
        #     overpot_i = self.overpotential_h2o2
        # else:
        #     overpot_i = self.overpotential

        # if key is None:
        #     prop_name = ""
        # else:
        #     prop_name = "_".join([str(i) for i in key])
        #
        # if opt_name is not None:
        #     name_i = opt_name + ": " + prop_name + \
        #         " (OP: " + str(round(overpot_i, 2)) + ")"
        #
        # else:
        #     name_i = prop_name + \
        #         " (OP: " + str(round(overpot_i, 2)) + ")"
        #__|

        #__|

    def __create_plotly_series__(self,
        energy_lst,
        name="TEMP",
        group="group1",
        hover_text=None,
        color="rgb(22, 96, 167)",
        plot_mode="all",
        smart_format=None,
        ):
        """
        Create a plotly series for the current instance.

        Args:
            energy_lst:
            name:
            group:
            color:
            plot_mode:
                "all"
                "states_only"
                "full_lines"
        """
        #| - create_plotly_series
        # e_list = self.__convert_to_plotting_list__(energy_lst)
        # x_dat = e_list[0]
        # y_dat = e_list[1]

        # e_list = self.__convert_to_plotting_list__(energy_lst)
        # x_dat = e_list[0]
        y_dat = self.__convert_to_plotting_list__(energy_lst)

        if hover_text is None:
            hover_text = [np.nan for i_ind in range(5)]

        #| - Parameters
        if plot_mode == "all":
            show_leg_2 = False
        elif plot_mode == "states_only":
            show_leg_2 = False
        elif plot_mode == "full_lines":
            show_leg_2 = True
        #__|

        #| - Adding Breaks in Data
        x_dat = self.rxn_x_coord_array

        new_x_dat = copy.copy(x_dat)
        new_y_dat = copy.copy(y_dat)

        cnt = 2
        for i_ind in range(int(len(x_dat) / 2 - 1)):
            fill = new_x_dat[cnt - 1]
            new_x_dat.insert(cnt, fill)
            new_y_dat.insert(cnt, None)
            cnt += 3
        #__|

        #| - Creating x-data in middle of states
        short_y = np.array(y_dat)[::2]

        xdat = list(set(new_x_dat))
        xdat.sort()

        cnt = 0
        short_x = []
        for i_ind in range(int(len(xdat) / 2)):
            short_x.append(xdat[cnt] + 0.5)  # FIXME Replace 0.5 with variable
            cnt += 2
        #__|

        #| - Smart Format Dict ************************************************

        #| - DICTS
        plot_parameter_dict = {
            "dash": None,
            }

        # smart_format = [
        #
        #     # [
        #     #     {"spinpol": True},
        #     #     {"dash": "dot"},
        #     #     ],
        #
        #     [
        #         {"system": "Fe_slab"},
        #         {"dash": "dashdot"},
        #         ],
        #
        #     [
        #         {"system": "N_graph_Fe"},
        #         {"dash": "dot"},
        #         ],
        #
        #     [
        #         {"system": "graph_Fe"},
        #         {"dash": "dash"},
        #         ],
        #
        #     [
        #         {"system": "graphene"},
        #         {"dash": None},
        #         ],
        #
        #     ]

        #__|

        if self.fe_df is not None and smart_format is not None:
            for format_i in smart_format:

                # format_i key:values
                df_col_name = list(format_i[0])[0]
                value_i = list(format_i[0].values())[0]
                setting_name = list(format_i[1])[0]
                setting_value = list(format_i[1].values())[0]


                if df_col_name in list(self.fe_df):

                    # *OOH, *O, *OH entries must have value_i as the
                    # column entries in the column df_col_name

                    df = self.fe_df
                    df_wo_bulk = df[df["adsorbate"] != "bulk"]

                    if all(df_wo_bulk[df_col_name] == value_i):
                        plot_parameter_dict.update(
                            {setting_name: setting_value},
                            )

                else:
                    print("Dataframe column " + df_col_name + " not present")

        #__| ******************************************************************

        #| - Series Color
        if "color" in list(plot_parameter_dict):
            color_out = plot_parameter_dict["color"]
        else:
            plot_parameter_dict["color"] = color
        #__|

        #| - Plotly Scatter Plot Instances

        #| - Thick horizontal state lines
        data_1 = Scatter(
            x=new_x_dat,
            y=new_y_dat,
            legendgroup=group,
            showlegend=True,
            name=name,
            hoverinfo="none",  # TEMP - 180317
            # text=hover_text,

            connectgaps=False,
            line=dict(
                # color=color,
                color=plot_parameter_dict["color"],
                width=2,
                # dash="dot",  # TEMP
                dash=plot_parameter_dict["dash"],  # TEMP
                ),
            mode="lines",
            )
        #__|

        #| - Full, thin line
        data_2 = Scatter(
            x=new_x_dat,
            y=new_y_dat,
            legendgroup=group,
            name=name,
            connectgaps=True,
            showlegend=show_leg_2,
            hoverinfo="none",
            text=hover_text,

            line=dict(
                # color=color,
                color=plot_parameter_dict["color"],
                width=1,
                ),
            mode="lines",
            )
        #__|

        #| - Points in middle of energy states
        data_3 = Scatter(
            x=short_x,
            y=short_y,
            legendgroup=group,
            name=name,
            showlegend=False,
            hoverinfo="y+text",
            text=hover_text,
            marker=dict(
                size=14,
                color=color,
                opacity=0.,
                ),
            mode="markers",
            )
        #__|

        #__|

        #| - Plot Mode (which data series to plot)
        if plot_mode == "all":
            data_lst = [data_1, data_2, data_3]
        elif plot_mode == "states_only":
            data_lst = [data_1, data_3]
        elif plot_mode == "full_lines":
            data_lst = [data_2, data_3]
        #__|

        return(data_lst)
        #__|


    #__| @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    #__|




#| - __old__

    #| - __old__
        # bias=0.,
        # opt_name=None,
        # properties=None,
        # color_list=None,
        # i_cnt=0,
        # hover_text_col=None,
        # plot_mode="all",
        # smart_format=None,
    #__|

    #| - __old__
    # def create_rxn_coord_array(self,
    #     rxn_steps,
    #     spacing=0,
    #     step_size=1,
    #     ):
    #     """
    #     Create a reaction coordinate array ([0, 1, 1, 2, 2, 3]) for plotting.
    #
    #     Args:
    #         rxn_steps: <type 'int'>
    #             Number of steps in reaction coordinate including initial
    #             and final state.
    #             Ex. A -> Intermediate -> C has 3 steps/states
    #
    #         spacing: <type 'float'>
    #             Spacing inbetween the energy levels. The default of 0 creates
    #             a free energy diagram that looks like steps
    #     """
    #     #| - create_rxn_coord_array
    #     lst = []
    #     for i in range(1, rxn_steps):
    #         if i == 1:
    #             lst.append(step_size)
    #             lst.append(step_size + spacing)
    #         if i != 1:
    #             lst.append(lst[-1] + step_size)
    #             lst.append(lst[-2] + step_size + spacing)
    #
    #     lst.insert(0, 0)
    #     lst.append(lst[-1] + step_size)
    #
    #     return(lst)
    #     #__|
    #__|

#__|
